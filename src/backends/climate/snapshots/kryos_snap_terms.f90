module kryos_snap_terms
    ! The terms of the climate forcing (yelmox climate_terms, snapshots
    ! branch): each term is one coefficient times one pattern, added to the
    ! reference climate on an explicit list of fields:
    !   snapshots -- the difference of two snapshots, snap1 - snap0 (or snap1
    !                itself without snap0), applied by the snapshot code at sea
    !                level (kryos_snap_snapshots reads these terms);
    !   uniform   -- a constant pattern: tas + f_ta*f (and pr * exp(beta_p*
    !                f_ta*f)), to + f_to*f, so + f_so*f, applied here;
    !   replace   -- the field takes the value of the coefficient at every depth
    !                in place of the reference (to, so); the other terms on the
    !                field add to it.
    ! The terms are the names of the instance key `terms`, each a group
    ! [term_<name>] from the par file, else the first library file that has
    ! it; unset keys from the template `term` of input/snapshots_defaults.nml.
    ! Their coefficients go into the adapter's private table
    ! (kryos_snap_coefs).
    !
    ! The value f of each term's coefficient is evaluated here once per
    ! update (terms_eval), kept in the term, read by the snapshot code and the
    ! applications below, and written for every term in one place
    ! (terms_write_step).
    !
    ! Order within a climate update (the adapter): terms_eval -> snapshots
    ! (with the snapshot terms) -> the uniform terms on the atmosphere
    ! (terms_apply_atm) -> the derived coefficients (set by the adapter from
    ! the composed atmosphere; terms_eval of the terms they drive) -> the ocean
    ! terms (terms_apply_ocn) -> taper -> bias. The ocean application is
    ! idempotent: it rebuilds the ocean from a copy of the snapshots' ocean
    ! (terms_save_ocean).
    !
    ! Contracts, checked at init (terms_check) against what the backend
    ! supplies (climate_caps_class), an error naming the term:
    !   - a uniform term on tas or pr needs an atmosphere that drives the surface;
    !   - a derived coefficient is computed after the atmosphere, so its term
    !     may touch to and so only;
    !   - a uniform term on to or so needs ocean profiles; a replace term: to
    !     or so only, one per field, ocean profiles.

    use precision,        only: wp
    use fesm_log,         only: LOG_INFO, log_msg, fesm_error
    use nml,              only: nml_read, nml_validate
    use ncio,             only: nc_write
    use kryos_snap_out,   only: climate_out_class
    use kryos_snap_coefs, only: snap_coef_table_t, coefs_add, coefs_value, coefs_is_derived
    use kryos_snap_nml,   only: SNAP_DEFAULTS, FILE_LEN, group_file

    implicit none

    private

    integer, parameter :: FLD_MAX = 4    ! tas, pr, to, so

    ! What the backend supplies, for the contracts (set by the adapter).
    type, public :: climate_caps_class
        character(len=16) :: name = ""           ! the backend, for messages
        logical :: atm_drives_srf = .false.      ! the surface is computed from the atmosphere of climate_out
        logical :: has_tsl_ann    = .false.      ! sea-level air temperature (now and ref)
        logical :: has_snapshots  = .false.      ! reads the snapshot terms
        logical :: has_ocn_prof   = .false.      ! ocean profiles (to_ann, so_ann)
    end type climate_caps_class

    ! One term: coefficient x pattern on a list of fields.
    type, public :: term_class
        character(len=64) :: name    = ""
        character(len=64) :: group   = ""
        character(len=64) :: coef    = ""          ! the coefficient; "" = the constant 1
        character(len=16) :: pattern = "snapshots" ! snapshots | uniform | replace
        character(len=64) :: snap0   = ""          ! snapshots: the snapshot at coefficient 0 ("" = none)
        character(len=64) :: snap1   = ""          ! snapshots: the snapshot at coefficient 1
        character(len=8)  :: fields(FLD_MAX) = ""
        logical  :: on_tas = .false., on_pr = .false., on_to = .false., on_so = .false.
        real(wp) :: f_ta = 1.0_wp                  ! [K]   uniform: tas + f_ta*f (pr * exp(beta_p*f_ta*f))
        real(wp) :: f_to = 1.0_wp                  ! [K]   uniform: to + f_to*f
        real(wp) :: f_so = 1.0_wp                  ! [psu] uniform: so + f_so*f
        logical  :: derived = .false.              ! the coefficient is derived (set by the model)
        real(wp) :: f = 0.0_wp                     ! the coefficient value (1 without one), terms_eval
    end type term_class

    type, public :: terms_class
        character(len=64) :: who = ""               ! the instance, for messages
        type(term_class), allocatable :: term(:)
        integer :: n = 0
        ! The snapshots' ocean of the last update (terms_apply_ocn is idempotent).
        logical :: ocn_saved = .false.
        real(wp), allocatable :: to_base(:,:,:), so_base(:,:,:)
    end type terms_class

    public :: caps_need_atm
    public :: terms_init, terms_check, terms_eval, terms_save_ocean, terms_apply_atm, terms_apply_ocn
    public :: terms_write_step, coef_label

contains

    subroutine terms_init(terms, par_file, library, names, coefs, who)
        ! Read each [term_<name>] of `names` (blank entries ignored) and add
        ! its coefficient to the table; the backend-dependent contracts follow
        ! in terms_check once the backend is initialised.
        type(terms_class),       intent(inout) :: terms
        character(len=*),        intent(in)    :: par_file, library(:), names(:), who
        type(snap_coef_table_t), intent(inout) :: coefs

        character(len=FILE_LEN) :: gfile
        integer :: k, j

        terms%who = who
        terms%n   = count(len_trim(names) > 0)
        if (.not. allocated(coefs%coef)) allocate(coefs%coef(0))    ! no term: an empty table (size() of an unallocated one is undefined)
        if (allocated(terms%term)) deallocate(terms%term)
        allocate(terms%term(terms%n))

        do k = 1, terms%n
            associate (t => terms%term(k))
            t%name = trim(names(k))
            if (any(names(1:k-1) == names(k))) call term_error(terms, t, "listed twice in the terms of "//trim(who)//".")
            t%group = "term_"//trim(t%name)
            gfile = group_file(par_file, library, trim(t%group), "a term of "//trim(who))
            call nml_validate(trim(gfile), SNAP_DEFAULTS, trim(t%group), defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "coef",    t%coef,    init=.true., &
                          defaults_file=SNAP_DEFAULTS, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "pattern", t%pattern, init=.true., &
                          defaults_file=SNAP_DEFAULTS, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "snap0",   t%snap0,   init=.true., &
                          defaults_file=SNAP_DEFAULTS, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "snap1",   t%snap1,   init=.true., &
                          defaults_file=SNAP_DEFAULTS, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "fields",  t%fields,  init=.true., &
                          defaults_file=SNAP_DEFAULTS, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "f_ta",    t%f_ta, &
                          defaults_file=SNAP_DEFAULTS, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "f_to",    t%f_to, &
                          defaults_file=SNAP_DEFAULTS, defaults_group="term")
            call nml_read(trim(gfile), trim(t%group), "f_so",    t%f_so, &
                          defaults_file=SNAP_DEFAULTS, defaults_group="term")

            select case (trim(t%pattern))
            case ("snapshots", "uniform", "replace")
            case default
                call term_error(terms, t, "pattern = '"//trim(t%pattern)//"'; must be snapshots, uniform or replace.")
            end select

            if (count(len_trim(t%fields) > 0) == 0) call term_error(terms, t, "fields is empty.")
            do j = 1, FLD_MAX
                select case (trim(t%fields(j)))
                case ("")
                case ("tas"); t%on_tas = .true.
                case ("pr");  t%on_pr  = .true.
                case ("to");  t%on_to  = .true.
                case ("so");  t%on_so  = .true.
                case default
                    call term_error(terms, t, "field '"//trim(t%fields(j))//"'; must be tas, pr, to or so.")
                end select
            end do

            if (len_trim(t%coef) > 0) then
                call coefs_add(coefs, par_file, library, trim(t%coef))
                t%derived = coefs_is_derived(coefs, t%coef)
            end if

            if (trim(t%pattern) == "snapshots") then
                if (len_trim(t%snap1) == 0) call term_error(terms, t, "pattern snapshots needs snap1.")
            else
                if (len_trim(t%snap0) > 0 .or. len_trim(t%snap1) > 0) &
                    call term_error(terms, t, "snap0/snap1 belong to the snapshots pattern.")
            end if
            if (trim(t%pattern) == "replace" .and. (t%on_tas .or. t%on_pr)) &
                call term_error(terms, t, "pattern replace is for the ocean (to, so) only.")
            if (t%derived .and. (t%on_tas .or. t%on_pr)) &
                call term_error(terms, t, "its coefficient "//trim(t%coef)//" is derived (computed after "// &
                                "the atmosphere), so it can only drive to or so.")
            if (len_trim(t%coef) == 0 .and. trim(t%pattern) /= "snapshots") &
                call term_error(terms, t, "a "//trim(t%pattern)//" term needs a coefficient.")

            call log_msg(LOG_INFO, "kryos_snap_terms", trim(t%group)//": "//trim(t%pattern)// &
                         ", coefficient "//trim(coef_label(t%coef))//", fields "//trim(fields_list(t)), &
                         "instance "//trim(who))
            end associate
        end do

        do k = 1, terms%n
            associate (t => terms%term(k))
            if (trim(t%pattern) /= "replace") cycle
            do j = 1, k-1
                if (trim(terms%term(j)%pattern) /= "replace") cycle
                if ((t%on_to .and. terms%term(j)%on_to) .or. (t%on_so .and. terms%term(j)%on_so)) &
                    call term_error(terms, t, "a second replace term on the same field as "// &
                                    trim(terms%term(j)%group)//".")
            end do
            end associate
        end do

    end subroutine terms_init

    subroutine terms_check(terms, caps)
        ! The contracts that depend on the backend (see the header).
        type(terms_class),        intent(in) :: terms
        type(climate_caps_class), intent(in) :: caps
        integer :: k

        do k = 1, terms%n
            associate (t => terms%term(k))
            select case (trim(t%pattern))
            case ("snapshots")
                if (.not. caps%has_snapshots) &
                    call term_error(terms, t, "pattern snapshots needs the snapshot backend; climate = "// &
                                    trim(caps%name)//" reads no snapshots.")
            case ("uniform")
                if (t%on_tas .or. t%on_pr) call caps_need_atm(caps, "["//trim(t%group)//"] (a uniform term on tas/pr)")
                if ((t%on_to .or. t%on_so) .and. .not. caps%has_ocn_prof) &
                    call term_error(terms, t, "lists to or so, but climate = "//trim(caps%name)// &
                                    " has no ocean profiles.")
            case ("replace")
                if (.not. caps%has_ocn_prof) &
                    call term_error(terms, t, "pattern replace needs ocean profiles; climate = "// &
                                    trim(caps%name)//" has none.")
            end select
            end associate
        end do

    end subroutine terms_check

    subroutine caps_need_atm(caps, who)
        ! The one predicate of the terms and transforms that change the
        ! atmosphere (tas, pr) after the backend: stop unless it drives the
        ! surface (caps%atm_drives_srf). `who` names the term or transform.
        type(climate_caps_class), intent(in) :: caps
        character(len=*),         intent(in) :: who

        if (caps%atm_drives_srf) return
        call fesm_error("kryos_snap_terms::caps_need_atm", trim(who)//" changes the atmosphere after the "// &
                        "backend, which cannot reach the surface", "climate = "//trim(caps%name))
    end subroutine caps_need_atm

    ! ===== application ======================================================

    subroutine terms_eval(terms, coefs, derived)
        ! The value f of each term's coefficient from the table (1 for a
        ! term without one): every term, once per update before the backend.
        ! With derived = .true., after the derived coefficients are set: only
        ! the uniform and replace terms with a derived coefficient, applied
        ! after that (terms_apply_ocn); a snapshot term keeps the value the
        ! backend applied.
        type(terms_class),       intent(inout) :: terms
        type(snap_coef_table_t), intent(in)    :: coefs
        logical,                 intent(in), optional :: derived

        logical :: only_derived
        integer :: k

        only_derived = .false.
        if (present(derived)) only_derived = derived

        do k = 1, terms%n
            associate (t => terms%term(k))
            if (only_derived .and. (.not. t%derived .or. trim(t%pattern) == "snapshots")) cycle
            t%f = 1.0_wp
            if (len_trim(t%coef) > 0) t%f = coefs_value(coefs, t%coef)
            end associate
        end do

    end subroutine terms_eval

    subroutine terms_apply_atm(terms, out)
        ! The uniform terms on the atmosphere: tas, ta_ann, tsl_ann + f_ta*f;
        ! pr, pr_ann * exp(beta_p*f_ta*f) (beta_p: the backend's sensitivity).
        type(terms_class),       intent(in)    :: terms
        type(climate_out_class), intent(inout) :: out

        real(wp) :: dTa
        integer  :: k, m

        do k = 1, terms%n
            associate (t => terms%term(k), now => out%now)
            if (trim(t%pattern) /= "uniform") cycle
            if (.not. (t%on_tas .or. t%on_pr)) cycle
            dTa = t%f_ta*t%f
            if (t%on_tas) then
                now%tas    = now%tas + dTa
                now%ta_ann = now%ta_ann + dTa
                if (allocated(now%tsl_ann)) now%tsl_ann = now%tsl_ann + dTa
            end if
            if (t%on_pr) then
                do m = 1, size(now%pr, 3)
                    now%pr(:,:,m) = now%pr(:,:,m) * exp(now%beta_p*dTa)
                end do
                if (allocated(now%pr_ann)) now%pr_ann = now%pr_ann * exp(now%beta_p*dTa)
            end if
            end associate
        end do

    end subroutine terms_apply_atm

    subroutine terms_save_ocean(terms, out)
        ! Keep the snapshots' ocean of this update, the base the ocean terms
        ! are applied to.
        type(terms_class),       intent(inout) :: terms
        type(climate_out_class), intent(in)    :: out

        terms%ocn_saved = .true.
        if (allocated(out%now%to_ann)) terms%to_base = out%now%to_ann
        if (allocated(out%now%so_ann)) terms%so_base = out%now%so_ann

    end subroutine terms_save_ocean

    subroutine terms_apply_ocn(terms, out)
        ! The ocean terms on the saved base: a replace term sets the field to
        ! its coefficient at every depth, the uniform ones add f_to*f / f_so*f
        ! (the climate has the ocean a term needs: terms_check). Idempotent.
        type(terms_class),       intent(in)    :: terms
        type(climate_out_class), intent(inout) :: out

        logical  :: any_ocn
        integer  :: k

        any_ocn = .false.
        do k = 1, terms%n
            if (trim(terms%term(k)%pattern) == "snapshots") cycle
            if (terms%term(k)%on_to .or. terms%term(k)%on_so) any_ocn = .true.
        end do
        if (.not. any_ocn) return

        if (.not. terms%ocn_saved) then
            call fesm_error("kryos_snap_terms::terms_apply_ocn", "no saved ocean (terms_save_ocean first)", &
                            "instance "//trim(terms%who))
        end if

        ! The base.
        if (allocated(terms%to_base)) out%now%to_ann = terms%to_base
        if (allocated(terms%so_base)) out%now%so_ann = terms%so_base

        ! Replace terms first.
        do k = 1, terms%n
            associate (t => terms%term(k))
            if (trim(t%pattern) /= "replace") cycle
            if (t%on_to) out%now%to_ann = t%f
            if (t%on_so) out%now%so_ann = t%f
            end associate
        end do

        ! Then the uniform ones.
        do k = 1, terms%n
            associate (t => terms%term(k))
            if (trim(t%pattern) /= "uniform") cycle
            if (.not. (t%on_to .or. t%on_so)) cycle
            if (t%on_to) out%now%to_ann = out%now%to_ann + t%f_to*t%f
            if (t%on_so) out%now%so_ann = out%now%so_ann + t%f_so*t%f
            end associate
        end do

    end subroutine terms_apply_ocn

    subroutine terms_write_step(terms, filename, ncid, n)
        ! Record n of every term into an open 2D output file: the value of
        ! its coefficient, term_<name> (1 for a snapshot term without one).
        type(terms_class), intent(in) :: terms
        character(len=*),  intent(in) :: filename
        integer,           intent(in) :: ncid, n
        integer :: k
        do k = 1, terms%n
            associate (t => terms%term(k))
            call nc_write(filename, "term_"//trim(t%name), t%f, dim1="time", &
                          start=[n], count=[1], ncid=ncid, units="1", &
                          long_name="Coefficient "//trim(coef_label(t%coef))//" of the "// &
                                    trim(t%pattern)//" term "//trim(t%name))
            end associate
        end do
    end subroutine terms_write_step

    ! ===== helpers ==========================================================

    function fields_list(t) result(s)
        type(term_class), intent(in) :: t
        character(len=48) :: s
        integer :: j
        s = ""
        do j = 1, FLD_MAX
            if (len_trim(t%fields(j)) > 0) s = trim(s)//" "//trim(t%fields(j))
        end do
        s = adjustl(s)
    end function fields_list

    function coef_label(name) result(s)
        ! The coefficient name, or "1" for none.
        character(len=*), intent(in) :: name
        character(len=64) :: s
        s = "1"
        if (len_trim(name) > 0) s = name
    end function coef_label

    subroutine term_error(terms, t, msg)
        type(terms_class), intent(in) :: terms
        type(term_class),  intent(in) :: t
        character(len=*),  intent(in) :: msg
        call fesm_error("kryos_snap_terms", "["//trim(t%group)//"] "//trim(msg), "instance "//trim(terms%who))
    end subroutine term_error

end module kryos_snap_terms
