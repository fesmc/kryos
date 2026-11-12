module kryos_ice_data
    ! ice backend data: the ice geometry read, not solved. The thickness, and optionally the bed and the
    ! sea surface, come from field groups (fesm_datasrc; instance keys H_ice, z_bed, z_sl), static or
    ! transient (an ice history such as PaleoMIST, a reference topography): every component downstream
    ! (the hub, a climate, a surface, the global ice field of the globe) sees them as an ice model's
    ! exports. Without H_ice the ice and without z_bed the bed are the domain's reference topography (a
    ! present-day run on it), without z_sl the sea surface 0.
    !
    ! Each domain has its own instance and its own sources; the hub of each domain takes its ice, and
    ! globe_ice merges every domain's into the background of the globe (a global climate, a global geo).
    ! The data are read at the instance's evaluation time (the forcing time with clock = "forcing").
    ! Stateless: no restart file (the fields are read again at the restored time).
    !
    ! Exports as plastic's: H_ice, z_srf (z_bed + H grounded, z_sl + (1 - rho_i/rho_w) H floating, the bed
    ! or the sea surface without ice), z_bed, z_sl, f_ice, f_grnd (the flotation rule with the run's
    ! densities), H_grnd, dHidt and mb_net (the change of the data over the step), smb = bmb = cmb = 0;
    ! static z_bed_ref, H_ice_pd (the reference topography), footprint (the whole domain).

    use precision,       only: dp
    use fesm_log,        only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields,     only: field_state_t, fields_init, fields_add, fields_lock, fields_put, &
                               FIELD_INSTANT, FIELD_INTENSIVE, FIELD_EXTENSIVE
    use fesm_datasrc,    only: datasrc_t, datasrc_init, datasrc_update, datasrc_get, datasrc_describe, datasrc_end
    use nml,             only: nml_read, nml_has_param
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t, component_optional, NAME_LEN, &
                               PASS_INIT, PASS_GEOMETRY, PASS_SPINUP, PASS_RESTART
    use kryos_ice,       only: ice_t, ice_opt_t

    implicit none

    private
    public :: data_ice_t, data_ice_constraints

    character(len=*), parameter :: DATASRC_DEFAULTS = "input/datasrc_defaults.nml"
    real(dp),         parameter :: KG_TO_GT = 1.0e-12_dp

    type, extends(ice_t) :: data_ice_t
        type(init_ctx_t)    :: ctx
        type(datasrc_t)     :: ds
        type(field_state_t) :: st                        ! the fields the datasrc serves (H_ice, z_bed, z_sl)
        logical  :: has_H = .false., has_zb = .false., has_zsl = .false., has_ds = .false.   ! has_ds: any bound
        real(dp) :: rho_i = 0.0_dp, rho_w = 0.0_dp, dx = 0.0_dp, dy = 0.0_dp
        real(dp), allocatable :: H(:,:), z_bed(:,:), z_sl(:,:), H_prev(:,:)
        real(dp), allocatable :: z_bed_ref(:,:), H_ice_pd(:,:), footprint(:,:)
    contains
        procedure :: init          => data_init_
        procedure :: init_state    => data_init_state_
        procedure :: update        => data_update_
        procedure :: restart_write => data_restart_write_
        procedure :: restart_read  => data_restart_read_
        procedure :: flotation     => data_flotation_
        procedure :: end           => data_end_
    end type data_ice_t

contains

    function data_ice_constraints() result(con)
        ! Regional, any grid, several instances; no backend restart file.
        type(constraints_t) :: con
        con = constraints_t(restart=.false.)
    end function data_ice_constraints

    subroutine data_init_(self, ctx)
        ! The field groups of the instance (H_ice, z_bed, z_sl; each optional), on the instance grid and the
        ! run's calendar; the reference topography; the run's densities.
        class(data_ice_t), intent(inout) :: self
        type(init_ctx_t),  intent(in)    :: ctx

        character(len=8) :: names(3)
        integer :: h, d2(2)

        self%ctx  = ctx
        self%prec = dp
        self%caps%has_flotation = .true.
        self%nx = ctx%grid%G%nx
        self%ny = ctx%grid%G%ny
        self%dx = ctx%grid%G%dx*ctx%grid%cs%xy_conv
        self%dy = ctx%grid%G%dy*ctx%grid%cs%xy_conv
        self%rho_i = ctx%cnst%rho_ice
        self%rho_w = ctx%cnst%rho_sw
        self%opt = ice_opt_t(z_bed_ref=.true., H_ice_ref=.true., footprint=.true., H_ice_pd=.true., H_grnd=.true.)
        call component_optional(self, "z_bed")
        call component_optional(self, "z_sl")
        call component_optional(self, "smb")
        call component_optional(self, "T_srf")
        call component_optional(self, "bmb_shlf")

        if (.not. (allocated(ctx%ref%z_bed) .and. allocated(ctx%ref%H_ice))) then
            call fesm_error("kryos_ice_data::init", "reference topography incomplete", "instance "// &
                            trim(self%name)//": z_bed and H_ice are needed on grid `"//trim(self%grid)//"`")
        end if
        self%z_bed_ref = ctx%ref%z_bed
        self%H_ice_pd  = ctx%ref%H_ice
        allocate(self%footprint(self%nx, self%ny), self%H(self%nx, self%ny), self%z_bed(self%nx, self%ny), &
                 self%z_sl(self%nx, self%ny), self%H_prev(self%nx, self%ny))
        self%footprint = 1.0_dp

        self%has_H   = bound(ctx, "H_ice")
        self%has_zb  = bound(ctx, "z_bed")
        self%has_zsl = bound(ctx, "z_sl")

        d2 = [self%nx, self%ny]
        call fields_init(self%st, trim(self%name)//":data")
        call fields_add(self%st, "H_ice", "m", FIELD_INSTANT, FIELD_EXTENSIVE, h, grid=trim(ctx%grid_name), dims=d2)
        call fields_add(self%st, "z_bed", "m", FIELD_INSTANT, FIELD_INTENSIVE, h, grid=trim(ctx%grid_name), dims=d2)
        call fields_add(self%st, "z_sl",  "m", FIELD_INSTANT, FIELD_INTENSIVE, h, grid=trim(ctx%grid_name), dims=d2)
        call fields_lock(self%st)
        names = [character(len=8) :: "H_ice", "z_bed", "z_sl"]
        self%has_ds = self%has_H .or. self%has_zb .or. self%has_zsl
        if (self%has_ds) then
            call datasrc_init(self%ds, trim(ctx%par_file), trim(ctx%instance), self%st, names, ctx%grid, &
                              DATASRC_DEFAULTS, required=pack(names, [self%has_H, self%has_zb, self%has_zsl]), &
                              var_defs=pack(ctx%library, len_trim(ctx%library) > 0), domain=trim(ctx%domain), &
                              time_ref=ctx%time_ref)
            call datasrc_describe(self%ds)
        end if
        call log_msg(LOG_INFO, "kryos_ice_data", "the ice geometry: H_ice "// &
                     merge("from data               ", "the reference topography", self%has_H)//"; z_bed "// &
                     merge("from data               ", "the reference topography", self%has_zb)// &
                     "; z_sl "//merge("from data", "0        ", self%has_zsl)//"; rho_i "//to_str(self%rho_i)// &
                     ", rho_w "//to_str(self%rho_w), "instance "//trim(self%name))
    end subroutine data_init_

    logical function bound(ctx, key)
        ! Does the instance group bind the field `key` to a field group?
        type(init_ctx_t), intent(in) :: ctx
        character(len=*), intent(in) :: key
        character(len=NAME_LEN) :: grp
        grp = ""
        if (nml_has_param(ctx%par_file, ctx%instance, key)) call nml_read(ctx%par_file, ctx%instance, key, grp)
        bound = len_trim(grp) > 0
    end function bound

    subroutine read_at(self, t)
        ! The fields at time t (model time; the datasrc maps it to the files' calendar).
        class(data_ice_t), intent(inout) :: self
        real(dp),          intent(in)    :: t

        if (self%has_ds) call datasrc_update(self%ds, t)
        if (self%has_H) then
            call datasrc_get(self%ds, "H_ice", self%H)
            self%H = max(self%H, 0.0_dp)
        else
            self%H = self%H_ice_pd
        end if
        if (self%has_zb) then
            call datasrc_get(self%ds, "z_bed", self%z_bed)
        else
            self%z_bed = self%z_bed_ref
        end if
        if (self%has_zsl) then
            call datasrc_get(self%ds, "z_sl", self%z_sl)
        else
            self%z_sl = 0.0_dp
        end if
    end subroutine read_at

    subroutine data_init_state_(self, st, pass)
        ! PASS_INIT: the data at t1 and every INSTANT export, the static ones included. PASS_GEOMETRY: the
        ! same again. PASS_SPINUP: the MEAN exports of the virtual previous step (no change). PASS_RESTART:
        ! the data at the restored time (the core restored the exports).
        class(data_ice_t), intent(inout) :: self
        type(step_time_t), intent(in)    :: st
        integer,           intent(in)    :: pass

        select case (pass)
        case (PASS_INIT)
            call read_at(self, st%t_eval)
            self%H_prev = self%H
            call fields_put(self%exports, self%h_z_bed_ref, self%z_bed_ref, t1=st%t1)
            call fields_put(self%exports, self%h_footprint, self%footprint, t1=st%t1)
            call fields_put(self%exports, self%h_H_ice_pd,  self%H_ice_pd,  t1=st%t1)
            call exports_geometry(self, st%t1)
        case (PASS_GEOMETRY)
            call read_at(self, st%t_eval)
            self%H_prev = self%H
            call exports_geometry(self, st%t1)
        case (PASS_SPINUP)
            call exports_rates(self, st)
        case (PASS_RESTART)
            call read_at(self, st%t_eval)
            self%H_prev = self%H
        end select
    end subroutine data_init_state_

    subroutine data_update_(self, st)
        ! The data at t1 (the evaluation time); the rates over [t0, t1].
        class(data_ice_t), intent(inout) :: self
        type(step_time_t), intent(in)    :: st

        self%H_prev = self%H
        call read_at(self, st%t_eval)
        call exports_geometry(self, st%t1)
        call exports_rates(self, st)
    end subroutine data_update_

    subroutine data_flotation_(self, H, z_bed, z_sl, f_grnd)
        ! Grounded where rho_i H >= rho_w (z_sl - z_bed), with the run's densities.
        class(data_ice_t), intent(in)  :: self
        real(dp),          intent(in)  :: H(:,:), z_bed(:,:), z_sl(:,:)
        real(dp),          intent(out) :: f_grnd(:,:)
        f_grnd = merge(1.0_dp, 0.0_dp, self%rho_i*H >= self%rho_w*(z_sl - z_bed))
    end subroutine data_flotation_

    subroutine exports_geometry(self, t1)
        ! INSTANT: H_ice, z_srf, z_bed, z_sl, f_ice, f_grnd, H_grnd, H_ice_ref, M_ice.
        class(data_ice_t), intent(inout) :: self
        real(dp),          intent(in)    :: t1

        real(dp), allocatable :: f_grnd(:,:), zs(:,:)

        allocate(f_grnd(self%nx, self%ny), zs(self%nx, self%ny))
        call data_flotation_(self, self%H, self%z_bed, self%z_sl, f_grnd)
        where (self%H > 0.0_dp .and. f_grnd > 0.5_dp)
            zs = self%z_bed + self%H
        elsewhere (self%H > 0.0_dp)
            zs = self%z_sl + (1.0_dp - self%rho_i/self%rho_w)*self%H
        elsewhere
            zs = max(self%z_bed, self%z_sl)
        end where
        call fields_put(self%exports, self%h_H_ice,  self%H,     t1=t1)
        call fields_put(self%exports, self%h_z_srf,  zs,         t1=t1)
        call fields_put(self%exports, self%h_z_bed,  self%z_bed, t1=t1)
        call fields_put(self%exports, self%h_z_sl,   self%z_sl,  t1=t1)
        call fields_put(self%exports, self%h_f_ice,  merge(1.0_dp, 0.0_dp, self%H > 0.0_dp), t1=t1)
        call fields_put(self%exports, self%h_f_grnd, f_grnd,     t1=t1)
        call fields_put(self%exports, self%h_H_grnd, self%H - (self%rho_w/self%rho_i)*(self%z_sl - self%z_bed), t1=t1)
        call fields_put(self%exports, self%h_H_ice_ref, self%H_ice_pd, t1=t1)
        call fields_put(self%exports, self%h_M_ice, self%rho_i*sum(self%H)*self%dx*self%dy*KG_TO_GT, t1=t1)
    end subroutine exports_geometry

    subroutine exports_rates(self, st)
        ! MEAN over [t0, t1]: dHidt, mb_net = rho_i dHidt; smb = bmb = cmb = 0 (the ice is imposed).
        class(data_ice_t), intent(inout) :: self
        type(step_time_t), intent(in)    :: st

        real(dp), allocatable :: dHdt(:,:)

        allocate(dHdt(self%nx, self%ny))
        if (st%t1 > st%t0) then
            dHdt = (self%H - self%H_prev)/(st%t1 - st%t0)
        else
            dHdt = 0.0_dp
        end if
        call fields_put(self%exports, self%h_dHidt,  dHdt,                t1=st%t1, t0=st%t0)
        call fields_put(self%exports, self%h_smb,    0.0_dp*dHdt,         t1=st%t1, t0=st%t0)
        call fields_put(self%exports, self%h_bmb,    0.0_dp*dHdt,         t1=st%t1, t0=st%t0)
        call fields_put(self%exports, self%h_cmb,    0.0_dp*dHdt,         t1=st%t1, t0=st%t0)
        call fields_put(self%exports, self%h_mb_net, self%rho_i*dHdt,     t1=st%t1, t0=st%t0)
    end subroutine exports_rates

    subroutine data_restart_write_(self, dir, t)
        ! Stateless: constraints restart = .false., so the core never calls this.
        class(data_ice_t), intent(inout) :: self
        character(len=*),  intent(in)    :: dir
        real(dp),          intent(in)    :: t
        associate (unused => t + len(dir) + self%nx); end associate
    end subroutine data_restart_write_

    subroutine data_restart_read_(self, dir, t)
        class(data_ice_t), intent(inout) :: self
        character(len=*),  intent(in)    :: dir
        real(dp),          intent(in)    :: t
        associate (unused => t + len(dir) + self%nx); end associate
    end subroutine data_restart_read_

    subroutine data_end_(self)
        class(data_ice_t), intent(inout) :: self
        if (self%has_ds) call datasrc_end(self%ds)
    end subroutine data_end_

end module kryos_ice_data
