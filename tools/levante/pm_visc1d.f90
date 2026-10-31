program pm_visc1d
   !! PaleoMIST test of the standalone VILMA2: the 1-D radial viscosity profile that
   !! pre_spinup_1d uses for a 3-D run (solid_earth_spinup in vilma_coupling.f90), from
   !! VILMA2's own routines. With the run's namelist: the Gauss grid of vilma.x (build_grid
   !! of vilma_drive), the FE radial nodes of the response (response_init_ve), the 3-D
   !! log10(eta) on grid x nodes (load_visc_3d: read, perturbation, clamp), then, per node,
   !! the unweighted mean over the grid points, as solid_earth_spinup does.
   !!   cd <3-D run folder>; pm_visc1d.x pm.nml <out.csv>
   !! Writes r_node [m], depth [m], log10 eta (the mean), and the min and max over the grid.
   !! Built against libvilma.a of the stack (pm_visc1d.sh).
   use vilma_precision,       only: wp
   use vilma_constants,       only: kyr
   use vilma_params,          only: vilma_param_class, vilma_par_load
   use vilma_control,         only: DEFAULTS_FILE
   use vilma_earth_structure, only: earth_model, build_earth, load_visc_3d
   use vilma_response,        only: response, response_init_ve
   use vilma_sht,             only: sht_grid, sht_grid_init
   implicit none

   type(vilma_param_class) :: p
   type(sht_grid)          :: sht
   type(earth_model)       :: em
   type(response)          :: resp
   real(wp), allocatable   :: visc_node(:,:)
   character(len=512)      :: cfg, out
   integer :: nlat, nphi, nh, r, u

   call get_command_argument(1, cfg);  call get_command_argument(2, out)
   if (len_trim(cfg) == 0 .or. len_trim(out) == 0) error stop 'usage: pm_visc1d.x <nml> <out.csv>'
   call vilma_par_load(p, trim(cfg), defaults_file=DEFAULTS_FILE)
   if (.not. p%l_visc_3d) error stop 'pm_visc1d: l_visc_3d is .false. in the namelist'

   ! the grid of vilma.x (vilma_drive build_grid)
   nlat = p%nlat;  if (nlat <= 0) nlat = 2*p%lmax + 2
   nphi = p%nphi;  if (nphi <= 0) nphi = 4*p%lmax
   if (p%mmax >= 0 .and. p%eps_polar > 0.0_wp) then
      call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mmax=p%mmax, mres=p%mres, eps_polar=p%eps_polar)
   else if (p%mmax >= 0) then
      call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mmax=p%mmax, mres=p%mres)
   else if (p%eps_polar > 0.0_wp) then
      call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mres=p%mres, eps_polar=p%eps_polar)
   else
      call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mres=p%mres)
   end if

   em = build_earth(p)
   call response_init_ve(resp, em, sht, 0.1_wp*kyr)            ! for the node radii resp%r only
   call load_visc_3d(p, sht, resp%r, visc_node)
   nh = size(visc_node, 1)
   write(*,'(a,i0,a,i0,a,i0,a,i0)') ' lmax=', p%lmax, '  grid ', sht%nphi, ' x ', sht%nlat, '  nodes ', size(resp%r)

   open(newunit=u, file=trim(out), status='replace', action='write')
   write(u,'(a)') 'r_node,depth,log10_visc_mean,log10_visc_min,log10_visc_max'
   do r = 1, size(visc_node, 2)
      write(u,'(es23.16,4(",",es23.16))') resp%r(r), 6371.0e3_wp - resp%r(r), &
           sum(visc_node(:, r)) / real(nh, wp), minval(visc_node(:, r)), maxval(visc_node(:, r))
   end do
   close(u)
   write(*,'(a)') ' wrote '//trim(out)
end program pm_visc1d
