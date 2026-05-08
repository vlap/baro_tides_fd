module config_mod
    use precisions, only: dp, wp
    implicit none

    ! Originally, these parameters were part of the 'control' module in control_file.f90.
    ! They have been moved here to centralize configuration and decouple it from IO logic.
    type params
        integer :: isave = 0
        integer :: cleanup = 1
        integer :: graphics = 0
        integer :: messages = 0
        character(len=100) :: cpts = 'm2'
        integer :: ncpts = 1
        character(len=17) :: gridid = '0'
        integer :: use_sol
        integer :: nph = 720
        integer :: coor = 1
        integer :: load_etopo = 0
        integer :: etopo_res = 5400
        integer :: bays = 1
        integer :: inlandseas = 1
        integer :: hmin_scheme = 1
        real(dp) :: hmin = 10
        integer :: gppw = 4
        integer :: fr_scheme = 2
        real(dp) :: cd = 0.0025
        real(dp) :: ubar = 1
        real(dp) :: Qbar = 100
        integer :: ndays = 0
        integer :: nppp = 0
        real(dp) :: sh_depth = 1000
        real(dp) :: coast_dist = 3e5
        integer :: sal_scheme = 1
        real(dp) :: beta0 = 0.085
        real(dp) :: betamin = 0
        real(dp) :: betamax = 0.3
        integer :: ntrunc = 360
        integer :: save_sht = 0
        real(dp) :: sal_avg = 10
        integer :: itd_scheme = 0
        integer :: N_form = 1
        real(dp) :: Ns = .02
        real(dp) :: Nl = 500
        real(dp) :: itd_coeff = 1
        integer :: trapped = 1
        integer :: sht_smooth_H = 720
        integer :: smooth_type = 1
        integer :: baro_on_smoothed = 0
        character(len=100) :: topo_file = 'topo_rot_2.0min_pole_15_-40.dat'
        character(len=50) :: N_data = 'woa05_1deg_pole_15_-40'
        character :: lib = 'm'
        integer :: omp_num_threads = 0
        integer :: mpi_num_nodes = 0
        integer :: blas_num_threads = 0
        character(len=7) :: solver = 'pardiso'
        character(len=7) :: gmres_prec = 'ilut'
        integer :: gmres_rest_it = 3
        integer :: gmres_tol = 6
        integer :: pardiso_iterative = 0
        integer :: pardiso_iter_tol = 6
        integer :: pardiso_symbolic = 1
        integer :: pardiso_ooc = 1
        integer :: pardiso_max_ram = 20*1024
        integer :: pardiso_max_swap = 0
        integer :: pardiso_keep_files = 0
        integer :: cvg_scheme = 1
        real(dp) :: cvg_dhbar = 0.01
        real(dp) :: cvg_dhmax = 0.1
        integer :: p_fric = 1
        real(dp) :: p_avg = 0.5
        real(dp) :: latP = 15
        real(dp) :: lonP = -40
        real(dp) :: re = 6.371e6
        real(dp) :: omega = 7.292115e-5
        real(dp) :: g = 9.80665
        real(dp) :: rhoe = 5515
        real(dp) :: rhoo = 1030
    end type params

    ! Directory structure configuration.
    ! Originally hardcoded in baro_v1.f90, now centralized for portability.
    type path_config
        character(len=256) :: data_dir = './data/'
        character(len=256) :: base_dir
        character(len=256) :: in_dir
        character(len=256) :: out_dir
        character(len=256) :: matlab_dir = './matlab/'
        character(len=256) :: etopo_dir
        character(len=256) :: etopo_file
        character(len=256) :: nocs_dir
        character(len=256) :: topo_dir_out
        character(len=256) :: topo_dir_in
        character(len=256) :: N_data_dir
        
        ! Runtime derived paths
        character(len=256) :: dir_cols, dir_grid, dir_mats, dir_sols, dir_global
    end type path_config

contains

    subroutine init_paths(paths)
        type(path_config), intent(inout) :: paths
        
        paths%base_dir = trim(paths%data_dir) // 'LAG/baro_fd/'
        paths%in_dir = trim(paths%base_dir) // 'in/'
        paths%out_dir = trim(paths%base_dir) // 'out/'
        paths%etopo_dir = trim(paths%data_dir) // 'ETOPO/'
        paths%etopo_file = trim(paths%etopo_dir) // 'ETOPO1_Ice_g_gmt4.grd'
        paths%nocs_dir = trim(paths%data_dir) // 'NOCS/'
        paths%topo_dir_out = trim(paths%in_dir) // 'topo/etopo/'
        paths%topo_dir_in = trim(paths%in_dir) // 'topo/nocs/'
        paths%N_data_dir = trim(paths%in_dir) // 'ocean_N/'
    end subroutine init_paths

    subroutine make_save_dirs(isave, now, paths)
        integer, intent(in) :: isave
        character(17), intent(out) :: now
        type(path_config), intent(inout) :: paths
        
        if (isave == 1) then
            call time_and_date(trim(paths%out_dir), now)
        else
            now = '0000_00_00__00_00'
        endif

        paths%dir_global = trim(paths%out_dir) // now // '/'
        paths%dir_grid = trim(paths%dir_global) // 'global/grid/'
        paths%dir_cols = trim(paths%dir_global) // 'global/cols/'
        paths%dir_mats = trim(paths%dir_global) // 'global/mats/'
        paths%dir_sols = trim(paths%dir_global) // 'global/sols/'

        if (isave == 1) then
            write(*, '("Saving data as ", a)') now
        else
            call system('rm -rf ' // trim(paths%dir_global))
        endif
        
        call system('mkdir -p ' // trim(paths%dir_grid))
        call system('mkdir -p ' // trim(paths%dir_grid) // 'temp/')
        call system('mkdir -p ' // trim(paths%dir_cols))
        call system('mkdir -p ' // trim(paths%dir_cols) // 'temp/')
        call system('mkdir -p ' // trim(paths%dir_mats))
        call system('mkdir -p ' // trim(paths%dir_mats) // 'temp/')
        call system('mkdir -p ' // trim(paths%dir_sols))
        call system('mkdir -p ' // trim(paths%dir_sols) // 'temp/')
    end subroutine make_save_dirs

    subroutine clean_files(level, paths)
        integer, intent(in) :: level
        type(path_config), intent(in) :: paths
        
        if (level > 0) then
            write(*, '("Cleaning up temporary files ")')
            call system('rm -rf ' // trim(paths%dir_cols) // 'temp/')
            call system('rm -rf ' // trim(paths%dir_grid) // 'temp/')
            call system('rm -rf ' // trim(paths%dir_mats) // 'temp/')
            call system('rm -rf ' // trim(paths%dir_sols) // 'temp/')
        endif

        if (level > 1) then
            call system('rm -rf ' // trim(paths%dir_mats))
        endif
    end subroutine clean_files

    ! Helper copied from save_load for path logic
    subroutine time_and_date(out_dir, now)
        character(len=*), intent(in) :: out_dir
        character(8)  :: date
        character(10) :: time
        character(17) :: now
        logical :: dir_e
        integer :: mins, delay=1
        call date_and_time(DATE=date,TIME=time)
        now = date(1:4) // '_' // date(5:6) // '_' // date(7:8) // '__' // time(1:2) //'_' // time(3:4)
        inquire( file=trim(out_dir)//now//'/.', exist=dir_e )
        do while (dir_e)
            read( time(3:4), '(i2.2)' ) mins
            mins = mins + delay
            write( time(3:4), '(i2.2)' ) mins
            now = date(1:4) // '_' // date(5:6) // '_' // date(7:8) // '__' // time(1:2) //'_' // time(3:4)
            inquire( file=trim(out_dir)//now//'/.', exist=dir_e )
            delay = delay + 1
        enddo
    end subroutine time_and_date

end module config_mod
