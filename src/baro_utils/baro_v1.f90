program baro_v1

use precisions, only: wp, cwp, dp
use baro_integrals
use baro_solver_mod
use generate_matrices
use generate_grid
use my_trigs
use my_sparse
use save_load
use control
use dispmodule


     implicit none

!    IMPORTANT Directories & Files
     character(len=*), parameter :: data_dir = '/home/amt6/amtvl/Desktop/scratch/Tides/data/'
     character(len=*), parameter :: base_dir = data_dir // 'LAG/baro_fd/'
     character(len=*), parameter :: in_dir = base_dir // 'in/'
     character(len=*), parameter :: out_dir = base_dir // 'out/'
     character(17)               :: save_gridid

     character(len = len(out_dir) + 17 + 1 + len('global/grid/'))     :: dir_cols, dir_grid, dir_mats, dir_sols
     character(len = len(out_dir) + 17)     :: dir_global

     character(len=*), parameter :: matlab_dir = '/home/amt6/amtvl/workspace/Matlab/tides/baro_fd/'

     character(len=*), parameter :: etopo_dir = data_dir // 'ETOPO/'
     character(len=*), parameter :: etopo_file = etopo_dir // 'ETOPO1_Ice_g_gmt4.grd'! original ETOPO NetCDF topo file
     character(len=*), parameter :: nocs_dir = data_dir // 'NOCS/'
!     character(len=*), parameter :: nocs_file = nocs_dir // 'nocs_etopo2.nc'! original NOCS NetCDF topo file
!	Processed topo files input/output
     character(len=*), parameter :: topo_dir_out = in_dir // 'topo/etopo/'
     character(len=*), parameter :: topo_dir_in = in_dir // 'topo/nocs/'
     character(len=*), parameter :: topo_file = 'topo_rot_2.0min_pole_15_-40.dat' ! topo_rot_2.0min_pole_10_-30.dat
                                                                                  ! topo_rot_2.0min_pole_15_-40.dat
     character(len=*), parameter :: N_data_dir = in_dir // 'ocean_N/'

     character(len=*), parameter :: P_file = 'control_file.txt'	! contains parameters of the problem and the solver
     character(len=*), parameter :: GD_file = 'grid_file.txt'	! contains grid parameters

     type(params)	:: P
     type(grid_dims):: GD

!**********************************************************************************************************
!**********************************************************************************************************
	 logical			:: dir_e
     integer			:: ccpt, ncpts
     character(len=2)	:: cpt

     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
     integer :: wall_t1, wall_t2, clock_rate, clock_max
     type(domain_integrals), allocatable :: di(:)

!********************************************************************

     call system('clear')
!%================================
!% load parameters of the problem:
!%================================
     call control_file(P_file, P)
!%========================
!% set up file structure:
!%========================
     call make_save_dirs(P%isave, save_gridid, len(out_dir), out_dir, dir_cols, dir_grid, dir_mats, dir_sols, dir_global)
	 call system('yes | cp control_file.txt ' //  dir_global)

if (len(trim(P%gridid)) < len('0000_00_00_00_00')) then
! nph, latP and omega should remain the same as in that run
	!%================================================
	!% 1) Loading/preparing topography file
	!% 2) Set up the global C-grid and mapping matrices
	!%================================================
		call generate_global_grid(etopo_file, topo_file, topo_dir_in, topo_dir_out, dir_grid, dir_cols, P, GD)
		! Write the parameters of the grid into a specified file
		call write_GD(dir_grid//GD_file, GD)
	!%================================================
	!% 1) Write H, ta on the u/v/h-grids
	!% 2) Generate and write sparse matrices
	!%================================================
		call generate_global_matrices(dir_grid, dir_cols, dir_mats, P, GD)
else
	inquire( file=out_dir//trim(P%gridid)//'/.', exist=dir_e )
	if  (.not. dir_e ) then
!		Use MATLAB functions to export the data in M-files into binaries
		call system('xterm -e matlab22 -nosplash -nodesktop -logfile remoteAutocode.log -r ' // &
		'"addpath(genpath(''' // matlab_dir // ''')); save_grid_binary(''' // trim(P%gridid) // '''); exit;" ')
!		print *, ' matlab22 -nosplash -nodesktop -logfile remoteAutocode.log -r ' // &
!		'"addpath(genpath(''' // matlab_dir // ''')); save_grid_binary(''2012_06_13__16_54'')" '

!		if the matlab script failed
		inquire( file=out_dir//trim(P%gridid)//'/.', exist=dir_e )
		if  (.not. dir_e ) then
			write(*,'(a)') 'Directory /'//trim(P%gridid)//' for the previously generated grid files doesn''t exist.'
			stop
		end if
	end if

	!%================================================
	!% Link to previously generated files (export from matlab is already done)
	!%================================================
     call system('ln -s ' // out_dir//trim(P%gridid)//'/'//'global/grid/'//'*.dat ' // dir_grid)
     call system('ln -s ' // out_dir//trim(P%gridid)//'/'//'global/cols/'//'*.dat ' // dir_cols)
     call system('ln -s ' // out_dir//trim(P%gridid)//'/'//'global/mats/'//'*.dat ' // dir_mats)
     	! Remove unnecessary links
     call system('rm -f ' // dir_mats // 'mat_*.dat')
		! Read the parameters of the grid into a specified file
     call system('ln -s ' // out_dir//trim(P%gridid)//'/'//'global/grid/'//'*.txt ' // dir_grid)
     call read_GD(dir_grid//GD_file, GD)
     ! Update P to match params of the uploaded grid
     P%nph = GD%nph
     P%coor = GD%coor

     ! If the available solution is going to be used as an initial guess then copy it
!	call system('cp ' // out_dir//trim(P%gridid)//'/'//'global/sols/'//'*.dat ' // dir_sols)

end if

!%================================================
!% 1) Collect all the sparse matrices together
!% 2) Write mat and bcdiag
!%================================================
	call write_baro_mats(N_data_dir, dir_cols, dir_grid, dir_mats, trim(P%cpts), P, GD%nu, GD%nv, GD%nh)

!%===================================
!% solve for a global barotropic tide
!%===================================
  call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

if ((P%fr_scheme <= 2).and.(P%sal_scheme <= 1)) then
  call baro_solver_linear(trim(P%cpts), P, GD, dir_grid, dir_cols, dir_mats, dir_sols)

elseif ((P%fr_scheme == 3).or.(P%sal_scheme >= 2)) then
  call baro_solver(trim(P%cpts), P, GD, dir_grid, dir_cols, dir_mats, dir_sols, matlab_dir, save_gridid) ! dir_global

else
      print *, "You must choose a friction scheme fr_scheme: 2 or 3"
      stop
end if

  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )

    print *, ""
    call disp ('===>>> Total time spent on solving the system: ' &
			    // trim(ADJUSTL(conv_secs(T2-T1))) // ' CPU, ' &
	            //trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) )))//'  Wall')
    print *, ""

!  %===========================================
!  % calculate integrals to check the solutions
!  % do for each component
!  %===========================================
  call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

	ncpts = len(trim(P%cpts))/2
	allocate(di(ncpts))
	do ccpt = 1, ncpts
		cpt=P%cpts(2*ccpt-1:2*ccpt)
	!	Print the integrals on every iteration
		call baro_domain_integrals(cpt, P, GD%nu, GD%nv, GD%nh, dir_grid, dir_cols, dir_mats, dir_sols, di(ccpt))
		call show_domain_integrals(cpt, di(ccpt))
	enddo
	call save_di(trim(P%cpts), di, 'di.txt', dir_sols)

  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )

    print *, ""
    call disp ('===>>> Calculation of domain integrals: ' &
			    // trim(ADJUSTL(conv_secs(T2-T1))) // ' CPU, ' &
	            //trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) )))//'  Wall')

!    %============================
!    % clean up temprorary files
!    %============================
!	First argument is clean up level: 0 - none, 1 - remove files in 'temp' dirs, 2 - remove ALL files in dir_mats
	call clean_files(1, dir_cols, dir_grid, dir_mats, dir_sols)

!==========================================================================================

end program
