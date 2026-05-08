program baro_v1

use precisions, only: wp, cwp, dp
use config_mod
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
     type(path_config) :: paths
     character(17)               :: save_gridid

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
     call init_paths(paths)
     call make_save_dirs(P%isave, save_gridid, paths)
	 call system('yes | cp control_file.txt ' //  paths%dir_global)

if (len(trim(P%gridid)) < len('0000_00_00_00_00')) then
! nph, latP and omega should remain the same as in that run
	!%================================================
	!% 1) Loading/preparing topography file
	!% 2) Set up the global C-grid and mapping matrices
	!%================================================
		call generate_global_grid(trim(paths%etopo_file), trim(P%topo_file), &
                                  trim(paths%topo_dir_in), trim(paths%topo_dir_out), &
                                  trim(paths%dir_grid), trim(paths%dir_cols), P, GD)
		! Write the parameters of the grid into a specified file
		call write_GD(trim(paths%dir_grid)//GD_file, GD)
	!%================================================
	!% 1) Write H, ta on the u/v/h-grids
	!% 2) Generate and write sparse matrices
	!%================================================
		call generate_global_matrices(trim(paths%dir_grid), trim(paths%dir_cols), &
                                      trim(paths%dir_mats), P, GD)
else
	inquire( file=trim(paths%out_dir)//trim(P%gridid)//'/.', exist=dir_e )
	if  (.not. dir_e ) then
!		Use MATLAB functions to export the data in M-files into binaries
		call system('xterm -e matlab22 -nosplash -nodesktop -logfile remoteAutocode.log -r ' // &
		'"addpath(genpath(''' // trim(paths%matlab_dir) // ''')); save_grid_binary(''' // trim(P%gridid) // '''); exit;" ')

!		if the matlab script failed
		inquire( file=trim(paths%out_dir)//trim(P%gridid)//'/.', exist=dir_e )
		if  (.not. dir_e ) then
			write(*,'(a)') 'Directory /'//trim(P%gridid)//' for the previously generated grid files doesn''t exist.'
			stop
		end if
	end if

	!%================================================
	!% Link to previously generated files (export from matlab is already done)
	!%================================================
     call system('ln -s ' // trim(paths%out_dir)//trim(P%gridid)//'/'//'global/grid/'//'*.dat ' // trim(paths%dir_grid))
     call system('ln -s ' // trim(paths%out_dir)//trim(P%gridid)//'/'//'global/cols/'//'*.dat ' // trim(paths%dir_cols))
     call system('ln -s ' // trim(paths%out_dir)//trim(P%gridid)//'/'//'global/mats/'//'*.dat ' // trim(paths%dir_mats))
     	! Remove unnecessary links
     call system('rm -f ' // trim(paths%dir_mats) // 'mat_*.dat')
		! Read the parameters of the grid into a specified file
     call system('ln -s ' // trim(paths%out_dir)//trim(P%gridid)//'/'//'global/grid/'//'*.txt ' // trim(paths%dir_grid))
     call read_GD(trim(paths%dir_grid)//GD_file, GD)
     ! Update P to match params of the uploaded grid
     P%nph = GD%nph
     P%coor = GD%coor

end if

!%================================================
!% 1) Collect all the sparse matrices together
!% 2) Write mat and bcdiag
!%================================================
	call write_baro_mats(trim(paths%N_data_dir), trim(paths%dir_cols), trim(paths%dir_grid), &
                         trim(paths%dir_mats), trim(P%cpts), P, GD%nu, GD%nv, GD%nh)

!%===================================
!% solve for a global barotropic tide
!%===================================
  call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

if ((P%fr_scheme <= 2).and.(P%sal_scheme <= 1)) then
  call baro_solver_linear(trim(P%cpts), P, GD, trim(paths%dir_grid), trim(paths%dir_cols), &
                          trim(paths%dir_mats), trim(paths%dir_sols))

elseif ((P%fr_scheme == 3).or.(P%sal_scheme >= 2)) then
  call baro_solver(trim(P%cpts), P, GD, trim(paths%dir_grid), trim(paths%dir_cols), &
                   trim(paths%dir_mats), trim(paths%dir_sols), trim(paths%matlab_dir), save_gridid) ! dir_global

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
		call baro_domain_integrals(cpt, P, GD%nu, GD%nv, GD%nh, trim(paths%dir_grid), &
                                   trim(paths%dir_cols), trim(paths%dir_mats), trim(paths%dir_sols), di(ccpt))
		call show_domain_integrals(cpt, di(ccpt))
	enddo
	call save_di(trim(P%cpts), di, 'di.txt', trim(paths%dir_sols))

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
	call clean_files(1, paths)

!==========================================================================================

end program
