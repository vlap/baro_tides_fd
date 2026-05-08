module control

     use precisions, only: dp, wp
     use config_mod, only: params
     use spline
     use save_load

     implicit none

!********************************************************************
! DEFAULT parameters of the problem and the solver.
! The defaults are over-written by the values from the control file
!********************************************************************
! Params type moved to config_mod.f90


     type tide_params
          real(wp) :: amp,lovef,omega0
     end type


    contains
!********************************************************************
subroutine control_file(file_name, P)

  implicit none

  character(len=*) :: file_name
  ! Control file variables
  type(params) :: P

  character(len=100) :: buffer, label
  integer :: pos, pos_
  integer :: fh = 15
  integer :: ios = 0
  integer :: line = 0



  open(fh, file=file_name)

  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.

  do while (ios == 0)
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        pos_ = scan(buffer, '!')
        pos = scan(buffer, '=')
        !print *, pos_, pos
        if ( (pos>0).and.((pos<pos_).or.(pos_==0)) ) then

	        label = buffer(1:pos-1)
	        if (pos_ > 0) then
	        	buffer = trim(buffer(pos+1:pos_-1))
	        else
	        	buffer = trim(buffer(pos+1:))
	        end if
		!print *, trim(label), trim(buffer)

	        select case (label)
!********************************************************************
! General specifications:
!********************************************************************
	        case ('cpts')
	           read(buffer, *, iostat=ios) P%cpts
	           !print *, 'Read cpts: ', P%cpts
	           P%ncpts = len(trim(P%cpts))/2
	        case ('graphics')
	           read(buffer, *, iostat=ios) P%graphics
	           !print *, 'Read graphics: ', P%graphics
	        case ('messages')
	           read(buffer, *, iostat=ios) P%messages
!	           print *, 'Read messages: ', P%messages
	        case ('isave')
	           read(buffer, *, iostat=ios) P%isave
	           !print *, 'Read isave: ', P%isave
	        case ('cleanup')
	           read(buffer, *, iostat=ios) P%cleanup
	           !print *, 'Read cleanup: ', P%cleanup
	        case ('gridid')
	           read(buffer, *, iostat=ios) P%gridid
	           !print *, 'Read gridid: ', P%gridid
	        case ('use_sol')
	           read(buffer, *, iostat=ios) P%use_sol
	           !print *, 'Read use_sol: ', P%use_sol
!********************************************************************
! Global grid:
!********************************************************************
	        case ('nph')
	           read(buffer, *, iostat=ios) P%nph
	           !print *, 'Read nph: ', P%nph
	        case ('coor')
	           read(buffer, *, iostat=ios) P%coor
	           !print *, 'Read coor: ', P%coor
	        case ('load_etopo')
	           read(buffer, *, iostat=ios) P%load_etopo
	           !print *, 'Read load_etopo: ', P%load_etopo
	        case ('etopo_res')
	           read(buffer, *, iostat=ios) P%etopo_res
	           !print *, 'Read etopo_res: ', P%etopo_res
	        case ('bays')
	           read(buffer, *, iostat=ios) P%bays
	           !print *, 'Read bays: ', P%bays
	        case ('inlandseas')
	           read(buffer, *, iostat=ios) P%inlandseas
	           !print *, 'Read inlandseas: ', P%inlandseas
	        case ('hmin_scheme')
	           read(buffer, *, iostat=ios) P%hmin_scheme
	           !print *, 'Read hmin_scheme: ', P%hmin_scheme
	        case ('hmin')
	           read(buffer, *, iostat=ios) P%hmin
	           !print *, 'Read hmin: ', P%hmin
	        case ('gppw')
	           read(buffer, *, iostat=ios) P%gppw
	           !print *, 'Read gppw: ', P%gppw
!********************************************************************
!  Bottom friction scheme:
!********************************************************************
	        case ('fr_scheme')
	           read(buffer, *, iostat=ios) P%fr_scheme
	           !print *, 'Read fr_scheme: ', P%fr_scheme
	        case ('cd')
	           read(buffer, *, iostat=ios) P%cd
	           !print *, 'Read cd: ', P%cd
	        case ('ubar')
	           read(buffer, *, iostat=ios) P%ubar
	           !print *, 'Read ubar: ', P%ubar
	        case ('Qbar')
	           read(buffer, *, iostat=ios) P%Qbar
	           !print *, 'Read Qbar: ', P%Qbar
	        case ('ndays')
	           read(buffer, *, iostat=ios) P%ndays
	           !print *, 'Read ndays: ', P%ndays
	        case ('nppp')
	           read(buffer, *, iostat=ios) P%nppp
	           !print *, 'Read nppp: ', P%nppp
	        case ('sh_depth')
	           read(buffer, *, iostat=ios) P%sh_depth
	           !print *, 'Read sh_depth: ', P%Qbar
	        case ('coast_dist')
	           read(buffer, *, iostat=ios) P%coast_dist
	           !print *, 'Read coast_dist: ', P%coast_dist
!********************************************************************
!  SAL scheme:
!********************************************************************
	        case ('sal_scheme')
	           read(buffer, *, iostat=ios) P%sal_scheme
	           !print *, 'Read sal_scheme: ', P%sal_scheme
	        case ('beta0')
	           read(buffer, *, iostat=ios) P%beta0
	           !print *, 'Read beta: ', P%beta0
	        case ('betamin')
	           read(buffer, *, iostat=ios) P%betamin
	           !print *, 'Read betamin: ', P%betamin
	        case ('betamax')
	           read(buffer, *, iostat=ios) P%betamax
	           !print *, 'Read betamax: ', P%betamax
	        case ('ntrunc')
	           read(buffer, *, iostat=ios) P%ntrunc
	           !print *, 'Read ntrunc: ', P%ntrunc
	        case ('save_sht')
	           read(buffer, *, iostat=ios) P%save_sht
	           !print *, 'Read save_sht: ', P%save_sht
	        case ('sal_avg')
	           read(buffer, *, iostat=ios) P%sal_avg
	           !print *, 'Read sal_avg: ', P%sal_avg
!********************************************************************
!  ITD scheme:
!********************************************************************
	        case ('itd_scheme')
	           read(buffer, *, iostat=ios) P%itd_scheme
!	           print *, 'Read itd_scheme: ', P%itd_scheme
	        case ('N_form')
	           read(buffer, *, iostat=ios) P%N_form
!	           print *, 'Read N_form: ', P%N_form
	        case ('Ns')
	           read(buffer, *, iostat=ios) P%Ns
!	           print *, 'Read N on the surface: ', P%Ns
	        case ('Nl')
	           read(buffer, *, iostat=ios) P%Nl
!	           print *, 'Read stratification lengthscale: ', P%Nl
	        case ('itd_coeff')
	           read(buffer, *, iostat=ios) P%itd_coeff
!	           print *, 'Read a multiplying coeff in the ITD parametrization: ', P%itd_coeff
	        case ('trapped')
	           read(buffer, *, iostat=ios) P%trapped
!	           print *, 'Read trapped at crit lat: ', P%trapped
	        case ('sht_smooth_H')
	           read(buffer, *, iostat=ios) P%sht_smooth_H
!	           print *, 'Read truncation order for SHT of H: ', P%sht_smooth_H
	        case ('smooth_type')
	           read(buffer, *, iostat=ios) P%smooth_type
!	           print *, 'Read smoothing type of H: ', P%smooth_type
            case ('baro_on_smoothed')
               read(buffer, *, iostat=ios) P%baro_on_smoothed
            case ('topo_file')
               read(buffer, *, iostat=ios) P%topo_file
!	        case ('smooth_save_load')
!	           read(buffer, *, iostat=ios) P%smooth_save_load
!	           print *, 'Read save/load type of the smoothed H: ', P%smooth_save_load

	        case ('N_data')
	           read(buffer, *, iostat=ios) P%N_data
!	           print *, 'Read N data file: ', P%N_data
!********************************************************************
!   Nonlinear solver iterative scheme:
!********************************************************************
	        case ('lib')
	           read(buffer, *, iostat=ios) P%lib
!	           print *, 'Read lib: ', P%lib
	        case ('omp_num_threads')
	           read(buffer, *, iostat=ios) P%omp_num_threads
!	           print *, 'Read omp_num_threads: ', P%omp_num_threads
	        case ('mpi_num_nodes')
	           read(buffer, *, iostat=ios) P%mpi_num_nodes
!	           print *, 'Read mpi_num_nodes: ', P%mpi_num_nodes
	        case ('blas_num_threads')
	           read(buffer, *, iostat=ios) P%blas_num_threads
!	           print *, 'Read blas_num_threads: ', P%blas_num_threads
	        case ('solver')
	           read(buffer, *, iostat=ios) P%solver
!	           print *, 'Read solver: ', P%solver
	        case ('gmres_prec')
	           read(buffer, *, iostat=ios) P%gmres_prec
!	           print *, 'Read gmres_prec: ', P%gmres_prec
	        case ('gmres_rest_it')
	           read(buffer, *, iostat=ios) P%gmres_rest_it
!	           print *, 'Read gmres_rest_it: ', P%gmres_rest_it
	        case ('gmres_tol')
	           read(buffer, *, iostat=ios) P%gmres_tol
!	           print *, 'Read gmres_tol: ', P%gmres_tol
	        case ('pardiso_symbolic')
	           read(buffer, *, iostat=ios) P%pardiso_symbolic
!	           print *, 'Read pardiso_symbolic: ', P%pardiso_symbolic
	        case ('pardiso_iterative')
	           read(buffer, *, iostat=ios) P%pardiso_iterative
!	           print *, 'Read pardiso_iterative: ', P%pardiso_iterative
	        case ('pardiso_iter_tol')
	           read(buffer, *, iostat=ios) P%pardiso_iter_tol
!	           print *, 'Read pardiso_iter_tol: ', P%pardiso_iter_tol
	        case ('pardiso_ooc')
	           read(buffer, *, iostat=ios) P%pardiso_ooc
!	           print *, 'Read pardiso_ooc: ', P%pardiso_ooc
	        case ('pardiso_max_ram')
	           read(buffer, *, iostat=ios) P%pardiso_max_ram
!	           print *, 'Read pardiso_max_ram: ', P%pardiso_max_ram
	        case ('pardiso_max_swap')
	           read(buffer, *, iostat=ios) P%pardiso_max_swap
!	           print *, 'Read pardiso_max_swap: ', P%pardiso_max_swap
	        case ('pardiso_keep_files')
	           read(buffer, *, iostat=ios) P%pardiso_keep_files
!	           print *, 'Read pardiso_keep_files: ', P%pardiso_keep_files

	        case ('cvg_scheme')
	           read(buffer, *, iostat=ios) P%cvg_scheme
!	           print *, 'Read cvg_scheme: ', P%cvg_scheme
	        case ('p_fric')
	           read(buffer, *, iostat=ios) P%p_fric
!	           print *, 'Read p_fric: ', P%p_fric
	        case ('cvg_dhbar')
	           read(buffer, *, iostat=ios) P%cvg_dhbar
!	           print *, 'Read domain averaged requirement on the residual: ', P%cvg_dhbar
	        case ('cvg_dhmax')
	           read(buffer, *, iostat=ios) P%cvg_dhmax
!	           print *, 'Read pointwise requirement on the residual: ', P%cvg_dhmax
	        case ('p_avg')
	           read(buffer, *, iostat=ios) P%p_avg
!	           print *, 'Read averaging between iterations: ', P%p_avg
!********************************************************************
! Physical Parameters:
!********************************************************************
	        case ('latP')
	           read(buffer, *, iostat=ios) P%latP
	           !print *, 'Read latP: ', P%latP
	        case ('lonP')
	           read(buffer, *, iostat=ios) P%lonP
	           !print *, 'Read lonP: ', P%lonP
	        case ('re')
	           read(buffer, *, iostat=ios) P%re
	           !print *, 'Read re: ', P%re
	        case ('omega')
	           read(buffer, *, iostat=ios) P%omega
	           !print *, 'Read omega: ', P%omega
	        case ('g')
	           read(buffer, *, iostat=ios) P%g
	           !print *, 'Read g: ', P%g
	        case ('rhoe')
	           read(buffer, *, iostat=ios) P%rhoe
	           !print *, 'Read rhoe: ', P%rhoe
	        case ('rhoo')
	           read(buffer, *, iostat=ios) P%rhoo
	           !print *, 'Read rhoo: ', P%rhoo
!********************************************************************
	        case default
	           !print *, 'Skipping invalid label at line', line
	        end select
        end if
     end if
  end do

end subroutine

!**********************************************************************************************************
function get_params(ncpts, cpts)

	implicit none

     integer, intent(in)				:: ncpts
     character(len=2*ncpts), intent(in)   :: cpts

	 integer 			:: ccpt
     character(len=2)	 :: cpt
     type (tide_params) :: get_params(ncpts)

	do ccpt = 1, ncpts
		cpt = cpts(2*ccpt-1:2*ccpt)
		get_params(ccpt) = get_pars(cpt)
	enddo

end function get_params

function get_pars(cpt)

	implicit none

     character(len=2)   :: cpt
     type (tide_params) :: get_pars

!% returns the tidal get_pars%amplitude, love factor, and frequency
!% for the specified tidal component cpt, which is one of
!%
!%  m2, k1, s2, o1, p1, n2, mf, k2, mm, q1

select case(cpt)

       case ('m2')
       get_pars%amp=0.242334
       get_pars%lovef=0.693
       get_pars%omega0=1.405189d-4

       case ('k1')
       get_pars%amp=0.141565
       get_pars%lovef=0.736
       get_pars%omega0=0.7292117d-4

       case ('s2')
       get_pars%amp=0.112743
       get_pars%lovef=0.693
       get_pars%omega0=1.454441d-4

       case ('o1')
       get_pars%amp=0.100661
       get_pars%lovef=0.695
       get_pars%omega0=0.6759774d-4

       case ('p1')
       get_pars%amp=0.046848
       get_pars%lovef=0.706
       get_pars%omega0=0.7252295d-4

       case ('n2')
       get_pars%amp=0.046397
       get_pars%lovef=0.693
       get_pars%omega0=1.378797d-4

       case ('mf')
       get_pars%amp=0.042041
       get_pars%lovef=0.693
       get_pars%omega0=0.053234d-4

       case ('k2')
       get_pars%amp=0.030684
       get_pars%lovef=0.693
       get_pars%omega0=1.458423d-4

       case ('mm')
       get_pars%amp=0.022191
       get_pars%lovef=0.693
       get_pars%omega0=0.026392d-4

       case ('q1')
       get_pars%amp=0.019273
       get_pars%lovef=0.695
       get_pars%omega0=0.6495854d-4

       case ('t2')
       get_pars%amp=0.1
       get_pars%lovef=0.693
       get_pars%omega0=0

       case default
       get_pars%amp=0
       get_pars%lovef=0
       get_pars%omega0=0
       write(*, '("** Need to enter a valid tidal component.")')
       write(*, '("** See function get_pars")')

end select

end function get_pars
!**********************************************************************************************************

subroutine get_load_love(n, nvals, h, k)
! given n, returns the load love numbers h'_n and k'_n.
	implicit none

	integer	:: n, nvals(n)
	real ( kind = 8 ), allocatable	:: nd(:), hd(:), kd(:), ypp(:)
	real ( kind = 8 )				:: ypval, yppval, zero_dp = 0

	real(wp)				:: h(n), k(n)
	real ( kind = 8 )		:: h_tmp(n), k_tmp(n)
	integer	:: d = size( (/ 1, 2, 3, 4, 5, 6, 8, 10, 18, 32, 56, 100, 180, 325, 550, 1000, 1800, 3000, 10000 /) )



	integer :: istat, j


    allocate(nd(d), hd(d), kd(d), ypp(d), stat = istat)
nd = (/ 1, 2, 3, 4, 5, 6, 8, 10, 18, 32, 56, 100, 180, 325, 550, 1000, 1800, 3000, 10000 /)
hd = - (/ 0.290, 1.001, 1.052, 1.053, 1.088, 1.147, 1.291, 1.433, 1.893, 2.379, 2.753, 3.058, 3.474, 4.107, &
		  4.629, 4.906, 4.953, 4.954, 4.956 /)
kd = -(/ 0., 0.615, 0.585, 0.528, 0.516, 0.535, 0.604, 0.682, 0.952, 1.240, 1.402, 1.461, 1.591, 1.928, 2.249,&
		  2.431, 2.465, 2.468, 2.469 /)


!call interp_lagrange ( dim_num, d, real(nd,wp), hd, n, real(nvals,wp), h )
!call interp_lagrange ( dim_num, d, real(nd,wp), kd, n, real(nvals,wp), k )

call spline_cubic_set ( d, nd, hd, 0, zero_dp, 0, zero_dp, ypp ) !0 - the spline should be a quadratic over the first and last intervals
do j = 1,n
	call spline_cubic_val ( d, nd, hd, ypp, real(nvals(j),kind=8), h_tmp(j), ypval, yppval )
end do

call spline_cubic_set ( d, nd, kd, 0,zero_dp, 0,zero_dp, ypp ) !0 - the spline should be a quadratic over the first and last intervals
do j = 1,n
	call spline_cubic_val ( d, nd, kd, ypp, real(nvals(j), kind=8), k_tmp(j), ypval, yppval )
end do
k_tmp = k_tmp/nvals

! back to the variables in our precision
h = h_tmp
k = k_tmp

!call save_vector(nd, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/out/0000_00_00__00_00/global/grid/'&
!     					 // 'nd' // '.dat')
!call save_vector(hd, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/out/0000_00_00__00_00/global/grid/'&
!     					 // 'hd' // '.dat')
!call save_vector(nvals, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/out/0000_00_00__00_00/global/grid/'&
!     					 // 'nvals' // '.dat')
!call save_vector(h, '/home/amsta/amtvl/Desktop/scratch/Tides/data/LAG/baro_fd/out/0000_00_00__00_00/global/grid/'&
!     					 // 'hvals' // '.dat')

end subroutine get_load_love


end module
