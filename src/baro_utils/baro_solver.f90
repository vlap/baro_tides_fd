module baro_solver_mod

!     use my_trigs
!     use generate_global_matrices
	 use generate_grid
	 use control
	 use baro_integrals
     use baro_projection
     use unsym_solvers
     use iter_solvers
     use my_sparse_aggregate
     use my_sparse
     use sal
     use itd
     use save_load
     use dispmodule
     use precisions, only: wp, cwp

     !use f90_kind      ! Kind-definition MODULE of the F-compiler
!     use sparse_utils  ! module for sparse matrix operations

!       External Subroutines
!      external   dzasum
!      real       dzasum

!       Solver
!     character(len=*), parameter, private :: solver='pardiso'! umfpack or pardiso
!     integer, parameter, private  :: p_fric = 1 ! 1 for friction; 2 for over-relaxation
!     integer, parameter, private  :: f = 1      ! 1 for 4-element f; 2/3 for iterative treatment
!     real(wp), parameter, private :: p_avg = 0.5! for averaging between iterations; 1 is no averaging
!     integer, parameter, private  :: cvg_scheme = 1! 1 for complex dh, 2 for abs(dh)
!     real(wp), parameter, private :: cvg_dhbar=0.01!  % domain averaged requirement
!     real(wp), parameter, private :: cvg_dhmax=0.1 !    % pointwise requirement

     contains

!==========================================================================================
subroutine baro_solver(cpts, P, GD, dir_grid, dir_cols, dir_mats, dir_sols, matlab_dir, save_gridid) ! dir_global

    implicit none
!%   Calculates a barotropic tide (iterative). Uses u/v/h + volume transport formulation
!%
!%   H_u : H on the u-grid points
!%   H_v : H on the v-grid points
!%   H_h : H on the h-grid points
!%   ta_u : tau on the u-grid points
!%   ta_v : tau on the v-grid points
!%   ta_h : tau on the h-grid points
type(disp_settings) ds

     integer ::     nph, nth, coor
     integer :: nu, nv, nh, np
     real(wp) :: latP, lonP
     real(wp) :: g, re, cdg, Q
     real(wp), allocatable :: beta(:), beta0(:)
!     complex(cwp) :: beta0_cmplx
!     complex(cwp), allocatable :: beta_cmplx(:)
     character(len = *) :: dir_grid, dir_cols, dir_mats, dir_sols, matlab_dir, save_gridid!, dir_global
     type (tide_params) :: pars
     type(params)	:: P
     type(grid_dims):: GD

	 logical :: 	 CGS_ready = .false.
     integer ::      j, istat, itn, cvgd, l,m,n
     real    ::      T1_itn, T2_itn, T1_cpt, T2_cpt ! for measuring CPU (NOT REAL TIME!)
     integer :: 	 wall_t1_itn, wall_t2_itn, wall_t1_cpt, wall_t2_cpt, clock_rate, clock_max
     integer, allocatable :: hp(:, :), up(:, :), vp(:, :)

     character(len=*) :: cpts
     character(len=2) :: cpt
     integer :: ncpts, ccpt
     integer :: iu(GD%nu), iv(GD%nv), ih(GD%nh)

     complex(cwp), allocatable :: tmp_u(:), tmp_v(:)
     complex(cwp), allocatable :: hsal_res(:)
     type (triplet) :: GHU, GHV
	 type (triplet) :: DUU, DVU, DUV, DVV

	 type (csr) :: tmp_csr
	 type (csr_cmplx) :: tmp_csr_cmplx

     type (triplet_cmplx) :: mat_coo_cmplx
     type (csr) :: mat_csr
     type (csr_cmplx) :: mat_csr_cmplx
     type (csc_cmplx) :: mat_csc_cmplx

     type (csr_cmplx) :: speye
     integer, allocatable :: bcdiag(:)
     complex(cwp), allocatable, dimension(:) :: rhs

     complex(cwp), allocatable :: uvh(:), u(:), v(:), h(:), hsal(:)

     complex(cwp), allocatable :: Du(:), Dv(:), Qu(:), Qv(:), tmp(:)
	 real(wp), allocatable :: KU(:), KV(:)
     real(wp), allocatable     :: H_u(:), H_v(:)!, H_h(:)

     real(wp), allocatable     :: delta(:, :, :) !delta(itn,type,ccpt)

     type(domain_integrals) :: di

!==========================================================================================
! DISP MODULE SETTINGS
!==========================================================================================
call tostring_set(rfmt='F12.1')
!==========================================================================================
! SHORTCUTS
!==========================================================================================
latP = P%latP
lonP = P%lonP
g = P%g
re = P%re
cdg = P%cd
Q = P%Qbar

nph = GD%nph
nth = GD%nta
nu = GD%nu
nv = GD%nv
nh = GD%nh
np = GD%np
!==========================================================================================

write(*, '("====================================")')
write(*, '("Welcome to baro_solver (nonlinear)")')
write(*, '("====================================")')

itn=1         ! iteration number
cvgd=0        ! switched to 1 when convergence occurs

ncpts=len(cpts)/2

     iu = (/ ( (j), j=1,nu ) /)
     iv = (/ ( (nu + j), j=1,nv ) /)
     ih = (/ ( (nu + nv + j), j=1,nh ) /)

if (P%sal_scheme >= 2) then
	allocate(beta0(ncpts), stat = istat)
	beta0 = P%beta0
endif

do while ( (cvgd == 0) .and. (itn < 100)) !.and. (maxval(delta(itn,3:6,:)) < 1000) )

	allocate(delta(itn, 8, len(cpts)/2), stat = istat)
	delta=0.
	if (itn > 1) then
	    call load_matrix_3d(delta(1:itn-1, :, :), dir_sols // 'delta.dat')
	endif

	  call disp('==============')
	  call disp(' Iteration ', itn)
	  call disp('==============')


	!%==========================
	!% solve for each component
	!%==========================
	!      call disp('Solving ' // tostring(mat_csr_cmplx%ni) //' x '// tostring(mat_csr_cmplx%nj)//' system with '// &
	!                 tostring(mat_csr_cmplx%nz)// ' entries for:' )
	call CPU_Time(T1_itn)
	call system_clock ( wall_t1_itn, clock_rate, clock_max )

	do ccpt = 1, ncpts

	    cpt=cpts(2*ccpt-1:2*ccpt)

	!%==================================================
	!% Load the basic sparse matrix generated by baro_uvhT_mat
	!% (includes terms: rotation, c.o.mass, itd and RHS (forcing))
	!%==================================================
	!     call load_alloc_cmat(mat, dir_mats // 'temp/mat_init.dat')
	     ! Convert mat from COO to CSR:
	!     call coo2csr(mat, mat_csr, P%lib)
	!     call dealloc_sparse(mat)
		if (P%itd_scheme /= 1) then
			call load_alloc_sparse(mat_csr, dir_mats // 'temp/' // 'mat_csr_init.dat')
		else
			call load_alloc_sparse(mat_csr, dir_mats // 'temp/' // cpt // '_mat_csr_init.dat')
		endif
			!  now transform into cmplx matrix
			call init_sparse_vals(mat_csr_cmplx,mat_csr%indi,mat_csr%indj,mat_csr%vals, &
	                                mat_csr%ni, mat_csr%nj, mat_csr%nz)
			call dealloc_sparse(mat_csr)

	!  %=============
	!  % impose bcs:
	!  %=============
	call load_alloc_vector(bcdiag, dir_mats // 'temp/bcdiag.dat')
	!mat=bcmat*mat; !fill the row of mat with zeroes for boundary variables (then will get -i*omega)
	!rhs=bcmat*rhs;
	     call dia_csr_mul(real(bcdiag,wp), mat_csr_cmplx, skit)
	     deallocate(bcdiag)

	!     mat = mat-i*omega0*speye(np)
	      pars = get_pars(cpt)
	      speye = speye_csr_cmplx(np, cmplx(0., -pars%omega0, kind=cwp) )
	      call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), speye, mkl)
	      call dealloc_sparse(speye)

	!  %=====================
	!  % Print matrix dimensions
	!  %=====================
	    if ((itn == 1).and.(ccpt == 1)) then
	    	call tostring_set(rfmt='ES0.2')
	      	call disp('Have [' // tostring(real(mat_csr_cmplx%ni)) //' x '// tostring(real(mat_csr_cmplx%nj))//'] sparse matrix ('// &
	                 'nnz = ' // tostring(real(mat_csr_cmplx%nz))//'). Solve with '//P%solver// '.' )
			call tostring_set_factory()
	    endif

	write(*, '(i2, ". for ", a, ": assembling the matrix and rhs, ")', advance='no') ccpt, cpt
	!%=======================
	!% calculate rhs forcing
	!%=======================
		if (itn == 1) then
			call baro_rhs(rhs, cpt, P, GD, dir_cols, dir_grid, dir_mats)
			call save_vector(rhs, dir_cols //'temp/'// cpt // '_rhs_0.dat')
		else
			call load_alloc_vector(rhs, dir_cols //'temp/'// cpt // '_rhs_0.dat')
		endif
	!==========================================================================================
	!    %==============================================================
	!    % add friction  (already included if fr_scheme = 0 or 1 or 2)
	!    %==============================================================
	     if (P%fr_scheme >= 3) then

			call load_alloc_vector(KU, dir_cols // 'KU.dat')
			call load_alloc_vector(KV, dir_cols // 'KV.dat')
			allocate(Qu(nu), Qv(nv), stat = istat)

			if (itn == 1) then ! Use linear friction option 2 (average value for transport Q)
				Q = P%Qbar
				call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
				call load_alloc_vector(H_v, dir_cols // 'H_v.dat')
				Qu = cdg*Q/H_u**2 / KU
				Qv = cdg*Q/H_v**2 / KV
				deallocate(H_u, H_v)
			else ! nonlinear bottom drag
			     call load_alloc_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
			     call load_alloc_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
			     call load_alloc_vector(Du, dir_sols // cpt // '_Du' // '_m0_drag' // '.dat')
			     call load_alloc_vector(Dv, dir_sols // cpt // '_Dv' // '_m0_drag' // '.dat')

			   where (u==0)
				   Qu=Du
			   elsewhere
				   Qu=Du/u
			    end where
			   where (v==0)
				   Qv=Dv
			   elsewhere
				   Qv=Dv/v
			    end where
	!		   Du=Du-P%p_fric*Qu*u
	!		   Dv=Dv-P%p_fric*Qv*v
	!		   Qu=P%p_fric*Qu
	!		   Qv=P%p_fric*Qv

			     rhs(iu) = rhs(iu) - KU*(Du - Qu*u)
			     rhs(iv) = rhs(iv) - KV*(Dv - Qv*v)

			   deallocate(u,v, Du, Dv)
			endif

	!       mat = mat + sparse(iu,iu,Qu./cosh(ta_u),np,np) + sparse(iv,iv,Qv./cosh(ta_v),np,np);
	!	    call init_sparse_vals(speye,(/ ( (j), j=1,nu+nv+1), ((nu+nv+1), j=nu+nv+1, np+1) /),(/ ( (j), j=1,nu+nv) /), &
	!		           (/ KU*Qu, KV*Qv/), np, np, nu+nv) ! Initialize diagonal CSR matrix
	!		      call dealloc_sparse(speye)
	!		call disp(speye%vals, orient='row')
			allocate(tmp(nh), stat=istat)
			tmp = 0.
			call csr_dia_add(mat_csr_cmplx, (/  KU*Qu, KV*Qv, tmp/))
			deallocate(tmp)

			deallocate(KU, KV, Qu, Qv)
		endif
	!==========================================================================================

	!================================================================================================
	!    % add self-attraction and loading (already included if sal_scheme = 0 or 1)
	!================================================================================================
		if (P%sal_scheme >= 2) then
			!	Define: beta and hsal_res
			allocate(hsal_res(nh), beta(nh), stat = istat)

			if (itn == 1) then
				hsal_res = 0
				beta = P%beta0
			else
				call load_alloc_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat')
				call load_vector(hsal_res, dir_sols //'temp/'// cpt // '_hsal' // '.dat') ! load hsal vector into hsal_res
				if (P%sal_scheme == 3) then
					call load_vector(beta, dir_sols //'temp/'// cpt // '_betasal' // '.dat')
				elseif (P%sal_scheme == 2) then
					beta = beta0(ccpt)
				endif

				hsal_res = hsal_res - beta * h

				deallocate(h)
		    endif
		!******************************************************************
		!		Add the corresponding terms into the linear matrix and RHS
		!		UVH-formulation only
		!******************************************************************
			!*********************************
			!	RHS
			!*********************************
			call load_alloc_sparse(GHU, dir_mats // 'GHU.dat')
			call load_alloc_sparse(GHV, dir_mats // 'GHV.dat')

			if (itn > 1) then
	      		allocate(tmp_u(nu), tmp_v(nv), stat = istat)

				call coo_vec_mul(GHU, hsal_res, P%lib, tmp_u)
				rhs(iu) = rhs(iu) + tmp_u
				call coo_vec_mul(GHV, hsal_res, P%lib, tmp_v)
				rhs(iv) = rhs(iv) + tmp_v

				deallocate(tmp_u, tmp_v)
			end if
			!*********************************
			!	Matrix
			!*********************************
			!	pressure gradients: d/dph terms in the matrix
			    call right_mult_diag(GHU, - beta, nh)
			    ! Transform GHU into a (np x np) matrix.
			    GHU%ni = np
			    GHU%nj = np
			    GHU%indj = GHU%indj + nu + nv ! (move elements to the "h" columns of the matrix)

				! Convert the matrix from COO to CSR:
			    call coo2csr(GHU, tmp_csr, P%lib)
!print *, "hi 1.5", tmp_csr%ni, tmp_csr%nj, tmp_csr%nz
			    call init_sparse_csr_cmplx(tmp_csr_cmplx,tmp_csr%indi,tmp_csr%indj,tmp_csr%vals, &
				                           tmp_csr%ni, tmp_csr%nj, tmp_csr%nz)
!print *, "hi 1.75"
			    call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), tmp_csr_cmplx, mkl)
!print *, "hi 1.95"
			    call dealloc_sparse(tmp_csr)
			    call dealloc_sparse(tmp_csr_cmplx)
!print *, "hi 2"
			!	pressure gradients: d/dta terms in the matrix
			    call right_mult_diag(GHV, - beta, nh)
			    ! Transform GHV into a (np x np) matrix.
			    GHV%ni = np
			    GHV%nj = np
			    GHV%indi = GHV%indi + nu ! (move elements to the "v" rows of the matrix)
			    GHV%indj = GHV%indj + nu + nv ! (move elements to the "h" columns of the matrix)

				! Convert the matrix from COO to CSR:
			    call coo2csr(GHV, tmp_csr, P%lib)
			    call init_sparse_csr_cmplx(tmp_csr_cmplx,tmp_csr%indi,tmp_csr%indj,tmp_csr%vals, &
				                           tmp_csr%ni, tmp_csr%nj, tmp_csr%nz)
			    call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), tmp_csr_cmplx, mkl)
			    call dealloc_sparse(tmp_csr)
			    call dealloc_sparse(tmp_csr_cmplx)
	!******************************************************************

		    call dealloc_sparse(GHU)
		    call dealloc_sparse(GHV)

			deallocate(hsal_res, beta)

		end if

	!==========================================================================================

	!==========================================================================================
	!    %================================================
	!    % add something to do with Coriolis matrices...
	!    %================================================

	!    if (f) > 1
	!      if itn > 1
	!        load([flags.dir.file,'/',domain,'/',cpt],'u','v');
	!        load(matfile,'v2uf','v2ufsp','u2vf','u2vfsp');
	!        rhs(iu)=rhs(iu)+(v2uf-v2ufsp)*v;
	!        rhs(iv)=rhs(iv)-(u2vf-u2vfsp)*u;
	!        clear u v v2uf v2ufsp
	!      end
	!    end

	!  %=============
	!  % impose bcs:
	!  %=============
	write(*, '("implementing b.c. ")')!, advance = 'no')

	call load_alloc_vector(bcdiag, dir_mats // 'temp/bcdiag.dat')
	     rhs = rhs * bcdiag
	!call save_vector(bcdiag, dir_cols // 'bcdiag' //'.dat')
	     deallocate(bcdiag)

	!write(*, '("DONE")')
	!==========================================================================================
	!     Save the final version of the matrix
	if (P%messages>=3) then
	      call csr2coo(mat_csr_cmplx, mat_coo_cmplx, P%lib)
	      call save_sparse(mat_coo_cmplx, dir_mats // 'mat_' // cpt //'.dat')
	      call dealloc_sparse(mat_coo_cmplx)
	endif

	if (P%messages >= 1) then
		call CPU_Time(T1_cpt)
		call system_clock ( wall_t1_cpt, clock_rate, clock_max )
	endif

!********************************************************************
	if (P%solver .eq. 'umfpack') then
	!	USE UMFPACK TO SOLVE THE SYSTEM

	      ! convert from csr to csc for umfpack solver (csc is convert to 0-based array)
	        call csr2csc(mat_csr_cmplx, mat_csc_cmplx, P%lib)
	        call dealloc_sparse(mat_csr_cmplx)
	!								   messages,filesave,loadsym,savenum
	        call umfpack_unsym(mat_csc_cmplx%ni, mat_csc_cmplx, rhs, uvh, .false.,1,.false., .false.)

	        call dealloc_sparse(mat_csc_cmplx)
	        deallocate(rhs)
	elseif (P%solver .eq. 'pardiso') then
	!	USE PARDISO TO SOLVE THE SYSTEM
		  if ((P%pardiso_iterative == 1).and.(itn > 1)) CGS_ready = ( max( delta(itn-1, 7, ccpt), delta(itn-1, 8, ccpt) ) < 0.1)
	      call pardiso_unsym(mat_csr_cmplx%ni, mat_csr_cmplx, rhs, uvh, P, itn*ccpt, CGS_ready, dir_mats // 'temp/') ! directory where description file (dir_global) and factors will be saved

!	      call dealloc_sparse(mat_csr_cmplx)
!	      deallocate(rhs)
	elseif (P%solver .eq. 'gmres') then
		  call mkl_gmres(mat_csr_cmplx%ni, mat_csr_cmplx, rhs, uvh, P, itn, cpt, dir_sols)
		  ! mat_csr_cmplx and rhs are DEALLOCATED inside mkl_gmres(...)
	else
	      print *, "You must choose a valid solver: pardiso or umfpack"
	      stop
	end if
!********************************************************************
		if (P%messages >= 1) then
			call CPU_Time(T2_cpt)
			call system_clock ( wall_t2_cpt, clock_rate, clock_max )
			call disp ('done (CPU: ' // trim(ADJUSTL(conv_secs(T2_cpt-T1_cpt))) //', Wall: '&
						//trim(ADJUSTL(conv_secs( real(wall_t2_cpt-wall_t1_cpt)/real(clock_rate) )))//')')
		endif
	!==========================================================================================
	!    %================
	!    % update solution
	!    %================
		call save_solution(nu, nv, nh, itn, cpt, uvh, P%p_avg, dir_sols)
		deallocate(uvh)
	!==========================================================================================
	!****************** CLEAN UP Memory and Hard Drive (Pardiso) ********
!	if  (P%solver .eq. 'pardiso') then
!		call pardiso_clean(mat_csr_cmplx%ni, P,  dir_mats // 'temp/')! directory where description file (dir_global) and factors will be saved
!	end if

	enddo
!===============================================================================================
	!  %==============================================
	!  % examine height differences, i.e. convergence
	!  %==============================================
	call calc_convergence(nu, nv, nh, itn, P, cpts, delta, dir_cols, dir_sols)

	  if (itn > 1) then
	    if ( (P%cvg_scheme == 1) .and. (maxval( delta(itn,4,:) ) < P%cvg_dhbar) .and. (maxval( delta(itn,3,:) ) < P%cvg_dhmax) ) then
	        cvgd = 1
	    elseif ((P%cvg_scheme == 2).and.(maxval( delta(itn,6,:) ) < P%cvg_dhbar).and.(maxval(delta(itn,5,:))<P%cvg_dhmax)) then
	        cvgd = 1
	    endif
	  endif
!==========================================================================================
	!    %================================
	!    % calculating updated h_SAL terms
	!    %================================
	if (P%sal_scheme >= 2) then
		call calc_sal(ncpts, cpts, P, GD, dir_grid, dir_cols, dir_sols, beta0)
	endif

!==========================================================================================
	!    %====================================
	!    % calculate updated nonlinear terms:
	!    %====================================
	call calc_dhat(cpts, nu, nv, P, dir_cols, dir_mats, dir_sols)

	if (P%messages >= 2 ) then
		do ccpt = 1, ncpts
			cpt=cpts(2*ccpt-1:2*ccpt)
		!	Print the integrals on every iteration
			call baro_domain_integrals(cpt, P, GD%nu, GD%nv, GD%nh, dir_grid, dir_cols, dir_mats, dir_sols, di)
			call show_domain_integrals(cpt, di)
		enddo
	call disp('---------------------------------------------------------------------------')
	endif
	!==========================================================================================

	call CPU_Time(T2_itn)
	call system_clock ( wall_t2_itn, clock_rate, clock_max )
	call disp ('===>>> Total time spent on Iteration ' // tostring(itn) //': ' &
	            // trim(conv_secs(T2_itn-T1_itn)) // ' CPU, ' &
	            //trim(conv_secs( real(wall_t2_itn-wall_t1_itn)/real(clock_rate) ))//'  Wall')
	!==========================================================================================

	!      %==============================
	!      % plot global tides and errors
	!      %==============================

		if  ( P%graphics == 1 ) then
	!		Use MATLAB functions to export the data in M-files into binaries
			call system('xterm -e matlab22 -nosplash -nodesktop -logfile remoteAutocode.log -r ' // &
			'"addpath(genpath(''' // matlab_dir // ''')); baro_v1_plot_conv(''' // &
			trim(save_gridid)//''', '// tostring(itn) // ', ''' // cpts // '''); waitforbuttonpress; exit;" ')
	!		print *, ' matlab22 -nosplash -nodesktop -logfile remoteAutocode.log -r ' // &
	!		'"addpath(genpath(''' // matlab_dir // ''')); save_grid_binary(''2012_06_13__16_54'')" '

	!		if the matlab script failed

	!		if  (.not. dir_e ) then
	!			write(*,'(a)') 'Directory /'//trim(P%gridid)//' for the previously generated grid files doesn''t exist.'
	!			stop
	!		end if
		end if
	!==========================================================================================

	itn = itn + 1

enddo

!=========================================================================================
!=================================DONE====================================================
!=========================================================================================
write(*, '("========================================================================")')
write(*, '("========================================================================")')


!    %==============================
!    % plot global tides and errors
!    %==============================
!
!do ccpt = 1, ncpts
!
!     cpt=cpts(2*ccpt-1:2*ccpt)
!
!     call load_alloc_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
!     call load_alloc_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
!     call load_alloc_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat')
!
!     write(*,*) cpt, " max tide height: ", maxval(abs(h))
!
!         call load_alloc_matrix(hp, dir_cols // 'hp.dat')
!         call load_alloc_matrix(up, dir_cols // 'up.dat')
!         call load_alloc_matrix(vp, dir_cols // 'vp.dat')
!
!		call save_sol4plot_cmplx(nph, nth, nh, hp, h, dir_sols // 'sol_h_'//cpt//'.dat')
!		call save_sol4plot_cmplx(nph, nth, nu, up, u, dir_sols // 'sol_u_'//cpt//'.dat')
!		call save_sol4plot_cmplx(nph, nth, nv, vp, v, dir_sols // 'sol_v_'//cpt//'.dat')
!
!		deallocate(hp, up, vp)
!
!enddo

end subroutine baro_solver
!==========================================================================================

!==========================================================================================
subroutine baro_solver_linear(cpts, P, GD, dir_grid, dir_cols, dir_mats, dir_sols)

    implicit none
!%   Calculates a barotropic tide (linear bottom friction). Uses u/v/h + volume transport formulation
!%
!%   H_u : H on the u-grid points
!%   H_v : H on the v-grid points
!%   H_h : H on the h-grid points
!%   ta_u : tau on the u-grid points
!%   ta_v : tau on the v-grid points
!%   ta_h : tau on the h-grid points

     integer :: nph, nth
     integer :: nu, nv, nh, np
     real(wp):: latP, lonP
     real(wp):: g, re!, cdg, Q, beta
     character(len = *) :: dir_grid, dir_cols, dir_mats, dir_sols
     type (tide_params) :: pars
     type(params)	:: P
     type(grid_dims):: GD

     integer ::      j, istat
     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
     integer :: wall_t1, wall_t2, clock_rate, clock_max
     integer, allocatable :: hp(:, :),up(:, :),vp(:, :)

     character(len=*) :: cpts
     character(len=2) :: cpt
     integer :: ncpts, ccpt
     integer :: iu(GD%nu), iv(GD%nv), ih(GD%nh)

     type (triplet) :: mat
     type (triplet_cmplx) :: mat_coo_cmplx
     type (csr) :: mat_csr
     type (csr_cmplx) :: mat_csr_cmplx
     type (csc_cmplx) :: mat_csc_cmplx

	 type (csr) :: tmp_csr
	 type (csr_cmplx) :: tmp_csr_cmplx
	 type (triplet) :: tempmat
     integer ::      nmat, filled
	 type (triplet) :: DUU, DVU, DUV, DVV
     complex(cwp), allocatable :: Qu(:), Qv(:)!, tmp(:), Du(:), Dv(:)
	 real(wp), allocatable :: KU(:), KV(:)
     integer, allocatable, dimension(:) :: ufactor,vfactor

     type (csr_cmplx) :: speye
     integer, allocatable :: bcdiag(:)
     complex(cwp), allocatable :: rhs(:)

     complex(cwp), allocatable :: uvh(:), u(:), v(:), h(:)

!==========================================================================================
! SHORTCUTS
!==========================================================================================
latP = P%latP
lonP = P%lonP
g = P%g
re = P%re

nph = GD%nph
nth = GD%nta
nu = GD%nu
nv = GD%nv
nh = GD%nh
np = GD%np
!==========================================================================================
call tostring_set(rfmt='F12.1')
write(*, '("====================================")')
write(*, '("Welcome to baro_solver_linear")')
write(*, '("====================================")')


ncpts=len(cpts)/2

     iu = (/ ( (j), j=1,nu ) /)
     iv = (/ ( (nu + j), j=1,nv ) /)
     ih = (/ ( (nu + nv + j), j=1,nh ) /)

!%==========================
!% solve for each component
!%==========================
do ccpt = 1, ncpts

     cpt=cpts(2*ccpt-1:2*ccpt)

     pars = get_pars(cpt)

!%==================================================
!% Load the sparse matrix generated by baro_uvhT_mat
!%==================================================
     call load_alloc_sparse(mat, dir_mats // 'temp/mat_init.dat')
     call load_alloc_vector(bcdiag, dir_mats // 'temp/bcdiag.dat')

! Convert mat from COO to CSR:
     call coo2csr(mat, mat_csr, P%lib)
     call dealloc_sparse(mat)

!  now transform into cmplx matrix
     call init_sparse_vals(mat_csr_cmplx,mat_csr%indi,mat_csr%indj,mat_csr%vals, &
                                mat_csr%ni, mat_csr%nj, mat_csr%nz)
      call dealloc_sparse(mat_csr)
!%==================================================

      call disp('Solving ' // tostring(mat_csr_cmplx%ni) //' x '// tostring(mat_csr_cmplx%nj)//' system with '// &
                 tostring(mat_csr_cmplx%nz)// ' entries for:' )
!==========================================================================================
!    %========================
!    % add internal tide drag ! should be incorporated into write_baro_mats
!    %========================
    if (P%itd_scheme == 1) then

		call load_alloc_vector(KU, dir_cols // 'KU.dat')
		call load_alloc_vector(KV, dir_cols // 'KV.dat')
		call load_alloc_sparse(DUU, dir_mats // 'DUU.dat')
		call load_alloc_sparse(DUV, dir_mats // 'DUV.dat')
		call load_alloc_sparse(DVU, dir_mats // 'DVU.dat')
		call load_alloc_sparse(DVV, dir_mats // 'DVV.dat')
		allocate(Qu(nu), Qv(nv), stat = istat)

		call calc_itd_prefactor(ufactor, vfactor, cpt, P, nu, nv, dir_cols, dir_grid)
		!***********************************************
		!	Add terms into the matrix
		!***********************************************
		filled = 0
		nmat = DUU%nz + DUV%nz + DVU%nz + DVV%nz
		call init_sparse_0(tempmat, np,np, nmat)

		! 1) DUU
		call left_mult_diag(DUU, ufactor*KU/pars%omega0, nu)
		call concatenate_sparse(tempmat, DUU, 0, 0, filled)
		! 2) DVU
		call left_mult_diag(DVU, ufactor*KU/pars%omega0, nu)
		call concatenate_sparse(tempmat, DVU, 0, nu, filled)
		! 3) DUV
		call left_mult_diag(DUV, vfactor*KV/pars%omega0, nv)
		call concatenate_sparse(tempmat, DUV, nu, 0, filled)
		! 4) DVV
		call left_mult_diag(DVV, vfactor*KV/pars%omega0, nv)
		call concatenate_sparse(tempmat, DVV, nu, nu, filled)

!		Convert the tempmat from COO to CSR:
	    call coo2csr(tempmat, tmp_csr, P%lib)
	    call init_sparse_csr_cmplx(tmp_csr_cmplx,tmp_csr%indi,tmp_csr%indj,tmp_csr%vals, &
			                           tmp_csr%ni, tmp_csr%nj, tmp_csr%nz)
!	Save for testing
!		call save_vector( ufactor*KU/pars%omega0, dir_cols // 'tmp_u.dat')
!		call save_vector( vfactor*KV/pars%omega0, dir_cols // 'tmp_v.dat')
!     Save the final version of the matrix
!		call csr2coo(tmp_csr_cmplx, mat_coo_cmplx, P%lib)
!		call save_sparse(mat_coo_cmplx, dir_mats // 'mat1_' // cpt //'.dat')
!		call dealloc_sparse(mat_coo_cmplx)

	    call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), tmp_csr_cmplx, mkl)
	    call dealloc_sparse(tmp_csr)
	    call dealloc_sparse(tmp_csr_cmplx)

		call deallocate_sparse(DUU)
		call deallocate_sparse(DUV)
		call deallocate_sparse(DVU)
		call deallocate_sparse(DVV)
		deallocate(KU, KV, ufactor, vfactor)

	end if
!==========================================================================================

!     Save the final version of the matrix
!      call csr2coo(mat_csr_cmplx, mat_coo_cmplx, P%lib)
!      call save_sparse(mat_coo_cmplx, dir_mats // 'mat_' // cpt //'.dat')
!      call dealloc_sparse(mat_coo_cmplx)

!%=======================
!% calculate rhs forcing
!%=======================
     call baro_rhs(rhs, cpt, P, GD, dir_cols, dir_grid, dir_mats)
     call save_vector(rhs, dir_mats // 'rhs_' // cpt //'.dat')

!  %=============
!  % impose bcs:
!  %=============
write(*, '("Implementing boundary conditions...")', advance = 'no')
!mat=bcmat*mat; !fill the row of mat with zeroes for boundary variables
!rhs=bcmat*rhs; ! bcmat is already the way it should be!!!!!
     call dia_csr_mul(real(bcdiag,wp), mat_csr_cmplx, skit)
     rhs = rhs * bcdiag
     deallocate(bcdiag)
print *, "............. COMPLETE"

!==========================================================================================
!     mat = mat-i*omega0*speye(np)
      speye = speye_csr_cmplx(np, cmplx(0., -pars%omega0, kind=cwp) )
      call csr_csr_add(mat_csr_cmplx, cmplx(1., 0., kind=cwp), speye, mkl)

!  %=======
!  % solve
!  %=======

  write(*, '(" ", a, "... " )', advance='no') cpt
  call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

if (P%solver .eq. 'umfpack') then
!     USE UMFPACK TO SOLVE THE SYSTEM

      ! convert from csr to csc for umfpack solver (csc is convert to 0-based array)
        call csr2csc(mat_csr_cmplx, mat_csc_cmplx, P%lib)
        call dealloc_sparse(mat_csr_cmplx)
!                                           messages,filesave,loadsym,savenum
        call umfpack_unsym(mat_csc_cmplx%ni, mat_csc_cmplx, rhs, uvh, .false.,1,.false., .false.)
        call dealloc_sparse(mat_csc_cmplx)

elseif (P%solver .eq. 'pardiso') then
!     USE PARDISO TO SOLVE THE SYSTEM

	  call pardiso_unsym(mat_csr_cmplx%ni, mat_csr_cmplx, rhs, uvh, P, ccpt, .false., dir_mats // 'temp/') ! dir where factors can be saved
      call dealloc_sparse(mat_csr_cmplx)

else
      print *, "You must choose a valid solver: pardiso or umfpack"
      stop
end if

  call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )
  write(*, '("done in ", f5.1, "s.")') T2-T1

!  %=================
!  % output solution
!  %=================
    allocate(u(nu), v(nv), h(nh), stat = istat)
    u = uvh(iu)
    v = uvh(iv)
    h = uvh(ih)
    deallocate(uvh)

     call save_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')

     call save_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
     call save_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat')

	deallocate(u, v, h)

enddo

!=========================================================================================
!=================================DONE====================================================
!=========================================================================================
write(*, '("========================================================================")')
write(*, '("========================================================================")')

do ccpt = 1, ncpts

     cpt=cpts(2*ccpt-1:2*ccpt)

     call load_alloc_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
     call load_alloc_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
     call load_alloc_vector(h, dir_sols // cpt // '_h' // '_m0' // '.dat')

     write(*,*) cpt, " max tide height: ", maxval(abs(h))

!  if flags.graphics == 1
!
!    %==============================
!    % plot global tides and errors
!    %==============================

     if (ccpt == 1) then

         call load_alloc_matrix(hp, dir_cols // 'hp.dat')
         call load_alloc_matrix(up, dir_cols // 'up.dat')
         call load_alloc_matrix(vp, dir_cols // 'vp.dat')

		call save_sol4plot_cmplx(nph, nth, nh, hp, h, dir_sols // 'sol_h_'//cpt//'.dat')
		call save_sol4plot_cmplx(nph, nth, nu, up, u, dir_sols // 'sol_u_'//cpt//'.dat')
		call save_sol4plot_cmplx(nph, nth, nv, vp, v, dir_sols // 'sol_v_'//cpt//'.dat')

		deallocate(hp, up, vp)

     endif

	deallocate(u,v,h)

enddo

end subroutine baro_solver_linear
!==========================================================================================

!==========================================================================================
subroutine baro_rhs(rhs, cpt, P, GD, dir_cols, dir_grid, dir_mats)

    implicit none
!
!% calculates the complex barotropic tidal forcing
!% for the given tidal component.

     character(len=2), intent(in)   :: cpt
     type(params), intent(in)		:: P
     type(grid_dims), intent(in)	:: GD

     integer  :: nu, nv, nh, coor
     real(wp) :: latP, lonP

     type (triplet) :: GHU, GHV

     integer ::      j
     integer ::      istatus
     integer, allocatable :: 	 iu(:), iv(:)
!     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
!     integer :: wall_t1, wall_t2, clock_rate, clock_max

     complex(cwp), allocatable, dimension(:) :: rhs, heq
!     real(wp), allocatable, dimension(:) :: tmp2, ph_h, ph_u

     character(len = *) :: dir_grid, dir_cols, dir_mats

coor = P%coor
latP =  P%latP
lonP = P%lonP

nu = GD%nu
nv = GD%nv
nh = GD%nh


	allocate(iu(nu), iv(nv), stat=istatus)
     iu = (/ ( (j), j=1,nu ) /)
     iv = (/ ( (j), j=1,nv ) /)
!%==============================================================
!% Load the necessary matrices:
!% GHU, GHV
!%==============================================================
     call load_alloc_sparse(GHU, dir_mats // 'GHU.dat')
     call load_alloc_sparse(GHV, dir_mats // 'GHV.dat')
!*************************************************************************
     allocate(rhs(nu+nv+nh), stat = istatus)
     rhs = 0.

     call calc_heq_h(heq, cpt, nh, latP, lonP, coor, dir_cols, dir_grid)
!     call save_vector(heq, dir_cols // 'heq.dat')

     call coo_vec_mul(GHU, heq, P%lib, rhs(1:nu))
     call dealloc_sparse(GHU)

     call coo_vec_mul(GHV, heq, P%lib, rhs(nu+1:nu+nv))
     call dealloc_sparse(GHV)

!     call save_vector(rhs, dir_cols // 'rhs.dat')

end subroutine baro_rhs

!**********************************************************************************************************
!**********************************************************************************************************

end module baro_solver_mod



