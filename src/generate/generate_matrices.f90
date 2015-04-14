module generate_matrices

     use precisions, only: dp, wp, cwp
     use dispmodule
     use my_trigs
     use my_sparse
     use my_sparse_aggregate
     use save_load
     use generate_grid
     use itd
	 use control

     contains

!    %========================================================
!    % Generate aggregate matrices
!    %========================================================
!==========================================================================================
subroutine write_baro_mats(N_data_dir, dir_cols, dir_grid, dir_mats, cpts, P, nu, nv, nh)

     character(len=*) :: dir_cols, dir_grid, dir_mats
     character(len = *) :: N_data_dir
     type(params) :: P
     real(wp)	  :: cdg, udg, Q
     integer :: nu, nv, nh, np, cooruv, nph, nth
     real(wp)	:: re

     integer, allocatable :: up(:, :), vp(:, :)
     integer, target  :: iu(nu), iv(nv)!, ih(nh), ip(np)
     real(wp), allocatable ::  ta_u(:), ta_v(:), ph_vg(:), ta_ug(:)!, ta_h(:)

!     real    ::      T1, T2
!     integer :: wall_t1, wall_t2, clock_rate, clock_max

     integer ::      j
     integer ::      istat!, statusj, statusval

     character(len=*) :: cpts
     character(len=2) :: cpt
     integer :: ncpts, ccpt
	 type (tide_params) :: pars

     real(wp), allocatable :: Qu(:), Qv(:)
     real(wp), allocatable :: H_u(:), H_v(:), H_h(:), KU(:), KV(:)
     real(wp), allocatable :: dHdph(:), dHdta(:), N_u(:), N_v(:)

	 type (triplet)			:: tempmat
     integer				:: nmat, filled ! points how many elements out of nmat are already filled
	 integer, allocatable	:: ufactor(:), vfactor(:)
	 type (csr) :: tmp_csr

     type (triplet) :: GHU, GHV, JUH, JVH, DUU, DVU, DUV, DVV
     type (triplet) :: u2vf, v2uf, u2v, v2u, h2uddph, h2vddta

     type (triplet) :: mat
	 type (csr) :: mat_csr, mat_csr_cpt
     integer, pointer :: bcdiag(:)

!%==============================
write(*, '("Preparing matrices: ")', advance = 'no')

!%=============================================
! shortcuts
!%=============================================
!g = P%g
re = P%re
!%====================================================================
!% differentiate between lat-lon/MERC, velocity/transport formulations
!%====================================================================
!	Only Transport-uvh formulation so far

if (P%coor == 1) then ! MERCATOR
!  if (P%uvform == 1) then
!    cooruv=11 ! VELOCITY
!  elseif (Puvform == 2)
    cooruv=12 ! TRANSPORT
!  endif
elseif (P%coor == 2) then ! LAT/LON
!  if (P%uvform == 1) then
!    cooruv=21 ! VELOCITY
!  elseif (Puvform == 2)
    cooruv=22 ! TRANSPORT
!  endif
endif
!%=============================================

!    %==========================
!    % set up the basic matrices
!    %==========================
      call baro_uvhT_mat(nu, nv, nh, P, nmat, dir_cols, dir_mats)

!    %=====================================
!    % rotation, c.o.mass and RHS (forcing)
!    %=====================================

     iu = (/ ( (j), j=1,nu ) /)
     iv = (/ ( (j), j=1,nv ) /)
!     ih = (/ ( (j), j=1,nh ) /)
!     ip = (/ ( (j), j=1,np ) /)

!%=============================================
!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!
!	Target:
!	nmat = v2uf%nz + u2vf%nz + h2uddph%nz + h2vddta%nz + u2hddph%nz + v2hddta%nz + nu + nv
!	After calling routine baro_uvhT_mat:
!	nmat = h2uddph%nz + h2vddta%nz + u2hddph%nz + v2hddta%nz+ nu + nv
!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!
	call load_alloc_sparse(v2uf, dir_mats // 'v2uf.dat')
	call load_alloc_sparse(u2vf, dir_mats // 'u2vf.dat')

	nmat = nmat + v2uf%nz + u2vf%nz
	np = nu + nv + nh
	filled = 0
	call init_sparse_0(mat, np,np, nmat)

!%===============
!	Add rotation
!%===============
!	1) v2uf
	v2uf%vals = - v2uf%vals ! multiply by -1
	call concatenate_sparse(mat, v2uf, 0, nu, filled)
	call dealloc_sparse(v2uf)

	write(*, '("v2uf, ")', advance = 'no')
!	2) u2vf
    call concatenate_sparse(mat, u2vf, nu, 0, filled)
	call dealloc_sparse(u2vf)

	write(*, '("u2vf, ")', advance = 'no')

!%==============================
!	Add already generated terms
!%==============================
!	Gravity (linear SAL is included if sal_scheme==1):
	call load_alloc_sparse(GHU, dir_mats // 'GHU.dat')
	call load_alloc_sparse(GHV, dir_mats // 'GHV.dat')

	if (P%sal_scheme == 1) then
		GHU%vals = (1 - P%beta0)*GHU%vals
		GHV%vals = (1 - P%beta0)*GHV%vals
	end if

	call concatenate_sparse(mat, GHU, 0, nu+nv, filled)
	call concatenate_sparse(mat, GHV, nu, nu+nv, filled)

	call dealloc_sparse(GHU)
	call dealloc_sparse(GHV)

!%===================================
!	Advection gradients (c.o. mass)
!%===================================
	call load_alloc_sparse(JUH, dir_mats // 'JUH.dat')
	call load_alloc_sparse(JVH, dir_mats // 'JVH.dat')

	call concatenate_sparse(mat, JUH, nu+nv, 0, filled)
	call concatenate_sparse(mat, JVH, nu+nv, nu, filled)

	call dealloc_sparse(JUH)
	call dealloc_sparse(JVH)

!%===================================
!	Linear bottom friction
!%===================================
	call load_alloc_vector(KU, dir_cols // 'KU.dat')
	call load_alloc_vector(KV, dir_cols // 'KV.dat')

	call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
	call load_alloc_vector(H_v, dir_cols // 'H_v.dat')

!	TRANSPORT FORMULATION
	allocate(Qu(nu), Qv(nv), stat = istat)
		cdg = P%cd
	if (P%fr_scheme == 1) then
		udg = P%ubar
		Qu = cdg*udg
		Qv = cdg*udg
	elseif (P%fr_scheme == 2) then
		Q = P%Qbar
		Qu = cdg*Q/H_u
		Qv = cdg*Q/H_v
	else
		Qu = 0
		Qv = 0
	end if

	if (P%coor == 1) then
		call load_alloc_vector(ta_u, dir_cols // 'ta_u.dat')
		call load_alloc_vector(ta_v, dir_cols // 'ta_v.dat')

		Qu=Qu * cosh(ta_u)
		Qv=Qv * cosh(ta_v)

		deallocate(ta_u, ta_v)
	endif
!      if (uvform == 2) then
        Qu=Qu/H_u
        Qv=Qv/H_v
!      endif


	call init_sparse(tempmat, iu, iu, KU*Qu, nu, nu, nu)
	call concatenate_sparse(mat, tempmat, 0, 0, filled)
	call dealloc_sparse(tempmat)

	call init_sparse(tempmat, iv, iv, KV*Qv, nv,nv, nv)
	call concatenate_sparse(mat, tempmat, nu, nu, filled)
	call dealloc_sparse(tempmat)
	write(*, '("Linear BBL, ")', advance = 'no')

	deallocate(Qu, Qv, H_u, H_v, KU, KV)

! Convert mat from COO to CSR:
     call coo2csr(mat, mat_csr, P%lib)
     call save_sparse(mat, dir_mats // 'temp/mat_init.dat')
     call dealloc_sparse(mat)

if (P%itd_scheme /= 1) then
     call save_sparse(mat_csr, dir_mats // 'temp/' // 'mat_csr_init.dat')
     call dealloc_sparse(mat_csr)
endif
!%======================
!	Internal tide drag	! Done separately bc it requires adding two filled sparse matrices
!%======================
if (P%itd_scheme == 1) then
	call load_alloc_sparse(h2uddph, dir_mats // 'h2uddph.dat')
	call load_alloc_sparse(h2vddta, dir_mats // 'h2vddta.dat')
	call load_alloc_sparse(v2u, dir_mats // 'v2u.dat')
	call load_alloc_sparse(u2v, dir_mats // 'u2v.dat')

	if (P%smooth_type < 1) then
		call load_alloc_vector(H_h, dir_cols // 'H_h.dat')
	else
		call load_alloc_vector(H_h, dir_cols // 'H_sht_h.dat')
	endif

	call load_alloc_vector(ta_u, dir_cols // 'ta_u.dat')
	call load_alloc_vector(ta_v, dir_cols // 'ta_v.dat')

    allocate(dHdph(nu),dHdta(nv), stat = istat)

	call coo_vec_mul(h2uddph, H_h, P%lib, dHdph)
	call coo_vec_mul(h2vddta, H_h, P%lib, dHdta)

	call dealloc_sparse(h2uddph)
	call dealloc_sparse(h2vddta)
	deallocate(H_h)

	call calc_N_on_grid(N_u, N_v, P, nu, nv, nh, N_data_dir, dir_cols, dir_grid, dir_mats)
!		call save_vector( N_u, dir_cols // 'N_u.dat')
!		call save_vector( N_v, dir_cols // 'N_v.dat')

	allocate(Qu(nu), Qv(nv), stat = istat)

	!***********************************************************************************
	select case (cooruv)
	!***********************************************************************************
		case (12) ! MERCATOR, TRANSPORT
			Qu = P%itd_coeff*N_u**2 * (dHdph*cosh(ta_u)/re)**2 * cosh(ta_u)
			call init_sparse(DUU, iu, iu, Qu, nu, nu, nu)

		    call right_mult_diag(v2u, dHdta*cosh(ta_v)**2/re, nv)
			Qu = P%itd_coeff*N_u**2 * (dHdph*cosh(ta_u)/re)
		    call left_mult_diag(v2u, Qu, nu)
			call init_sparse(DVU, v2u%indi, v2u%indj, v2u%vals, v2u%ni, v2u%nj, v2u%nz)

		    call right_mult_diag(u2v, dHdph*cosh(ta_u)**2/re, nu)
			Qv = P%itd_coeff*N_v**2 * (dHdta*cosh(ta_v)/re)
		    call left_mult_diag(u2v, Qv, nv)
			call init_sparse(DUV, u2v%indi, u2v%indj, u2v%vals, u2v%ni, u2v%nj, u2v%nz)

			Qv = P%itd_coeff*N_v**2 * (dHdta*cosh(ta_v)/re)**2 * cosh(ta_v)
			call init_sparse(DVV, iv, iv, Qv, nv, nv, nv)

		case (22) ! LAT-LON, TRANSPORT
			Qu = P%itd_coeff*N_u**2 * (dHdph/(re*cos(ta_u)))**2
			call init_sparse(DUU, iu, iu, Qu, nu, nu, nu)

		    call right_mult_diag(v2u, dHdta/re, nv)
			Qu = P%itd_coeff*N_u**2 * (dHdph/(re*cos(ta_u)))
		    call left_mult_diag(v2u, Qu, nu)
			call init_sparse(DVU, v2u%indi, v2u%indj, v2u%vals, v2u%ni, v2u%nj, v2u%nz)

		    call right_mult_diag(u2v, dHdph/(re*cos(ta_u)), nu)
			Qv = P%itd_coeff*N_v**2 * (dHdta/re)
		    call left_mult_diag(u2v, Qv, nv)
			call init_sparse(DUV, u2v%indi, u2v%indj, u2v%vals, u2v%ni, u2v%nj, u2v%nz)

			Qv = P%itd_coeff*N_v**2 * (dHdta/re)**2
			call init_sparse(DVV, iv, iv, Qv, nv, nv, nv)

	end select

	! Save dHdph, dHdta to compare smoothness
!        call load_alloc_matrix(up, dir_cols // 'up.dat')
!        call load_alloc_matrix(vp, dir_cols // 'vp.dat')
!		call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat', nph)
!		call load_alloc_vector(ta_ug, dir_grid // 'ta_ug.dat', nth)
!	select case (cooruv)
!	!***********************************************************************************
!		case (12) ! MERCATOR, TRANSPORT
!			call save_sol4plot(nph, nth, nu, up, dHdph*cosh(ta_u)/re, dir_mats // 'dHdph.dat')
!			call save_sol4plot(nph, nth, nv, vp, dHdta*cosh(ta_v)/re, dir_mats // 'dHdta.dat')
!		case (22) ! LAT-LON, TRANSPORT
!			call save_sol4plot(nph, nth, nu, up, dHdph/(re*cos(ta_u)), dir_mats // 'dHdph.dat')
!			call save_sol4plot(nph, nth, nv, vp, dHdta/re, dir_mats // 'dHdta.dat')
!	end select
!
!		deallocate(ph_vg, ta_ug, up, vp)


	deallocate(Qu, Qv, dHdph, dHdta, N_u, N_v, ta_u, ta_v)
	call deallocate_sparse(v2u)
	call deallocate_sparse(u2v)

	call save_sparse(DUU, dir_mats // 'DUU.dat')
	call save_sparse(DUV, dir_mats // 'DUV.dat')
	call save_sparse(DVU, dir_mats // 'DVU.dat')
	call save_sparse(DVV, dir_mats // 'DVV.dat')

	call deallocate_sparse(DUU)
	call deallocate_sparse(DUV)
	call deallocate_sparse(DVU)
	call deallocate_sparse(DVV)

!==========================================================================================

		call load_alloc_vector(KU, dir_cols // 'KU.dat')
		call load_alloc_vector(KV, dir_cols // 'KV.dat')

	ncpts=len(cpts)/2

	do ccpt = 1, ncpts

		cpt=cpts(2*ccpt-1:2*ccpt)
		pars = get_pars(cpt)

		call calc_itd_prefactor(ufactor, vfactor, cpt, P, nu, nv, dir_cols, dir_grid)
		!***********************************************
		!	Add terms into the matrix
		!***********************************************
		call load_alloc_sparse(DUU, dir_mats // 'DUU.dat')
		call load_alloc_sparse(DUV, dir_mats // 'DUV.dat')
		call load_alloc_sparse(DVU, dir_mats // 'DVU.dat')
		call load_alloc_sparse(DVV, dir_mats // 'DVV.dat')

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

		deallocate(ufactor, vfactor)
		call deallocate_sparse(DUU)
		call deallocate_sparse(DUV)
		call deallocate_sparse(DVU)
		call deallocate_sparse(DVV)

		! Convert the tempmat from COO to CSR:
	    call coo2csr(tempmat, tmp_csr, P%lib)
	    call csr_csr_add(mat_csr, real(1., kind=wp), tmp_csr, mkl, mat_csr_cpt)
	    call dealloc_sparse(tmp_csr)
	    call dealloc_sparse(tempmat)

!     Save the COO version of the matrix
!		call csr2coo(mat_csr_cpt, mat, P%lib)
!		call save_sparse(mat, dir_mats // 'temp/' // cpt // '_mat_init.dat')
!		call dealloc_sparse(mat)

!%================================================
!% Store the initial matrix for each cpt in a file
!%================================================
     call save_sparse(mat_csr_cpt, dir_mats // 'temp/' // cpt // '_mat_csr_init.dat')
     call dealloc_sparse(mat_csr_cpt)

	end do

    call dealloc_sparse(mat_csr)
	deallocate(KU, KV)

end if
!==========================================================================================

	write(*, '("Linear ITD (1 matrix per constituent), ")', advance = 'no')
!==========================================================================================

!%==================================================
!% define the b.c. diagonal matrix:
!%==================================================
! load up, vp, hp
     call load_alloc_matrix(up, dir_cols // 'up.dat')
     call load_alloc_matrix(vp, dir_cols // 'vp.dat')

  allocate(bcdiag(np), stat = istat)
  bcdiag(:) = (/ 1-up(:,3), 1-vp(:,3), (1, j=1,nh) /)

write(*, '("b.c. ")', advance = 'no')

      call save_vector(bcdiag, dir_mats // 'temp/bcdiag.dat')
      deallocate(up, vp, bcdiag)
!==========================================================================================
write(*, '("DONE")')


end subroutine

!**********************************************************************************************************
!**********************************************************************************************************

subroutine baro_uvhT_mat(nu, nv, nh, P, nmat, dir_cols, dir_mats)

    implicit none

!% Produces a linear inversion/propagator matrix for a
!% column uvh using the TRANSPORT VELOCITY formulation.
!% Includes LINEAR friction
!% DOES NOT include scalar SAL
!% No frequency dependent d/dt or internal tide drag terms.
     character(len=*) :: dir_cols, dir_mats
!     character(len = *) :: N_data_dir
     type(params) :: P
     integer, intent(in) :: nu, nv, nh
     integer, intent(out) :: nmat

	integer		:: cooruv !differentiate between lat-lon/MERC, velocity/transport formulations
!     integer, allocatable :: up(:, :), vp(:, :)
     real(wp)			  :: g, re
     real(wp), allocatable ::  ta_u(:), ta_v(:), ta_h(:), ones(:)
     real(wp), allocatable :: H_u(:), H_v(:)!, H_h(:)
     integer	:: istat

     type (triplet) :: h2uddph, h2vddta, u2hddph, v2hddta
!     type (triplet) :: GHU, GHV, JUH, JVH, DUU, DVU, DUV, DVV
!     real(wp), allocatable :: KU(:), KV(:)

!%=============================================
! shortcuts
!%=============================================
g = P%g
re = P%re
!%====================================================================
!% differentiate between lat-lon/MERC, velocity/transport formulations
!%====================================================================
!	Only Transport-uvh formulation so far

if (P%coor == 1) then ! MERCATOR
!  if (P%uvform == 1) then
!    cooruv=11 ! VELOCITY
!  elseif (Puvform == 2)
    cooruv=12 ! TRANSPORT
!  endif
elseif (P%coor == 2) then ! LAT/LON
!  if (P%uvform == 1) then
!    cooruv=21 ! VELOCITY
!  elseif (Puvform == 2)
    cooruv=22 ! TRANSPORT
!  endif
endif
!%=============================================



!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!
!	nmat = h2uddph%nz + h2vddta%nz + u2hddph%nz + v2hddta%nz
!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!
!	realized in steps 1-2 and 5-6
	nmat = nu + nv

	call load_alloc_vector(ta_u, dir_cols // 'ta_u.dat')
	call load_alloc_vector(ta_v, dir_cols // 'ta_v.dat')
	call load_alloc_vector(ta_h, dir_cols // 'ta_h.dat')

select case (cooruv)
!***********************************************************************************
	case (12) ! MERCATOR, TRANSPORT
!***********************************************************************************
	!%================================================
	!	Pressure gradients (gravity like potentials)
	!%================================================
	!	1) GHU
		call load_alloc_sparse(h2uddph, dir_mats // 'h2uddph.dat')
		call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
	    call left_mult_diag(h2uddph, (g/re)*H_u, nu)

		call save_sparse(h2uddph, dir_mats // 'GHU.dat')
		nmat = nmat + h2uddph%nz
		call dealloc_sparse(h2uddph)
		deallocate(H_u)

		write(*, '("GHU, ")', advance = 'no')

	!	2) GHV
		call load_alloc_sparse(h2vddta, dir_mats // 'h2vddta.dat')
		call load_alloc_vector(H_v, dir_cols // 'H_v.dat')
	    call left_mult_diag(h2vddta, (g/re)*H_v, nv)

		call save_sparse(h2vddta, dir_mats // 'GHV.dat')
		nmat = nmat + h2vddta%nz
		call dealloc_sparse(h2vddta)
		deallocate(H_v)

		write(*, '("GHV, ")', advance = 'no')

	!%================================================
	!	Metrics multipliers
	!%================================================
	!	3) KU
		call save_vector(1/cosh(ta_u), dir_cols // 'KU.dat')
		write(*, '("KU, ")', advance = 'no')
	!	4) KV
		call save_vector(1/cosh(ta_v), dir_cols // 'KV.dat')
		write(*, '("KV, ")', advance = 'no')

	!%================================================
	!	Advection gradients (c.o. mass)
	!%================================================
	!	5) JUH
		call load_alloc_sparse(u2hddph, dir_mats // 'u2hddph.dat')
	    call left_mult_diag(u2hddph, (1/re)*cosh(ta_h)**2, nh)

		call save_sparse(u2hddph, dir_mats // 'JUH.dat')
		nmat = nmat + u2hddph%nz
		call dealloc_sparse(u2hddph)

	!	6) JVH
		call load_alloc_sparse(v2hddta, dir_mats // 'v2hddta.dat')
		call left_mult_diag(v2hddta, (1/re)*cosh(ta_h)**2, nh)

		call save_sparse(v2hddta, dir_mats // 'JVH.dat')
		nmat = nmat + v2hddta%nz
		call dealloc_sparse(v2hddta)

!***********************************************************************************
	case (22) ! LAT-LON, TRANSPORT
!***********************************************************************************
	!%================================================
	!	Pressure gradients (gravity like potentials)
	!%================================================
	!	1) GHU
		call load_alloc_sparse(h2uddph, dir_mats // 'h2uddph.dat')
		call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
	    call left_mult_diag(h2uddph, (g/re)*H_u/cos(ta_u), nu)

		call save_sparse(h2uddph, dir_mats // 'GHU.dat')
		nmat = nmat + h2uddph%nz
		call dealloc_sparse(h2uddph)
		deallocate(H_u)

		write(*, '("GHU, ")', advance = 'no')

	!	2) GHV
		call load_alloc_sparse(h2vddta, dir_mats // 'h2vddta.dat')
		call load_alloc_vector(H_v, dir_cols // 'H_v.dat')
	    call left_mult_diag(h2vddta, (g/re)*H_v, nv)

		call save_sparse(h2vddta, dir_mats // 'GHV.dat')
		nmat = nmat + h2vddta%nz
		call dealloc_sparse(h2vddta)
		deallocate(H_v)

		write(*, '("GHV, ")', advance = 'no')

	!%================================================
	!	Metrics multipliers
	!%================================================
	!	3) KU
		allocate(ones(nu), stat = istat)
		ones = 1
		call save_vector(ones, dir_cols // 'KU.dat')
		deallocate(ones)
		write(*, '("KU, ")', advance = 'no')
	!	4) KV
		allocate(ones(nv), stat = istat)
		ones = 1
		call save_vector(ones, dir_cols // 'KV.dat')
		deallocate(ones)
		write(*, '("KV, ")', advance = 'no')

	!%================================================
	!	Advection gradients (c.o. mass)
	!%================================================
	!	5) JUH
		call load_alloc_sparse(u2hddph, dir_mats // 'u2hddph.dat')
	    call left_mult_diag(u2hddph, (1/re)*1/cos(ta_h), nh)

		call save_sparse(u2hddph, dir_mats // 'JUH.dat')
		nmat = nmat + u2hddph%nz
		call dealloc_sparse(u2hddph)

	!	6) JVH
		call load_alloc_sparse(v2hddta, dir_mats // 'v2hddta.dat')
		call left_mult_diag(v2hddta, (1/re)*1/cos(ta_h), nh)
		call right_mult_diag(v2hddta, cos(ta_v), nv)

		call save_sparse(v2hddta, dir_mats // 'JVH.dat')
		nmat = nmat + v2hddta%nz
		call dealloc_sparse(v2hddta)

end select

	deallocate(ta_h, ta_u, ta_v)

end subroutine baro_uvhT_mat

!==========================================================================================
!==========================================================================================

subroutine generate_global_matrices(dir_grid, dir_cols, dir_mats, P, GD)

     character(len=*) :: dir_grid, dir_cols, dir_mats
     type(params) :: P
     type(grid_dims):: GD

     logical :: flag_itd_smooth, flag_baro_smooth

!     real(dp), allocatable :: H_hg(:, :), H_sht_hg(:, :)


!%=======================================
!% write H, ta, heq etc. on the columns:
!%=======================================
	 flag_itd_smooth = (P%itd_scheme>0).and.(P%smooth_type>0).and.(P%sht_smooth_H >= 0)
	 flag_baro_smooth = (P%smooth_type>0).and.(P%baro_on_smoothed == 1)
     call write_global_ccols_col(GD%nph, GD%nu, GD%nv, GD%nh, dir_grid, dir_cols, flag_itd_smooth, flag_baro_smooth)
!     deallocate(H_hg,H_sht_hg) ! don't need it anymore, free space
!%====================================
!% generate and write sparse matrices
!%====================================
     call write_global_cmats_col(GD%nph, GD%nu, GD%nv, GD%nh, P%latP, P%omega, dir_grid, dir_cols, dir_mats)

!     deallocate(ta_ug, ta_vg, ph_ug, ph_vg, H_hg)
!     deallocate(up, vp, hp, iu_ug, iv_vg, ih_hg)

     write(*, '("done.")')

end subroutine

!    %========================================================
!    % BACIS matrices
!    %========================================================
!==========================================================================================
subroutine write_global_cmats_col(nph, nu, nv, nh, latP, omega, dir_grid, dir_cols, dir_mats)
!%===============================================================
!% load the grid: H_hg, ta_ug, ta_vg, ph_ug, ph_vg, th_ug, th_vg,
!%                iu_ug, iv_vg, ih_hg
!%===============================================================
    implicit none

!% Generate some sparse matrices for operations on an Arakawa C grid.
!%
!%   h2uddph, h2u
!%   h2vddta, h2v
!%   u2hddph, u2h
!%   v2hddta, v2h
!%   u2v, u2vf, u2vfsp
!%   v2u, v2uf, v2ufsp

     integer                   ::     nph!, nth
     real(dp), allocatable     ::     th_vg(:), ta_vg(:), ph_ug(:)

     integer, allocatable :: up(:, :), vp(:, :), hp(:, :)
     integer, intent(in)  ::     nu, nv, nh!, np
     integer              :: nmax
     integer, allocatable :: iu_ug(:, :), iv_vg(:, :), ih_hg(:, :)
     real(dp), allocatable :: H_u(:), H_v(:)!, H_h(:)


     integer, pointer :: ivals(:), jvals(:), ispvals(:), jspvals(:)
     real(dp), pointer :: vals(:), vals1(:), fvals(:), fspvals(:)

     real(dp)  ::      dph, dta
     integer ::      istat, statusj, status1, status2
     integer ::      cu, cv, ch, cth, cph
     integer ::      cphl, cphr, ctha, cthb, nvs
     integer ::      hpl, hpr, hpa, hpb
     integer ::      upl, upr
     integer ::      vpa, vpb
     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
     integer :: wall_t1, wall_t2, clock_rate, clock_max

     real(dp), intent(in) :: omega, latP
     real(dp)			  :: th0

     character(len = *) :: dir_grid, dir_cols, dir_mats

!     real(dp), dimension(4)   :: thc, phc

     type (triplet_dp) :: h2uddph, h2u, h2vddta, h2v, u2hddph, u2h, v2hddta, v2h
     type (triplet_dp) :: u2v, u2vf, u2vfsp, v2u, v2uf, v2ufsp

!%=============================================
!% Load the necessary grid matrices and vectors
!%=============================================
! load up, vp, hp
     call load_alloc_matrix(up, dir_cols // 'up.dat')
     call load_alloc_matrix(vp, dir_cols // 'vp.dat')
     call load_alloc_matrix(hp, dir_cols // 'hp.dat')

     call load_alloc_vector(ph_ug, dir_grid // 'ph_ug.dat')
     call load_alloc_vector(ta_vg, dir_grid // 'ta_vg.dat') ! in Mercator coords it is ta
     														! in lat-lon coords it is th
     call load_alloc_vector(th_vg, dir_grid // 'th_vg.dat')
! load iu_ug, iv_vg, ih_hg
     call load_alloc_matrix(iu_ug, dir_grid // 'iu_ug.dat')
     call load_alloc_matrix(iv_vg, dir_grid // 'iv_vg.dat')
     call load_alloc_matrix(ih_hg, dir_grid // 'ih_hg.dat')
! Load columns H_u, H_v
     call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
     call load_alloc_vector(H_v, dir_cols // 'H_v.dat')
!%=================
!% make shortcuts:
!%=================

th0 = d2r(latP)
!ph0 = d2r(lonP)

!dph=2*pi/nph;    ! the distance between ph gridpoints
dph = ph_ug(2)-ph_ug(1)
!dta=dph;    ! the distance between ta gridpoints
dta = ta_vg(2)-ta_vg(1)

!     SECOND ORDER FINITE DIFFERENCE SCHEME

write(*, '("Making matrices:")')

nmax = max(nu, nv, nh)
    if (associated(ivals)) deallocate(ivals) ! 4*nmax for 4th order scheme
    allocate(ivals(2*nmax), stat = istat)
    if (associated(jvals)) deallocate(jvals)
    allocate(jvals(2*nmax), stat = statusj)
    if (associated(vals)) deallocate(vals)
    allocate(vals(2*nmax), stat = status1)
    if (associated(vals1)) deallocate(vals1)
    allocate(vals1(2*nmax), stat = status2)

    ivals = 0
    jvals = 0
    vals = 0
    vals1 = 0

    nvs = 0


write(*, '(" - d/dph and 1 for h-grid functions, evaluating on the u-grid ")', advance = 'no')

!%=============
!% initialize:
!%=============
!!!!!!!!!!!!!!       h2uddph, h2u
!%=============

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do cu = 1,nu

       cph=up(cu,1)
       cth=up(cu,2)

      cphr=cph
      cphl=1+modulo((cph-1)-1,nph)
       hpl=ih_hg(cphl,cth)
       hpr=ih_hg(cphr,cth)

       if ((hpl /= 0).and.(hpr /= 0)) then !% we are not at a Western/Eastern boundary

      !% 4th order FD
!              cphll=1+mod((cphl-1)-1,nph)
!              cphrr=1+mod((cphr+1)-1,nph);
!
!              hpll=ih_hg(cth,cphll); hprr=ih_hg(cth,cphrr);
!              if hpll > 0 & hprr > 0 & fdflag == 4
!                nvs=nvs+4;
!                ivals(nvs-3:nvs)=cu;
!                jvals(nvs-3:nvs)=[hpll hpl hpr hprr];
!                vals(nvs-3:nvs)=[1/24 -9/8 9/8 -1/24]/dph;
!                vals1(nvs-3:nvs)=[-1 9 9 -1]/16;
!              else
                nvs=nvs+2
                ivals(nvs-1:nvs)=cu
                jvals(nvs-1:nvs)=(/hpl, hpr/)
                vals(nvs-1:nvs)=(/-1., 1./)/dph
                vals1(nvs-1:nvs)=(/0.5, 0.5/)
       endif

    enddo

call init_sparse_dp(h2uddph,ivals,jvals,vals,nu,nh,nvs)
call init_sparse_dp(h2u,ivals,jvals,vals1,nu,nh,nvs)

call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%=============
!!!!!!!!!!!!!!       h2vddta, h2v
!%=============

write(*, '(" - d/dta and 1 for h grid functions, evaluating on the v-grid ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do cv = 1,nv

       cph=vp(cv,1)
       cth=vp(cv,2)

      ctha=cth
      cthb=cth-1
       hpa=ih_hg(cph,ctha)
       hpb=ih_hg(cph,cthb)

       if ((hpa /= 0).and.(hpb /= 0)) then !% we are not at a Northern/Southern boundary

      !% 4th order FD
!         cthaa=ctha+1; cthbb=cthb-1;
!         hpaa=ih_hg(cthaa,cph); hpbb=ih_hg(cthbb,cph);
!         if hpaa > 0 && hpbb > 0 && fdflag == 4
!           nvs=nvs+4;
!           ivals(nvs-3:nvs)=cv;
!           jvals(nvs-3:nvs)=[hpbb hpb hpa hpaa];
!           vals(nvs-3:nvs)=[1/24 -9/8 9/8 -1/24]/dta;
!           vals1(nvs-3:nvs)=[-1 9 9 -1]/16;
!         else
                nvs=nvs+2
                ivals(nvs-1:nvs)=cv
                jvals(nvs-1:nvs)=(/hpb, hpa/)
                vals(nvs-1:nvs)=(/-1., 1./)/dta
                vals1(nvs-1:nvs)=(/0.5, 0.5/)
       endif

    enddo

call init_sparse_dp(h2vddta,ivals,jvals,vals,nv,nh,nvs)
call init_sparse_dp(h2v,ivals,jvals,vals1,nv,nh,nvs)


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')


!%=============
!!!!!!!!!!!!!!       u2hddph, u2h
!%=============

write(*, '(" - d/dph and 1 for u-grid functions, evaluating on the h-grid ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do ch = 1,nh

       cph=hp(ch,1)
       cth=hp(ch,2)

      cphl=cph
      cphr=1+modulo((cph+1)-1,nph)
       upl=iu_ug(cphl,cth)
       upr=iu_ug(cphr,cth)


!       if ((hpa /= 0).and.(hpb /= 0)) then !% we are not at a Northern/Southern boundary

      !% 4th order FD
!      cphll=1+mod((cphl-1)-1,nph); cphrr=1+mod((cphr+1)-1,nph);
!      upll=iu_ug(cth,cphll); uprr=iu_ug(cth,cphrr);
!       if upll > 0 & uprr > 0 & fdflag == 4
!         nvs=nvs+4;
!         ivals(nvs-3:nvs)=ch;
!         jvals(nvs-3:nvs)=[upll upl upr uprr];
!         vals(nvs-3:nvs)=[1/24 -9/8 9/8 -1/24]/dph;
!         vals1(nvs-3:nvs)=[-1 9 9 -1]/16;
!       else
                nvs=nvs+2
                ivals(nvs-1:nvs)=ch
                jvals(nvs-1:nvs)=(/upl, upr/)
                vals(nvs-1:nvs)=(/-1., 1./)/dph
                vals1(nvs-1:nvs)=(/0.5, 0.5/)
!       endif

    enddo

call init_sparse_dp(u2hddph,ivals,jvals,vals,nh,nu,nvs)
call init_sparse_dp(u2h,ivals,jvals,vals1,nh,nu,nvs)



call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%=============
!!!!!!!!!!!!!!       v2hddta, v2h
!%=============

write(*, '(" - d/dta and 1 for v-grid functions, evaluating on the h-grid ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do ch = 1,nh

       cph=hp(ch,1)
       cth=hp(ch,2)

      ctha=cth+1
      cthb=cth
       vpa=iv_vg(cph,ctha)
       vpb=iv_vg(cph,cthb)

!       cthaa=ctha+1; cthbb=cthb-1;
!       vpaa=iv_vg(cthaa,cph); vpbb=iv_vg(cthbb,cph);
!       if vpaa > 0 & vpbb > 0 & fdflag == 4
!         nvs=nvs+4;
!         ivals(nvs-3:nvs)=ch;
!         jvals(nvs-3:nvs)=[vpbb vpb vpa vpaa];
!         vals(nvs-3:nvs)=[1/24 -9/8 9/8 -1/24]/dta;
!         vals1(nvs-3:nvs)=[-1 9 9 -1]/18;
!       else
                nvs=nvs+2
                ivals(nvs-1:nvs)=ch
                jvals(nvs-1:nvs)=(/vpb, vpa/)
                vals(nvs-1:nvs)=(/-1., 1./)/dta
                vals1(nvs-1:nvs)=(/0.5, 0.5/)
!       endif

    enddo

call init_sparse_dp(v2hddta,ivals,jvals,vals,nh,nv,nvs)
call init_sparse_dp(v2h,ivals,jvals,vals1,nh,nv,nvs)


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!**********************************************************************************************************
!**********************************************************************************************************
    deallocate(ivals, jvals, vals, vals1)
    allocate(ivals(4*nmax), jvals(4*nmax), vals(4*nmax), stat = status1)

    if (associated(fvals)) deallocate(fvals)
    allocate(fvals(4*nmax), stat = status2)

    if (associated(ispvals)) deallocate(ispvals) ! nmax for....
    allocate(ispvals(2*nmax), stat = istat)
    if (associated(jspvals)) deallocate(jspvals)
    allocate(jspvals(2*nmax), stat = statusj)
    if (associated(fspvals)) deallocate(fspvals)
    allocate(fspvals(2*nmax), stat = status1)

    ivals = 0
    jvals = 0
    vals = 0

    ispvals = 0
    jspvals = 0
    fspvals = 0
!**********************************************************************************************************
!**********************************************************************************************************
!%=============
!!!!!!!!!!!!!!       v2u, v2uf, v2ufsp
!%=============

write(*, '(" - v & fv,     			     evaluating on the u-grid ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do cu = 1,nu

       cph=up(cu,1)
       cth=up(cu,2)

       if (up(cu,3) == 0) then

       ctha=cth+1
       cthb=cth
       cphr=cph
       cphl=1+modulo((cph-1)-1,nph)

                nvs=nvs+4
                ivals(nvs-3:nvs)=cu
                jvals(nvs-3:nvs)=(/iv_vg(cphr,ctha), iv_vg(cphr,cthb), iv_vg(cphl,ctha), iv_vg(cphl,cthb)/)
                vals(nvs-3:nvs)=(/0.25, 0.25, 0.25, 0.25/)

                fvals(nvs-3:nvs)=0.25*calc_f_rot((/th_vg(ctha), th_vg(cthb), th_vg(ctha), th_vg(cthb)/), &
                                                 (/ph_ug(cph), ph_ug(cph), ph_ug(cph), ph_ug(cph)/), &
                                                 th0, omega)

!		    if flags.baro.num.f(2) == 2
!		      cf=1;
!		      while vp(jvals(nvs-4+cf)) == 1
!		        cf=cf+1;
!		      end
!		      ispvals(nvs/4)=cu;
!		      jspvals(nvs/4)=jvals(nvs-4+cf);
!		      fspvals(nvs/4)=4*fvals(nvs-4+cf);
!    elseif flags.baro.num.f(2) == 3
      ispvals(nvs/2-1:nvs/2)=cu;
      jspvals(nvs/2-1:nvs/2)=(/jvals(nvs-3), jvals(nvs)/)
      fspvals(nvs/2-1:nvs/2)=2*(/fvals(nvs-3), fvals(nvs)/)
!    end
       endif

    enddo

call init_sparse_dp(v2u,ivals,jvals,vals,nu,nv,nvs)
call init_sparse_dp(v2uf,ivals,jvals,fvals,nu,nv,nvs)

!if flags.baro.num.f(2) == 2
!  ispvals=ispvals(1:nvs/4);
!  jspvals=jspvals(1:nvs/4);
!  fspvals=fspvals(1:nvs/4);
!  v2ufsp=sparse(ispvals,jspvals,fspvals,nu,nv);
!  save(matfile,'-append','v2ufsp');
!  clear v2ufsp
!elseif flags.baro.num.f(2) == 3

call init_sparse_dp(v2ufsp,ispvals,jspvals,fspvals,nu,nv,nvs/2) ! nvs/2 elements!

!  save(matfile,'-append','v2ufsp');
!  clear v2ufsp
!end


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%=============
!!!!!!!!!!!!!!       u2v, u2vf, u2vfsp
!%=============

write(*, '(" - u & fu,     			     evaluating on the v-grid ")', advance = 'no')

    nvs = 0

call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    do cv = 1,nv

       cph=vp(cv,1)
       cth=vp(cv,2)

       if (vp(cv,3) == 0) then

       ctha=cth
       cthb=cth-1
       cphr=1+modulo((cph+1)-1,nph)
       cphl=cph


                nvs=nvs+4
                ivals(nvs-3:nvs)=cv
                jvals(nvs-3:nvs)=(/iu_ug(cphr,ctha), iu_ug(cphr,cthb), iu_ug(cphl,ctha), iu_ug(cphl,cthb)/)
                vals(nvs-3:nvs)=(/0.25, 0.25, 0.25, 0.25/)
                fvals(nvs-3:nvs)=0.25*calc_f_rot((/th_vg(cth), th_vg(cth), th_vg(cth), th_vg(cth)/), &
                                                 (/ph_ug(cphr), ph_ug(cphr), ph_ug(cphl), ph_ug(cphl)/), &
                                                 th0, omega)

!    if flags.baro.num.f(2) == 2
!      cf=1;
!      while up(jvals(nvs-4+cf)) == 1
!        cf=cf+1;
!      end
!      ispvals(nvs/4)=cv;
!      jspvals(nvs/4)=jvals(nvs-4+cf);
!      fspvals(nvs/4)=4*fvals(nvs-4+cf);
!              elseif flags.baro.num.f(2) == 3
                ispvals(nvs/2-1:nvs/2)=cv;
                jspvals(nvs/2-1:nvs/2)=(/jvals(nvs-3), jvals(nvs)/)
                fspvals(nvs/2-1:nvs/2)=2*(/fvals(nvs-3), fvals(nvs)/)
!              end
       endif

    enddo

call init_sparse_dp(u2v,ivals,jvals,vals,nv,nu,nvs)
call init_sparse_dp(u2vf,ivals,jvals,fvals,nv,nu,nvs)

!if flags.baro.num.f(2) == 2
!  ispvals=ispvals(1:nvs/4);
!  jspvals=jspvals(1:nvs/4);
!  fspvals=fspvals(1:nvs/4);
!  u2vfsp=sparse(ispvals,jspvals,fspvals,nv,nu);
!  save(matfile,'-append','u2vfsp');
!  clear u2vfsp
!elseif flags.baro.num.f(2) == 3

call init_sparse_dp(u2vfsp,ispvals,jspvals,fspvals,nv,nu,nvs/2) ! nvs/2 elements!

!  save(matfile,'-append','u2vfsp');
!  clear u2vfsp
!end


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%================================================
!% Deallocate what's possible
!%================================================
      deallocate(iu_ug, iv_vg, ih_hg, ph_ug, th_vg, ta_vg)

!%===================================================================================================================
!% save sparse matrices h2uddph, h2u, h2vddta, h2v, u2hddph, u2h, v2hddta, v2h, u2v, u2vf, u2vfsp, v2u, v2uf, v2ufsp
!%===================================================================================================================
write(*, '("Saving matrices: ")', advance = 'no')
      call save_cmat_dp(h2uddph, dir_mats // 'h2uddph.dat')
      call save_cmat_dp(h2u, dir_mats // 'h2u.dat')
      call save_cmat_dp(h2vddta, dir_mats // 'h2vddta.dat')
      call save_cmat_dp(h2v, dir_mats // 'h2v.dat')
      call save_cmat_dp(u2hddph, dir_mats // 'u2hddph.dat')
      call save_cmat_dp(u2h, dir_mats // 'u2h.dat')
      call save_cmat_dp(v2hddta, dir_mats // 'v2hddta.dat')
      call save_cmat_dp(v2h, dir_mats // 'v2h.dat')
      call save_cmat_dp(u2v, dir_mats // 'u2v.dat')
      call save_cmat_dp(v2u, dir_mats // 'v2u.dat')

      call save_cmat_dp(u2vfsp, dir_mats // 'u2vfsp.dat')
      call save_cmat_dp(v2ufsp, dir_mats // 'v2ufsp.dat')

!%================================================
!% Adjust the Coriolis matrices to conserve energy
!%================================================
      call write_modified_fmats_col(u2vf, v2uf, H_u, H_v, nu, nv)

      call save_cmat_dp(u2vf, dir_mats // 'u2vf.dat')
      call save_cmat_dp(v2uf, dir_mats // 'v2uf.dat')


end subroutine write_global_cmats_col

!==========================================================================================
!==========================================================================================

subroutine write_modified_fmats_col(u2vf, v2uf, H_u, H_v, nu, nv)
!%====================================================================
!% Different when solving in velocity or volume transport formulations
!%====================================================================
    implicit none

     integer, intent(in)   ::     nu, nv
     real(dp), dimension(:) :: H_u(nu), H_v(nv)

     type (triplet_dp)   :: u2vf, v2uf

!     integer ::      cu, cv, cz

!	Formulation in terms of velocity
!if flags.baro.num.uvform(domflag) == 1
!  v2uf=spdiags(1./sqrt(H_u),0,nu,nu)*v2uf*spdiags(sqrt(H_v),0,nv,nv);
!  u2vf=spdiags(1./sqrt(H_v),0,nv,nv)*u2vf*spdiags(sqrt(H_u),0,nu,nu);
!elseif flags.baro.num.uvform(domflag) == 2
!	Formulation in terms of volume transport

!  v2uf=spdiags(sqrt(H_u),0,nu,nu)*v2uf*spdiags(1./sqrt(H_v),0,nv,nv)
   call left_right_mult_diag_dp(v2uf, (H_u)**0.5, nu, 1/((H_v)**0.5), nv)

!  u2vf=spdiags(sqrt(H_v),0,nv,nv)*u2vf*spdiags(1./sqrt(H_u),0,nu,nu)
   call left_right_mult_diag_dp(u2vf, (H_v)**0.5, nv, 1/((H_u)**0.5), nu)

!end

end subroutine write_modified_fmats_col

!==========================================================================================
!==========================================================================================

function calc_f_rot(thc,phc,th0,omega)

implicit none

real(dp), intent(in)           :: thc(:), phc(:), th0, omega
real(dp), dimension(size(phc)) :: calc_f_rot

!% calculates the Coriolis parameters at a point on the
!% computational lat-lon grid. The computational grid is
!% shifted by ph0 in lon, and dropped down th0 from the
!% pole (in the direction of ph0).

calc_f_rot = 2*omega*(-sin(th0)*cos(thc)*cos(phc)+cos(th0)*sin(thc))

end function calc_f_rot

!**********************************************************************************************************
!**********************************************************************************************************

subroutine write_global_ccols_col(nph, nu, nv, nh, dir_grid, dir_cols, flag_itd_smooth, flag_baro_smooth)

    implicit none

!% Generate some useful vectors when using an Arakawa C grid:
!%
!%   H_u : H on the u-grid points
!%   H_v : H on the v-grid points
!%   H_h : H on the h-grid points
!%   H_sht_h : H on the h-grid points
!%   ta_u : tau on the u-grid points
!%   ta_v : tau on the v-grid points
!%   ta_h : tau on the h-grid points

     integer, intent(in) ::     nph!, nth
     integer, intent(in) ::     nu, nv, nh!, np
     real(dp), allocatable    ::     H_hg(:, :), H_sht_hg(:, :)
     real(dp), allocatable  :: ta_ug(:), ta_vg(:)!, ph_vg(:)
     integer, allocatable   :: up(:, :), vp(:, :), hp(:, :)

     real(dp), pointer :: ta_u(:), ta_v(:), ta_h(:)
     real(dp), pointer :: H_u(:), H_v(:), H_h(:), H_sht_h(:)

     integer ::      statusu, statusv, statush
     integer ::      cu, cv, ch, cth, cph, cphl, cphr, ctha, cthb
     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
     integer :: wall_t1, wall_t2, clock_rate, clock_max

     character(len=*) :: dir_grid, dir_cols
     logical 		:: flag_itd_smooth, flag_baro_smooth
!%=============================================
!% Load the necessary grid matrices and vectors
!%=============================================
     call load_alloc_vector(ta_ug, dir_grid // 'ta_ug.dat')
     call load_alloc_vector(ta_vg, dir_grid // 'ta_vg.dat')
     call load_alloc_matrix(up, dir_cols // 'up.dat')
     call load_alloc_matrix(vp, dir_cols // 'vp.dat')
     call load_alloc_matrix(hp, dir_cols // 'hp.dat')

	 call load_alloc_matrix(H_hg, dir_grid // 'H_hg.dat')

!write(*, '("Making vectors: ")')
write(*, '("Making vectors: ")', advance = 'no')
call tostring_set(rfmt='F12.1')

!write(*, '(" - tau,          evaluating on the grid-points...... ")', advance = 'no')
write(*, '("tau, ")', advance = 'no')
!%=============
!% initialize:
!%=============

!     tau
call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    if (associated(ta_u)) deallocate(ta_u)
    allocate(ta_u(nu), stat = statusu)
    if (associated(ta_v)) deallocate(ta_v)
    allocate(ta_v(nv), stat = statusv)
    if (associated(ta_h)) deallocate(ta_h)
    allocate(ta_h(nh), stat = statush)

ta_u=ta_ug(up(:,2)) ! 2: theta-index
!write(*, '(" 2 ")')
!print *, nv, shape(vp), maxval(vp(:,2)), nth, size(ta_vg)
ta_v=ta_vg(vp(:,2))
!write(*, '(" 3 ")')

ta_h=ta_ug(hp(:,2))
!write(*, '(" 4 ")', advance = 'no')

!call CPU_Time(T2)
!     call system_clock ( wall_t2, clock_rate, clock_max )
!     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!     H on u-grid
!write(*, '(" - H,            evaluating on the u-grid ")', advance = 'no')
write(*, '("H, evaluating on the u-, ")', advance = 'no')

!call CPU_Time(T1)
!call system_clock ( wall_t1, clock_rate, clock_max )

    if (associated(H_u)) deallocate(H_u)
    allocate(H_u(nu), stat = statusu)

do cu = 1, nu

      cph=up(cu,1)
      cth=up(cu,2)

      cphr=cph
      cphl=1+modulo((cph-1)-1,nph)

      if (up(cu,3) == 0) then
        H_u(cu) = (H_hg(cphl,cth)+H_hg(cphr,cth))/2
      else
        H_u(cu) = max( H_hg(cphl,cth), H_hg(cphr,cth) )
      endif

enddo

!call CPU_Time(T2)
!     call system_clock ( wall_t2, clock_rate, clock_max )
!     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!     H on v-grid
!write(*, '(" - H,            evaluating on the v-grid ")', advance = 'no')
write(*, '("v-, ")', advance = 'no')

!call CPU_Time(T1)
!call system_clock ( wall_t1, clock_rate, clock_max )

    if (associated(H_v)) deallocate(H_v)
    allocate(H_v(nv), stat = statusv)

do cv = 1, nv

      cph=vp(cv,1)
      cth=vp(cv,2)

      cthb=cth-1
      ctha=cth

      if (vp(cv,3) == 0) then
        H_v(cv) = (H_hg(cph,cthb)+H_hg(cph,ctha))/2
      else
        H_v(cv) = max( H_hg(cph,ctha), H_hg(cph,cthb) )
      endif

enddo

!call CPU_Time(T2)
!     call system_clock ( wall_t2, clock_rate, clock_max )
!     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!     H on h-grid
!write(*, '(" - H,            evaluating on the h-grid ")', advance = 'no')
write(*, '("h-grid ")', advance = 'no')

!call CPU_Time(T1)
!call system_clock ( wall_t1, clock_rate, clock_max )

    if (associated(H_h)) deallocate(H_h)
    allocate(H_h(nh), stat = statush)
    if (associated(H_sht_h)) deallocate(H_sht_h)
    allocate(H_sht_h(nh), stat = statush)

do ch = 1, nh
      cph=hp(ch,1)
      cth=hp(ch,2)

      H_h(ch) = H_hg(cph, cth)
enddo

!	SHT-smoothed bottom topography
if ( flag_itd_smooth ) then
    call load_alloc_matrix(H_sht_hg, dir_grid // 'H_sht_hg.dat')

	do ch = 1, nh
	      cph=hp(ch,1)
	      cth=hp(ch,2)

	      H_sht_h(ch) = H_sht_hg(cph, cth)
	enddo

    call save_vector(H_sht_h, dir_cols // 'H_sht_h.dat')
endif

!	SHT-smoothed is also used for the barotropic solver
!if ( flag_baro_smooth ) then
!    call load_alloc_matrix(H_sht_hg, dir_grid // 'H_orig_hg.dat')
!
!	do ch = 1, nh
!	      cph=hp(ch,1)
!	      cth=hp(ch,2)
!
!	      H_sht_h(ch) = H_sht_hg(cph, cth)
!	enddo
!
!    call save_vector(H_sht_h, dir_cols // 'H_orig_h.dat')
endif


call CPU_Time(T2)
     call system_clock ( wall_t2, clock_rate, clock_max )


     call disp ('(CPU: '//tostring(T2-T1)//'s, Wall:'//tostring(real(wall_t2-wall_t1)/real(clock_rate))//'s)')

!%=============================================
!% save columns H_u, H_v, H_h, ta_u, ta_v, ta_h
!%=============================================
     call save_vector(H_u, dir_cols // 'H_u.dat')
     call save_vector(H_v, dir_cols // 'H_v.dat')
     call save_vector(H_h, dir_cols // 'H_h.dat')
     call save_vector(ta_u, dir_cols // 'ta_u.dat')
     call save_vector(ta_v, dir_cols // 'ta_v.dat')
     call save_vector(ta_h, dir_cols // 'ta_h.dat')

end subroutine write_global_ccols_col

!==========================================================================================

end module generate_matrices



