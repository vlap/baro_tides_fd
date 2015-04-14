module itd

     use precisions, only: wp, cwp
     use my_trigs
     use my_sparse
     use my_sparse_aggregate
     use control
!     use generate_grid
     use save_load
     use interpolate
!     use spline

contains

!**********************************************************************************************************
!**********************************************************************************************************




!**********************************************************************************************************

subroutine calc_N_on_grid(N_u, N_v, P, nu, nv, nh, N_data_dir, dir_cols, dir_grid, dir_mats)

implicit none

	type(params) :: P
     real(wp), allocatable	:: ph_vg(:), th_ug(:)
     real(wp), allocatable	:: lon(:), lat(:), z(:), H_hg(:,:)
     integer, allocatable	:: hp(:,:)
     real(wp), allocatable	:: N_topo(:,:), N_hg(:,:)
     real(wp), allocatable	:: N_topo_3d(:,:,:), N_hg_3d(:,:,:)
     type (triplet)			:: h2v, h2u

     integer, intent(in)    :: nu, nv, nh
     integer            :: j, cz, istat, nlon, nlat, nz

     character(len = *)		:: dir_grid, dir_cols, dir_mats, N_data_dir
     character(len = 100)	:: N_data_file

     real(wp), allocatable, dimension(:) :: N_u, N_v, N_h


if (P%N_data .eq. '0') then !no data present
	N_u = P%Ns
	N_v = P%Ns
else
! load hp, ph_vg, ta_h
     call load_alloc_matrix(hp, dir_cols // 'hp.dat')
     call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
     call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat')

! load sparse matrices
    call load_alloc_sparse(h2u, dir_mats // 'h2u.dat')
	call load_alloc_sparse(h2v, dir_mats // 'h2v.dat')

! allocate N_u and N_v
    allocate(N_u(nu),N_v(nv), N_h(nh), stat = istat)

! load N data on the given grid
	N_data_file = N_data_dir // trim(adjustl(P%N_data)) // '.dat'
	! Check if WOA09 is being loaded (includes a 3D array with levels)
	if (scan(P%N_data,'9') > 0) then
		call load_alloc_topo_3d(N_data_file, nlon, nlat, nz, lon, lat, z, N_topo_3d, .false.)
		allocate(N_hg_3d(size(ph_vg), size(th_ug), nz), stat = istat)
		! interpolate horizontally onto our grid
		do cz = 1, nz
    		call bilinear_2d(nlon, nlat, d2r(lon), d2r(lat), N_topo_3d(:,:,cz), &
    						 size(ph_vg), size(th_ug), ph_vg, th_ug, N_topo, P%lib)
			N_hg_3d(:,:,cz) = N_topo
    	enddo
    	! now interpolate vertically
		call load_alloc_matrix(H_hg, dir_grid // 'H_hg.dat')
    	call interp_z_3d_to_2d(size(ph_vg), size(th_ug), nz, z, H_hg, N_hg_3d, N_hg)
    	deallocate(N_topo, N_topo_3d, N_hg_3d, H_hg)
	else
	! Assume WOA05 is being loaded (vertically averaged, 2D array)
		call load_alloc_topo(N_data_file, nlon, nlat, lon, lat, N_topo, .false.)
		! interpolate on our grid
    	call bilinear_2d(nlon, nlat, d2r(lon), d2r(lat), N_topo, size(ph_vg), size(th_ug), ph_vg, th_ug, N_hg, P%lib)
    	deallocate(N_topo)
	endif

! Because of spline interpolation there is a few negative values of N_hg
! Just turn them positive. So few that it doesn't matter
	N_hg = abs(N_hg)
!call save_matrix(N_hg, dir_grid // 'N_hg.dat')

 ! write N_u, N_v, N_h in columns
     N_u = 0
     N_v = 0
     N_h = 0

    do j = 1, nh
		N_h(j) = N_hg(hp(j,1),hp(j,2))
	end do
	call coo_vec_mul(h2u, N_h, P%lib, N_u)
	call coo_vec_mul(h2v, N_h, P%lib, N_v)

    call deallocate_sparse(h2u)
    call deallocate_sparse(h2v)
    deallocate(hp, ph_vg, th_ug, N_h)

end if

end subroutine calc_N_on_grid
!**********************************************************************************************************
subroutine interp_z_3d_to_2d(nlon, nlat, nz, z, H_hg, N_topo_3d, N_topo)
! given n, returns the load love numbers h'_n and k'_n.
	implicit none

	integer	:: nlon, nlat, nz
	real ( kind = 8 ), allocatable	:: ypp(:)
	real ( kind = 8 )				:: ypval, yppval, zero_dp = 0

	real(wp), intent(in)	:: z(nz), N_topo_3d(:,:,:), H_hg(:,:)
	real(wp), allocatable	:: N_topo(:,:)
	real ( kind = 8 )		:: N_tmp(nlon,nlat)
	integer :: istat, clon, clat


    allocate(ypp(nz), stat = istat)

do clon = 1,nlon
	do clat = 1,nlat
!	 set the spline to have 0 derivative at the right boundary
!						and 0 second derivative at the left endpoint
		call spline_cubic_set ( nz, real(z,kind=8), real(N_topo_3d(clon,clat,:),kind=8), 2, zero_dp, 1, zero_dp, ypp ) !0 - the spline should be a quadratic over the first and last intervals
		call spline_cubic_val ( nz, real(z,kind=8), real(N_topo_3d(clon,clat,:),kind=8), ypp, &
														real(H_hg(clon,clat),kind=8), N_tmp(clon,clat), ypval, yppval )
	enddo
enddo

! back to the variables in our precision
allocate(N_topo(nlon,nlat), stat = istat)
!overwrite extrapolated values for the deep ocean with deepest given N
where (H_hg > z(nz))
	N_topo = N_topo_3d(:,:,nz)
elsewhere
	N_topo = N_tmp
end where

end subroutine interp_z_3d_to_2d

!**********************************************************************************************************

subroutine calc_itd_prefactor(ufactor,vfactor, cpt, P, nu, nv, dir_cols, dir_grid)

implicit none

	type(params) :: P
     real(wp), allocatable	:: ta_vg(:), ph_vg(:), ta_ug(:), ph_ug(:)
     real(wp), allocatable	:: ta_v(:), ph_v(:), ta_u(:), ph_u(:), th_v(:), th_u(:)
     real(wp), allocatable	:: ph_u_nonrot(:), th_u_nonrot(:), ph_v_nonrot(:), th_v_nonrot(:)
	 real(wp)           :: th0, ph0
	 integer, allocatable   :: up(:,:), vp(:,:)

     integer, intent(in)    :: nu, nv!, nh
     integer            :: j, istat

	 character(len=2)   :: cpt
	 type (tide_params) :: pars
     character(len = *) :: dir_grid, dir_cols!, dir_mats, N_data_dir

     integer, allocatable, dimension(:) :: ufactor,vfactor


allocate(ufactor(nu), vfactor(nv), stat = istat)
ufactor = 0
vfactor = 0

if (P%trapped == 0) then !no do not cut off ITD conversion above crit lats
	ufactor = 1
	vfactor = 1
else
! load hp, ph_vg, ta_h
     call load_alloc_matrix(up, dir_cols // 'up.dat')
     call load_alloc_matrix(vp, dir_cols // 'vp.dat')

     call load_alloc_vector(ph_ug, dir_grid // 'ph_ug.dat')
     call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
     call load_alloc_vector(ta_ug, dir_grid // 'ta_ug.dat')
     call load_alloc_vector(ta_vg, dir_grid // 'ta_vg.dat')

    allocate(ta_v(nv), ph_v(nv), ta_u(nu), ph_u(nu), th_v(nv), th_u(nu), stat = istat)
! and fill them

	ph_u = (/ ( (ph_ug(up(j, 1))), j=1,nu ) /)
	ph_v = (/ ( (ph_vg(vp(j, 1))), j=1,nv ) /)
	ta_u = (/ ( (ta_ug(up(j, 2))), j=1,nu ) /)
	ta_v = (/ ( (ta_vg(vp(j, 2))), j=1,nv ) /)

	th_u=ta2th(ta_u, P%coor)
	th_v=ta2th(ta_v, P%coor)

	deallocate(ph_ug, ph_vg, ta_ug, ta_vg, ta_u, ta_v)
! poles
th0 = d2r(P%latP)
ph0 = d2r(P%lonP)
! transform to the non-rotated coords
	call calc_latlon_nonrot(ph_u, th_u, ph0, th0, ph_u_nonrot, th_u_nonrot)
	call calc_latlon_nonrot(ph_v, th_v, ph0, th0, ph_v_nonrot, th_v_nonrot)

	pars = get_pars(cpt)

	where ((abs(th_u_nonrot)) < asin(.5*pars%omega0/P%omega)) ufactor = 1
	where ((abs(th_v_nonrot)) < asin(.5*pars%omega0/P%omega))  vfactor = 1

	call save_vector(ufactor, dir_cols // 'ufactor_'//cpt//'.dat')
	call save_vector(vfactor, dir_cols // 'vfactor_'//cpt//'.dat')

!	print *, sum(ufactor), size(ufactor)
!	print *, sum(vfactor), size(vfactor)

	deallocate(ph_u, th_u,ph_v, th_v, ph_u_nonrot, th_u_nonrot, ph_v_nonrot, th_v_nonrot)

end if

end subroutine calc_itd_prefactor

!**********************************************************************************************************

end module itd
