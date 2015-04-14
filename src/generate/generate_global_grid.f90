module generate_grid

     use precisions, only: dp, eps
     use my_trigs
     use err_manager
     use dispmodule
	 use control
     use save_load
	 use read_etopo
	 use interpolate
     use spherepack_iface

     use testshrotate

     implicit none

!********************************************************************

type grid_dims

     integer ::     coor ! 1 for Mercator, 2 for lat-lon
     integer ::     nph, nta
     integer ::     nu=0, nv=0, nh=0, np=0
     real(dp)::		lonP, latP

end type

     contains
!***********************************************************************************************************
subroutine generate_global_grid(etopo_file, topo_file, topo_dir_in, topo_dir_out, dir_grid, dir_cols, P, GD)

     character(len=*) :: etopo_file, topo_file, topo_dir_in, topo_dir_out, dir_grid, dir_cols
     type(params)	  :: P
     type(grid_dims)  :: GD

     integer :: numLons, numLats
     real(dp), allocatable :: xValues(:), yValues(:)
     integer, allocatable  :: zValues(:, :)
     real(dp), allocatable :: ph(:), th(:)!, topo(:)

     real(dp), allocatable :: ph_vg(:), th_ug(:) !ta_ug(:), ta_vg(:),
     real(dp), allocatable :: H_hg(:, :), H_sht_hg(:,:)

     real(dp), allocatable	:: smooth(:)
     integer				:: j,k
     real(dp), allocatable	:: H_sht(:,:), topo_sht(:,:)
     real(dp), allocatable	:: lon_sht(:), lat_sht(:)
     integer			    :: nlon_sht, nlat_sht, istat

     character(len=10)	:: str_tmp
     character(len=50)	:: filename
     character(len=100) :: dirname
     logical			:: load_smooth_file

     real    ::      T1, T2! for measuring CPU (NOT REAL TIME!)
     integer :: 	 wall_t1, wall_t2, clock_rate, clock_max

!%================================================
!% Loading/preparing topography file
!%================================================
!%	Test rotate SH coeffcients
!	call test_sh_rotate(trim(topo_dir_out), trim(topo_file), P)
!	stop
!%================================================

     if (P%load_etopo==1) then
		numLons = P%etopo_res;
		numLats = P%etopo_res/2;
		call prepare_topo(etopo_file, numLons, numLats, xValues, yValues, zValues, &
                            real(P%lonP, kind=dp), real(P%latP, kind=dp))
		! add to the filename the resolution in minutes
		write (str_tmp, '(g10.2)') 360*60/real(numLons-1)
		filename = 'topo_rot_'//trim(adjustl(str_tmp))//'min_pole_'//tostring(P%latP)//'_'//tostring(P%lonP)//'.dat'
		dirname = topo_dir_out
		!write to the file
		call save_topo(topo_dir_out // trim(filename), numLons, numLats, xValues, yValues, zValues)

     else
		  filename = topo_file
		  dirname = topo_dir_in
          call load_alloc_topo(topo_dir_in//topo_file, numLons, numLats, xValues, yValues, zValues, (P%messages>=1))
     endif

! convert xValues and yValues to radians
    allocate(ph(size(xValues)), th(size(yValues)), stat = istat)
    ph=d2r(xValues)
    th=d2r(yValues)
	deallocate(xValues, yValues)

!%================================================
!% set up the global C-grid and operator matrices
!%================================================
	call gen_global_grid(zValues, ph, th, numLons-1, H_hg, dir_grid, P)
      ! (ph_ug, ph_vg, ta_ug, ta_vg) are saved into /global/cols

  GD%coor = P%coor
  GD%lonP = P%lonP
  GD%latP = P%latP
  GD%nph = P%nph
  GD%nta = size(H_hg,2)

      if ((minval(H_hg(:,1)) < 0.).or.(minval(H_hg(:,size(H_hg,2))) < 0.)) then
	      write(*, '("****** Ocean is in the first/last tau index of H_hg. Removed ******")')
    	    where (H_hg(:,1)<0) H_hg(:,1)=real(0,dp)! new depth is 0
	        where (H_hg(:,size(H_hg,2))<0) H_hg(:,size(H_hg,2))=real(0,dp)! new depth is 0
	  end if

      ! set the minimum permitted water depth on the global grid:
      call adjust_min_depth(H_hg, GD%nph, GD%nta, P, dir_grid)
!call save_matrix(H_hg, dir_grid // 'H_hg_raw.dat')

      ! bays
      if (P%bays == 1) then
      	call remove_bays(H_hg, GD%nph, GD%nta, real(0,dp))!, dir_grid)  ! new depth is 0
      end if
!call save_matrix(H_hg, dir_grid // 'H_hg_raw1.dat')
      ! seas
      if (P%inlandseas == 1) then
      	call remove_seas(H_hg, GD%nph, GD%nta, real(0,dp), GD%nph/10)  ! new depth is 0
      end if
!call save_matrix(H_hg, dir_grid // 'H_hg_raw2.dat')
!     Make H_hg POSITIVE, until here below the sea level H_hg was NEGATIVE.
      H_hg = -H_hg
!     save H on the C-grid
      call save_matrix(H_hg, dir_grid // 'H_hg.dat')

!%============================================
! Smooth bottom topography for the ITD scheme
!%============================================
if ( (P%itd_scheme>0).and.(P%smooth_type>0).and.(P%sht_smooth_H >= 0) ) then
	call CPU_Time(T1)
	call system_clock ( wall_t1, clock_rate, clock_max )

	write(*, '("Smoothing topography for the linear ITD scheme: ")', advance = 'no')
!		Attempt to load the smoothed topography

		inquire( file=trim(dirname)//'smoothed_t'//tostring(P%smooth_type)//'_trunc'//tostring(P%sht_smooth_H)// '_'// &
   	               				 trim(filename), exist=load_smooth_file )
		if (load_smooth_file) then
   			call load_alloc_topo(trim(dirname)//'smoothed_t'//tostring(P%smooth_type)//'_trunc'//tostring(P%sht_smooth_H)// '_'// &
   	               				 trim(filename), nlon_sht, nlat_sht, lon_sht, lat_sht, H_sht, (P%messages>=1))
		else
	     	call sht_smooth_H(ph, th,real(zValues,dp), P, lon_sht, lat_sht, H_sht, trim(dirname), trim(filename))
	    endif
		! DEALLOCATE TOPO
		deallocate(ph, th, zValues)

     	! shortcuts
		nlon_sht = size(lon_sht)
		nlat_sht = size(lat_sht)
		call load_alloc_vector(ph_vg, dir_grid // 'ph_vg.dat')
		call load_alloc_vector(th_ug, dir_grid // 'th_ug.dat')
	! Now restrict the smoothed topography (for ITD) onto the grid topo H_hg (interpolate)
	    call bilinear_2d(nlon_sht, nlat_sht, lon_sht, lat_sht, H_sht, &
	                                         GD%nph, GD%nta, ph_vg, th_ug, H_sht_hg, P%lib)
		deallocate(lon_sht, lat_sht, H_sht, ph_vg, th_ug)

		where (H_hg <= 0) H_sht_hg = 0
        call save_matrix(H_sht_hg, dir_grid // 'H_sht_hg.dat')
        deallocate(H_sht_hg)

	call CPU_Time(T2)
	call system_clock ( wall_t2, clock_rate, clock_max )
	call disp ('done (CPU: ' // trim(ADJUSTL(conv_secs(T2-T1))) //', Wall: '&
						//trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) )))//')')
endif

!     print *,  P%nph, GD%nta, P%nph*GD%nta, count(H_hg>0), count(H_hg>0)/(P%nph*GD%nta)
     call write_global_cgrid(H_hg, GD%nph, GD%nta, GD%nu, GD%nv, GD%nh, GD%np, dir_grid, dir_cols)

end subroutine

!***********************************************************************************************************

subroutine write_global_cgrid(H, nph, nth, nu, nv, nh, np, dir_grid, dir_cols)

    implicit none
!% for a height field H(nph, nth), defined on the h-grid,
!% generates appropriate u, v and h grids in an Arakawa C manner.
!% Returns :
!%
!%  up(1:nu,1:4) : phi-index, theta-index, boundary-flag, opt. position
!%  vp(1:nv,1:4) : phi-index, theta-index, boundary-flag, opt. position
!%  hp(1:nh,1:4) : phi-index, theta-index,             0, opt. position
!%
!%  upmat(nph, nth): zero if no grid-point, otherwise the u-th point
!%  vpmat(nph, nth+1) : zero if no grid-point, otherwise the v-th point
!%  hpmat(nph, nth) : zero if no grid-point, otherwise the h-th point

     integer, intent(in)     ::     nph, nth
     real(dp), allocatable:: H(:, :)
     integer, allocatable :: up(:, :), vp(:, :), hp(:, :)

     integer, allocatable :: upmat(:, :), vpmat(:, :), hpmat(:, :)
     integer ::     nu, nv, nh, np

     integer ::      statusu, statusv, statush
     integer ::      cth, cph, cphl, cphr

     character(len=*) :: dir_grid, dir_cols

write(*, '("Allocating grid points: ")', advance = 'no')

    if (allocated(up)) deallocate(up)
    allocate(up(nph*nth,4), stat = statusu)
    if (allocated(vp)) deallocate(vp)
    allocate(vp(nph*nth,4), stat = statusv)
    if (allocated(hp)) deallocate(hp)
    allocate(hp(nph*nth,4), stat = statush)

up = 0
vp = 0
hp = 0

    if (allocated(upmat)) deallocate(upmat)
    allocate(upmat(nph,nth), stat = statusu)
    if (allocated(vpmat)) deallocate(vpmat)
    allocate(vpmat(nph,nth+1), stat = statusv)
    if (allocated(hpmat)) deallocate(hpmat)
    allocate(hpmat(nph,nth), stat = statush)

upmat = 0
vpmat = 0
hpmat = 0

do cph = 1,nph
!print *, ""
!write(*, '(i3)') cph
     do cth = 1,nth
!write(*, '(i3,  " ")', advance = 'no') cth
           cphl=1+modulo((cph-1)-1,nph)
           cphr=1+modulo((cph+1)-1,nph)

         if (H(cph, cth) > 0) then
         if ( (H(cphr, cth) > 0) .or. (H(cphl, cth) > 0) .or. &
                (H(cph, cth-1) > 0) .or. (H(cph, cth+1) > 0) ) then

             !% allot an h-grid point:
             nh=nh+1
             np=np+1

             hp(nh,1)=cph
             hp(nh,2)=cth
             hp(nh,4)=np
             hpmat(cph, cth)=nh

             !% ... and some u-grid points:
             nu=nu+1
             np=np+1
             up(nu,1)=cph
             up(nu,2)=cth
             up(nu,4)=np
             upmat(cph, cth)=nu

             if (H(cphl, cth) <= 0) then
               up(nu,3)=1
             endif

             if (H(cphr, cth) <= 0) then
               nu=nu+1
               np=np+1
               up(nu,1)=cphr
               up(nu,2)=cth
               up(nu,3)=1
               up(nu,4)=np
               upmat(cphr, cth)=nu
             endif

             !% .. and some v-grid points:
             nv=nv+1
             np=np+1
             vp(nv,1)=cph
             vp(nv,2)=cth
             vp(nv,4)=np
             vpmat(cph,cth)=nv
             if (H(cph, cth-1) <= 0) then
               vp(nv,3)=1
             endif
             if (H(cph, cth+1) <= 0) then
               nv=nv+1
               np=np+1
               vp(nv,1)=cph
               vp(nv,2)=cth+1
               vp(nv,3)=1
               vp(nv,4)=np
               vpmat(cph,cth+1)=nv
             endif

         endif
         endif

     enddo
enddo

write(*, '( a, " for u, ", a, " for v, ", a, "  for h.")') tostring(nu),tostring(nv),tostring(nh)

!%=============================================
!% save the forward and backward C grid indexes
!%=============================================
! save iu_ug, iv_vg, ih_hg
     call save_matrix(upmat, dir_grid // 'iu_ug.dat')
     call save_matrix(vpmat, dir_grid // 'iv_vg.dat')
     call save_matrix(hpmat, dir_grid // 'ih_hg.dat')
! save up, vp, hp
     call save_matrix(up(1:nu,1:4), dir_cols // 'up.dat')
     call save_matrix(vp(1:nv,1:4), dir_cols // 'vp.dat')
     call save_matrix(hp(1:nh,1:4), dir_cols // 'hp.dat')

end subroutine write_global_cgrid

!***********************************************************************************************************

subroutine gen_global_grid(topo, ph, th, nlon, H_hg, dir_grid, P)!, lonP, latP)

    implicit none
!% receives the topography data
!%
!%  topo(nlon+1, nlon/2+1), lat(nlon/2+1), lon(nlon+1), latp, lonp
!%
!% where topo is a lat x lon grid of height (m) above sea level,
!% lat and lon are the corresponding arbitrarily spaced vectors,
!% and latp and lonp give the position of the rotated pole
!% (all in degrees).
!%
!% Generates whatever grids are necessary, and updates flags as appropriate.
!% Uses the input resolution flags.nph.

     integer, intent(in)::     nlon
     integer			::     nph
     real(dp), allocatable :: ph(:), th(:)
     integer, allocatable   :: topo(:, :)
!     real(dp)     ::     latP, lonP
     type(params) :: P
     real(dp)  ::      dph


     real(dp), allocatable :: ph_ug(:), ta_ug(:), th_ug(:), ph_vg(:), ta_vg(:), th_vg(:)
     real(dp), allocatable :: H_hg(:, :)

     integer ::      statusx, statusy!, statusz
     integer ::      j, js, jn
!     real(dp)  ::      dx, dy, dxR, dyR

     integer ::     nn, ns, nta
     real(dp)  ::      ta_s, ta_n
     integer ::     numLons, numLats

     character(len=*) :: dir_grid

!	real(dp), allocatable :: ph1(:), th1(:), h1(:,:)


!'nph','nta','ta_ug','ta_vg','ph_ug','ph_vg','th_ug','th_vg'
!%======================
!% make some short-cuts
!%======================
nph = P%nph
numLons = nlon+1
numLats = nlon/2+1

    dph=2*pi/nph;    ! the distance between gridpoints

!%==================================================================
!% calculate js (the southernmost index where ocean first appears):
!%==================================================================

js=1;
do while (minval(topo(:, js)) >= 0)
  js=js+1;
end do

!% calculate the corresponding (negative) value of tau:
ta_s=th2ta(th(js), P%coor);

!print *, "th(js), ta_s, th2ta(ta_s)", th(js), ta_s, ta2th(ta_s, P%coor)

!% this will lie in the ns-th coarse tau box south of the equator, where
ns=ceiling(-ta_s/dph);
!%==================================================================
!% calculate jn (the northernmost index where ocean first appears):
!%==================================================================

jn=numLats;
do while (minval(topo(:, jn)) >= 0)
  jn=jn-1;
end do

!% calculate the corresponding value of tau:
ta_n=th2ta(th(jn), P%coor);
!print *, "th, ta_n, th2ta(ta_n)", th(jn), ta_n, th2ta(ta_n, P%coor)

!% this will lie in the nn-th coarse tau box north of the equator, where
nn=ceiling(ta_n/dph);
!print *, "ns", ns
!print *, "nn", nn
!print *, r2d(ta_n), r2d(ta_s)
!%===================================
!% Generate global grid coordinates:
!%===================================
    if (allocated(ph_ug)) deallocate(ph_ug)
    allocate(ph_ug(nph), stat = statusx)
    if (allocated(ph_vg)) deallocate(ph_vg)
    allocate(ph_vg(nph), stat = statusx)

    nta=nn+ns+2

    if (allocated(ta_ug)) deallocate(ta_ug)
    allocate(ta_ug(nta), stat = statusy)
    if (allocated(ta_vg)) deallocate(ta_vg)
    allocate(ta_vg(nta+2), stat = statusy)
    if (allocated(th_ug)) deallocate(th_ug)
    allocate(th_ug(nta), stat = statusy)
    if (allocated(th_vg)) deallocate(th_vg)
    allocate(th_vg(nta+2), stat = statusy)

ph_ug=(/ ( (2*pi*j/nph), j=0,nph-1 ) /)
ph_vg=ph_ug+pi/nph

!ta_vg=(/ ( (j*(2*pi/nph)), j=-ns-2,nn+1 ) /)  ! -ns-1 doesn't always really work...
!th_vg=(/ ( (ta2th(ta_vg(j+ns+2))), j=-ns-2,nn+1 ) /) ! +ns+2 to adjust the lower index
ta_vg=(/ ( ((j-ns-3)*(2*pi/nph)), j=1,nta+2 ) /)  ! -ns-1 doesn't always really work...
!th_vg=(/ ( (ta2th(ta_vg(j), P%coor)), j=1,nn+ns+4 ) /) ! +ns+2 to adjust the lower index
th_vg = ta2th(ta_vg, P%coor)

ta_ug=(/ ( (ta_vg(j)+pi/nph), j=1,nta ) /)
!th_ug=(/ ( (ta2th(ta_ug(j), P%coor)), j=1,nta ) /)
th_ug = ta2th(ta_ug, P%coor)

!%=============================================
!% Interpolate bathymetry to global grid H_hg
!%=============================================
!     call interp_H_grid(numLons, numLats, ph, th, topo, nph, nta, ph_vg, th_ug, H_hg)
	write(*, '("Interpolating topo to ", f5.2, " min resolution global grid: ")', advance = 'no') 360*60/real(nph)
    call bilinear_2d(numLons, numLats, ph, th, real(topo,dp), nph, nta, ph_vg, th_ug, H_hg, P%lib)
	print *, "done."

!%==========================
!% save the coarse C grid
!%==========================
     call save_vector(ph_ug, dir_grid // 'ph_ug.dat')
     call save_vector(ph_vg, dir_grid // 'ph_vg.dat')
     call save_vector(ta_ug, dir_grid // 'ta_ug.dat')
     call save_vector(ta_vg, dir_grid // 'ta_vg.dat')
     call save_vector(th_ug, dir_grid // 'th_ug.dat')
     call save_vector(th_vg, dir_grid // 'th_vg.dat')

!! test bilinear_grid2grid
!	do i = 1, nph
!     H_hg(i,:) = 2*i
!    end do
!	do i = 1, nta
!     H_hg(:,i) = H_hg(:,i)-4*i
!    end do

!	SAVE GENERATED GRID
!     call save_topo(dir_grid // 'topo_orig.dat', nph, nta, ph_vg, th_ug, H_hg)


!    allocate(ph1(nph*10), stat = statusx)
!    allocate(th1(nta*10), stat = statusx)
!    ph1=(/ ( (2*pi*i/(nph*10))+pi/(nph*10), i=0,nph*10-1 ) /)
!    th1=(/ ( ( (th_ug(nta)-th_ug(1))/( nta*10-1 )*i + th_ug(1) ), i=0,nta*10-1 ) /)
!
!      call bilinear_grid2grid(nph, nta, ph_vg, th_ug, real(H_hg,dp), nph*10, nta*10, ph1, th1, h1)!, dir_grid)
!
!      !write to the file
!      call save_topo(dir_grid // 'topo_interp.dat', nph*10, nta*10, ph1, th1, h1)

end subroutine gen_global_grid
!***********************************************************************************************************

subroutine remove_seas(H_hg, nph, nta, Hstar, cmin)

    implicit none
!% receives the topography data H_hg(nph, nta)
!%
!% cmin gives the mimum number of connected points for "sea" to be significant
!% the removed seas are assigned H value Hstar

     integer ::     nta, nph
     real(dp) :: H_hg(nph, nta)
     logical :: in_class(nph, nta)
     integer :: sea_class(nph, nta) ! -1 for land, 0 for unassigned, 1,2,3... for diff seas
     integer :: n_class, cur_class, n_removed
     integer, intent(in) ::      cmin
     real(dp), intent(in) ::      Hstar
     integer, allocatable :: size_class(:)

     integer ::      status
     integer ::      i
     integer :: cph, cta, cph_tmp, cta_tmp

n_class = 0

!ind_b = (H_hg .le. 0 .and. ( (cshift(H_hg, shift=1, dim=1) .ge. 0) .or. &
!                             (cshift(H_hg, shift=-1, dim=1) .ge. 0) .or. &
!                             (cshift(H_hg, shift=1, dim=2) .ge. 0) .or. &
!                             (cshift(H_hg, shift=-1, dim=2) .ge. 0) )
!
!n_b = count(ind_b)

write(*, '("Removing inland seas/lakes: ")', advance = 'no')

sea_class = 0
cur_class = 0

do cph = 1, nph
     do cta = 1, nta

         if (H_hg(cph,cta) .ge. 0) then
              sea_class(cph,cta) = -1

         elseif (sea_class(cph,cta) .eq. 0) then
!         print *, cph, cta, H_hg(cph,cta), sea_class(cph,cta)
              cur_class = cur_class + 1
              sea_class(cph,cta) = cur_class
              cph_tmp = cph
              cta_tmp = cta
              call mark_connected(H_hg, nph, nta, cph_tmp, cta_tmp, sea_class, cur_class, n_class)
         endif
    enddo
enddo

write(*, '("found ", a, " connected components;")',advance='no') tostring(n_class)

if (allocated(size_class)) deallocate(size_class)
    allocate(size_class(n_class+2), stat = status)

n_removed = -1

do i = -1, n_class

    in_class = sea_class .eq. i
    size_class(i+2) = count(in_class)

    if ((i == -1).or.(size_class(i+2) < cmin)) then
         where (in_class) H_hg = Hstar ! assign height on the land to Hstar
         n_removed = n_removed + 1
    endif

enddo
write(*, '(" removed ", a, " of size < ", a, " cells.")') tostring(n_removed), tostring(cmin)
!write(*,'(200i8)') ( size_class(i+2), i=-1,n_class )

!        print *, " done."

end subroutine remove_seas

!***********************************************************************************************************

subroutine mark_connected(H_hg, nph, nta, cph, cta, sea_class, cur_class, n_class)
    implicit none

     integer ::     nta, nph
     real(dp) :: H_hg(nph, nta)
     integer :: sea_class(nph, nta) ! -1 for land, 0 for unassigned, 1,2,3... for diff seas
     integer :: cur_class, n_class
     integer :: cph, cta, cphl, cphr
     integer :: buf(nph*nta*3/4, 2) ! FIFO buffer for neighbours' coordinates
     integer :: bufp, bufe

n_class = n_class + 1

!write( *, '(4i6)', advance ='no' )  cph, cta, H_hg(cph,cta), n_class,cur_class

bufp = 1
bufe = 1

call push(buf, cph, cta, bufe)

do while (pop(buf, cph, cta, bufe, bufp))

     if (cta < nta) then
       if ((H_hg(cph, cta+1) < 0) .and. (sea_class(cph, cta+1) .eq. 0)) then
         !write( *, '(" 1")', advance = 'no')
         sea_class(cph, cta+1) = cur_class
         !print *, " 1", cph, cta+1
         call push(buf, cph, cta+1, bufe)
       endif
     endif

     if (cta > 1) then
       if ((H_hg(cph,cta-1) < 0) .and. (sea_class(cph,cta-1) .eq. 0)) then
!         write( *, '(" 4")', advance = 'no')
         sea_class(cph,cta-1) = cur_class
!         print *, " 4", cph, cta-1
         call push(buf, cph, cta-1, bufe)
       endif
     endif

     cphl=1 + modulo((cph-1)-1,nph)
     cphr=1 + modulo((cph+1)-1,nph)
!     print *, cphl, cphr

     if ((H_hg(cphr,cta) < 0) .and. (sea_class(cphr,cta) .eq. 0)) then
 !    print *, " 2"
         sea_class(cphr,cta) = cur_class
         call push(buf, cphr, cta, bufe)
     endif

     if ((H_hg(cphl,cta) < 0) .and. (sea_class(cphl,cta) .eq. 0)) then
     !print *, " 3"
         sea_class(cphl,cta) = cur_class
         call push(buf, cphl, cta, bufe)
     endif

enddo

!print *, bufe, bufp

end subroutine mark_connected

!***********************************************************************************************************

subroutine push(buf, cph, cta, bufe)
    implicit none

     integer :: cph, cta
     integer :: buf(:,:) !(nph*nta*2/3, 2)
     integer :: bufe

buf(bufe, 1) = cph
buf(bufe, 2) = cta
bufe = bufe + 1

end subroutine push

logical function pop(buf, cph, cta, bufe, bufp)
    implicit none

     integer :: cph, cta
     integer :: buf(:,:) !(nph*nta*2/3, 2)
     integer :: bufe, bufp
if (bufe == bufp) then
     pop = .false.
     return
else
     cph = buf(bufp, 1)
     cta = buf(bufp, 2)
     bufp = bufp + 1
     pop = .true.
     return
endif

end function pop

!***********************************************************************************************************

subroutine remove_bays(H_hg, nph, nta, Hstar)!,dir_grid)

    implicit none
!% receives the topography data H_hg(nph, nta) and removes bays
!% and inlets which are just one grid cell wide.
!% the removed seas are assigned H value Hstar

     integer, intent(in)  ::     nta, nph
     real(dp)			  :: H_hg(nph, nta)
     real(dp), intent(in) :: Hstar
     ! true if boundary is present on the  left, right, above, below
     logical :: l_b(nph, nta), r_b(nph, nta), u_b(nph, nta), b_b(nph, nta)
!     integer :: tmp(nph, nta)
     logical :: ind_3b(nph, nta) ! true when 3 boundaries are around a water cell
     integer :: n_removed = 1

!	character(len=*) :: dir_grid

write(*, '("Removing bays and inlets: ")', advance='no')

do while (n_removed > 0)

! find position of cells which has 3 sides on the boundary
r_b = (H_hg < 0.) .and. (cshift(H_hg, shift=1, dim=1) >= 0.)
l_b = (H_hg < 0.) .and. (cshift(H_hg, shift=-1, dim=1) >= 0.)
b_b = (H_hg < 0.) .and. (cshift(H_hg, shift=1, dim=2) >= 0.)
u_b = (H_hg < 0.) .and. (cshift(H_hg, shift=-1, dim=2) >= 0.)

ind_3b = (l_b .and. r_b .and. u_b) .or. (l_b .and. r_b .and. b_b) .or. &
         (l_b .and. b_b .and. u_b) .or. (r_b .and. b_b .and. u_b)

    where (ind_3b) H_hg = Hstar ! assign height on the land to Hstar

    n_removed = count(ind_3b)
    write(*, '(a)', advance = 'no') tostring(n_removed)
    if (n_removed == 0) then
    	write(*, '(".")')
    else
    	write(*, '(" + ")', advance = 'no')
    endif

end do

!call save_matrix(cshift(H_hg, shift=1, dim=1), dir_grid // 'l_shift.dat')
!
!tmp = 0
!where (l_b) tmp=1
!call save_matrix(tmp, dir_grid // 'l_b.dat')
!tmp = 0
!where (r_b) tmp=1
!call save_matrix(tmp, dir_grid // 'r_b.dat')
!tmp = 0
!where (u_b) tmp=1
!call save_matrix(tmp, dir_grid // 'u_b.dat')
!tmp = 0
!where (b_b) tmp=1
!call save_matrix(tmp, dir_grid // 'b_b.dat')

!        print *, "Done."

end subroutine remove_bays

!***********************************************************************************************************

subroutine adjust_min_depth(H_hg, nph, nta, P, dir_grid)

    implicit none
!% receives the topography data H_hg(nph, nta) and
!% adjusts the H_hg by setting the min depth
! scheme 1: a simple minimum depth
! scheme 2: min number of Grid Points Per Wavelength

     integer               ::     nta, nph
     real(dp)              :: H_hg(nph, nta)
     real(dp)              :: ta_ug(nta)

     type(params), intent(in)	:: P
     type (tide_params)	   		:: t_params(P%ncpts)

     integer   :: scheme, gppw, coor
     real(dp)  ::      hmin
     real(dp)  ::      g                 ! surface gravity, m/s^2
     real(dp)  ::      re                ! average earth's radius, m
!     real(dp)  ::      omega             ! angular velocity, rad/s

     real(dp) :: dx(nph, nta), kmax(nph, nta), hmin_m(nph, nta)
     integer  :: i
     real(dp) :: omegamax

     character(len=*) :: dir_grid

coor = P%coor
scheme = P%hmin_scheme
hmin =  P%hmin
gppw = P%gppw
g = P%g
re = P%re
!omega = P%omega

if (scheme == 1) then

     where ((H_hg < 0.).and.(H_hg > -hmin)) H_hg = -hmin

elseif (scheme == 2) then

!    omegamax = 2 * omega ! assume no higher constituents (m2 has the highest frequency)
	t_params = get_params(P%ncpts, P%cpts)
	omegamax = maxval(t_params%omega0)

     if (coor == 1) then ! Mercator

          call load_vector(ta_ug, dir_grid // 'ta_ug.dat')

          do i = 1, nph
              dx(i, :) = (2*pi*re/nph)*(1/cosh(ta_ug))
          enddo

     else if  (coor == 2) then ! lat-lon

          dx = 2*pi*re/nph

     endif

  kmax = 2*pi/(gppw*dx)   ! from min_wavelength = 2*pi/kmax = *MUST BE* = dx*gppw;

  hmin_m = g * (omegamax/kmax)**2

  where (hmin_m < hmin) hmin_m = hmin ! adjust for global minimum

  where ((H_hg > -hmin_m).and.(H_hg < 0.)) H_hg = -hmin_m ! now apply to H

endif

end subroutine adjust_min_depth
!***********************************************************************************************************

!********************************************************************
! Write the parameters of the grid into a specified file
!********************************************************************
subroutine write_GD(file_name, GD)

  implicit none

  character(len=*) :: file_name
  type(grid_dims):: GD

  integer :: fh = 15

  open(fh, file=file_name, action='write', status='replace')

  write(fh, '(a)') '! Basic parameters of the generated grid'
  write(fh, '(a)') 'coor = ' // tostring(GD%coor)
  write(fh, '(a)') 'nph = ' // tostring(GD%nph)
  write(fh, '(a)') 'nta = ' // tostring(GD%nta)
  write(fh, '(a)') 'nu = ' // tostring(GD%nu)
  write(fh, '(a)') 'nv = ' // tostring(GD%nv)
  write(fh, '(a)') 'nh = ' // tostring(GD%nh)
  write(fh, '(a)') 'np = ' // tostring(GD%np)
  write(fh, '(a)') 'lonp = ' // tostring(GD%lonP)
  write(fh, '(a)') 'latp = ' // tostring(GD%latP)

  close(fh)

end subroutine

!********************************************************************
! Read the parameters of the grid into a specified file
!********************************************************************
subroutine read_GD(file_name, GD)

  implicit none

  character(len=*) :: file_name
  type(grid_dims):: GD

  character(len=100) :: buffer, label
  integer :: pos, pos_
  integer :: fh = 15
  integer :: ios = 0
  integer :: line = 0


  open(fh, file=file_name, action='read')

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

	        select case (label)
!********************************************************************
	        case ('coor')
	           read(buffer, *, iostat=ios) GD%coor
	           !print*, 'Read coor: ', GD%coor
	        case ('nph')
	           read(buffer, *, iostat=ios) GD%nph
	           !print*, 'Read nph: ', GD%nph
	        case ('nta')
	           read(buffer, *, iostat=ios) GD%nta
	           !print*, 'Read nta: ', GD%nta
	        case ('nu')
	           read(buffer, *, iostat=ios) GD%nu
	           !print*, 'Read nu: ', GD%nu
	        case ('nv')
	           read(buffer, *, iostat=ios) GD%nv
	           !print*, 'Read nv: ', GD%nv
	        case ('nh')
	           read(buffer, *, iostat=ios) GD%nh
	           !print*, 'Read nh: ', GD%nh
	        case ('np')
	           read(buffer, *, iostat=ios) GD%np
	           !print*, 'Read np: ', GD%np
	        case ('lonp')
	           read(buffer, *, iostat=ios) GD%lonP
	           !print*, 'Read lonP: ', GD%lonP
	        case ('latp')
	           read(buffer, *, iostat=ios) GD%lonP
	           !print*, 'Read lonP: ', GD%lonP
!********************************************************************
	        case default
	           !print*, 'Skipping invalid label at line', line
	        end select
        end if
     end if
  end do

! Fix for the case when np is absent
if (GD%np == 0) then
	GD%np = GD%nu + GD%nv + GD%nh
end if

close(fh)

end subroutine

!**********************************************************************************************************

subroutine sht_smooth_H(lon, lat, H, P, lon_sht, lat_sht, H_smoothed, topo_dir, topo_file)

!% Given H on a (phi,th) grid, calculates H_smoothed

implicit none

     real(dp), intent(in)		:: lat(:), lon(:), H(:,:)
     character(len=*), intent(in)  :: topo_dir, topo_file
     type(params), intent(in)	:: P

     real(dp), allocatable	:: H_smoothed(:,:), smooth(:)
     integer			    :: nlon, nlat

     real(dp), allocatable	:: lon_sht(:), lat_sht(:)
     integer			    :: nlon_sht, nlat_sht
     real(dp)               :: dlon_sht, dlat_sht

     integer            :: ntrunc, istat, j, deln
     character(len = *), parameter :: gridtype = 'REG'	! (regularly spaced grid)
!     real    ::      T1, T2
	 logical, parameter :: save_coeffs = .true. ! choose between a slower method that saves SHT coeffs for the future usage and quicker method that doesn't

! shortcuts
nlon = size(lon)
nlat = size(lat)

if ( (nlon /= size(H, 1)).or.(nlat /= size(H, 2)) ) then
	print *, 'Dimensions of lon, lat and H do not agree'
	stop
end if

!***************************************************
! On a regular grid spanning -pi/2<th<pi/2 and 0<ph<2*pi
! the condition (ntrunc <= nlat-1) is equiv. to (2*ntrunc <= nlon)
	ntrunc = nlon/2

!***************************************************
! Construct a regular spaced lon-lat grid for the SHT routine
! as the basis take the lon partition of [0, 2*pi) interval (lon_sht = lon)
! Lat coord lat_sht must span [-pi/2, pi/2] and must start at the north pole

nlon_sht = 2*ntrunc
dlon_sht = 2*pi/nlon_sht
nlat_sht = nlon_sht/2 + 1
dlat_sht = pi/(nlat_sht-1)

if (dlon_sht /= dlat_sht) then
	print *, 'dlon_sht /= dlat_sht in SHT subroutine calc_hsal_gridded'
	stop
end if
!print *, "hi1", nlon, nlon_sht, nlat, nlat_sht

! initialises the sht grid structure
    allocate(lon_sht(nlon_sht), stat = istat)
    lon_sht = (/ ( (2*pi/nlon_sht*j), j=0,nlon_sht-1 ) /)	! the grid for the SHT routine spans longitude points
! Th must be the colatitude needed by the SHT routine (increases from the North pole southwards)
    allocate(lat_sht(nlat_sht), stat = istat)
    lat_sht = (/ ( (pi/2 - pi/(nlat_sht-1)*j), j=0,nlat_sht-1 ) /)

!print *, "hi2", lon_sht(1), lon_sht(nlon_sht), lon(1), lon(nlon_sht)
!print *, "hi3", lat_sht(1), lat_sht(nlat_sht), lat(1), lat(nlat_sht)

! Check if the original array satisfies the requirements or if interpolation has to be applied
if ( (lon_sht(1) == lon(1)).and.(lon_sht(nlon_sht) == lon(nlon_sht)) ) then
	if ( (lat_sht(1) == lat(1)).and.(lat_sht(nlat_sht) == lat(nlat_sht)) ) then
		allocate(H_smoothed(nlon_sht, nlat_sht), stat = istat)
		H_smoothed(:,:) = H(1:nlon_sht, 1:nlat_sht)
	elseif ( (lat_sht(1) == lat(nlat_sht)).and.(lat_sht(nlat_sht) == lat(1)) ) then
		allocate(H_smoothed(nlon_sht, nlat_sht), stat = istat)
		H_smoothed(:,:) = H(1:nlon_sht, nlat_sht:1:-1)
		print *, "the right thing"
	else
!		print *, "hi4"
		call bilinear_2d(nlon, nlat, lon, lat, H, nlon_sht, nlat_sht, lon_sht, lat_sht, H_smoothed, P%lib)
	endif
else
! otherwise interpolate H to the grid
    call bilinear_2d(nlon, nlat, lon, lat, H, nlon_sht, nlat_sht, lon_sht, lat_sht, H_smoothed, P%lib)
endif

!	call save_topo(topo_dir // topo_file, nlon_sht, nlat_sht, lon_sht, lat_sht, H_smoothed)

!Smoothing as a function of total wavenumber
allocate(smooth(nlat_sht), stat = istat)
smooth = 0.

if ((P%sht_smooth_H == 0).or.(P%sht_smooth_H >= nlat_sht)) then
		deln = nlat_sht
else
		deln = P%sht_smooth_H + 1
endif

if (P%smooth_type==1) then
	do j = 1, deln
		smooth(j) = 1. ! straightforward truncation smoothing;
	enddo
elseif (P%smooth_type==2) then
	do j = 1, nlat_sht
		smooth(j) = exp(-(float(j-1)/deln)**2)!Gaussian spectral smoothing, where deln is the e-folding width of the smoother in total wavenumber space.
	enddo
else
	print *, "Unexpected SHT smoothing type. Error in function generate_global_grid."
	stop
endif

! Do the main calculation to find smoothed H on the grid
!	print *, maxval(H_smoothed),minval(H_smoothed)
if (save_coeffs) then
	call specsmooth(H_smoothed,smooth,gridtype, (P%messages > 0), topo_dir//'sht_coeff_'//topo_file) ! spec_coeff_file
else
	call sht_project(H_smoothed,P%sht_smooth_H)
endif
!	print *, maxval(H_smoothed),minval(H_smoothed)

! Clean up
	call cleanup(gridtype)

if (P%smooth_type==1) then
   	call save_topo(topo_dir//'smoothed_t'//tostring(P%smooth_type)//'_trunc'//tostring(deln-1)// '_'// &
   	               topo_file, nlon_sht, nlat_sht, lon_sht, lat_sht, H_smoothed)
endif

end subroutine sht_smooth_H

!********************************************************************
!-------------------- OLD VERSION, OBSOLETE -------------------------
!********************************************************************
    subroutine interp_H_grid_old(numLons, numLats, xValues, yValues, zValues, &
        numLonsI, numLatsI, xValuesI, yValuesI, zValuesI)


     integer ::      status
     integer ::      i, j
     integer, allocatable :: indi(:), indj(:)
     real(dp)  ::      dx, dy!, dxI, dyI

     integer     ::     numLons, numLats
     real(dp), allocatable, intent(in) :: xValues(:), yValues(:)
     integer, allocatable, intent(in) :: zValues(:, :)

     real(dp)     ::     weight1, weight2
     real(dp)     ::     ztmp1, ztmp2
     integer     ::     numLonsI, numLatsI
     real(dp), allocatable, intent(in) :: xValuesI(:), yValuesI(:)
     integer, allocatable :: zValuesI(:, :)


write(*, '("Interpolating topo to ", f5.2, " min resolution global grid: ")', advance = 'no') 360*60/real(numLonsI)

If (numLons <= numLonsI .or. numLats <= numLatsI .or. numLonsI <= 0 .or. numLatsI <= 0 ) then
     Write (*, *) 'Bad Dimenisons...subroutine interp_H_grid 1'
     print *, numLons, numLonsI, " and ", numLats, numLatsI
     stop
end if

If (xValues(1) >= xValues(numLons) .or. yValues(1) >= yValues(numLats) ) then
     Write (*, *) 'Bad Inputs...subroutine interp_H_grid 2'
     stop
end if

if (allocated(ZValuesI)) deallocate(ZValuesI)
    allocate(ZValuesI(numLonsI, numLatsI), stat = status)
!    if (allocated(XValuesI)) deallocate(XValuesI)
!    allocate(XValuesI(numLonsI), stat = statusx)
!    if (allocated(YValuesI)) deallocate(YValuesI)
!    allocate(YValuesI(numLatsI), stat = statusy)

!    if (statusx /= 0 .or. statusy /= 0 .or. statusz /= 0) then
    if (status /= 0) then
         print *, "[Topo] Can't allocate the memory."
        stop
    end if

    dx = (xValues(numLons) - xValues(1))/(numLons - 1) !the topography is periodical, i.e. map from -180 to 180
    dy = (yValues(numLats) - yValues(1))/(numLats - 1) !and data for xValues(1) and xValues(numLons) is the same

!    dxI = (xValues(numLons) - xValues(1))/numLonsI !grid points are NOT on the -180~180 boundary
!    dyI = (yValues(numLats) - yValues(1))/(numLatsI - 1)

!XValuesI = (/ ( (xValues(1) + dxI*(i-1)), i=1,numLonsI ) /)
!yValuesI = (/ ( (yValues(1) + dyI*(j-1)), j=1,numLatsI ) /)

!Index to know which points of the original data use for interpolation
allocate(indi(numLonsI), stat = status)
allocate(indj(numLatsI), stat = status)

do i = 1, numLonsI
indi(i) = (xValuesI(i)-xValues(1))/dx + 1
enddo
do j = 1, numLatsI
indj(j) = (yValuesI(j)-yValues(1))/dy + 1
enddo

    do i = 1, numLonsI
         do j = 1, numLatsI
         !     First horizontal averaging
          weight1 = (xValuesI(i) - xValues(indi(i)))/dx
          weight2 = (xValues(indi(i)+1) - xValuesI(i))/dx
          if ( (weight1<0 - eps) .or. (weight1>1 + eps) .or. (weight2<0 - eps) .or. (weight2>1 + eps) )  then
                call handle_av_err('Incorrect averaging 1', real(weight1,kind=dp), real(weight2,kind=dp));
                stop
          end if

          ztmp1 = ZValues(indi(i) + 1, indj(j))*weight1 + &
                              ZValues(indi(i), indj(j))*weight2
          ztmp2 = ZValues(indi(i) + 1, indj(j)+1)*weight1 + &
                              ZValues(indi(i), indj(j)+1)*weight2

         !     Then vertical
          weight1 = (yValuesI(j) - yValues(indj(j)))/dy
          weight2 = (yValues(indj(j)+1) - yValuesI(j))/dy
          if (weight1<0 - eps .or. weight1>1 + eps .or. weight2<0 - eps .or. weight2>1 + eps )  then
                call handle_av_err('Incorrect averaging 2', real(weight1,kind=dp), real(weight2,kind=dp));
                stop
          end if

          ZValuesI(i, j) = ztmp2*weight1 + ztmp1*weight2
         enddo
    enddo

        print *, "done."

end subroutine
!***********************************************************************************************************

end module generate_grid



