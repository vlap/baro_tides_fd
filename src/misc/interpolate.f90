module interpolate

     use precisions, only: sp, dp, wp, cwp
     use my_trigs
     use my_sparse
     use save_load
!     use my_sparse_aggregate

!==========================================================================================
!**************************** INTERFACES **************************************************
!==========================================================================================
!------------------- BILINEAR INTERPOLATION (GRID TO GRID)-------------------------------
     interface bilinear_2d
          module procedure bilinear_grid2grid_4
          module procedure bilinear_grid2grid_8
          module procedure bilinear_grid2grid_cmplx
     end interface

     contains

!====================================================================
!---------- INTERPOLATEION ON A SPHERE (RANDOM TO GRID) -----------
!====================================================================
!subroutine save_sol4plot_cmplx(nph, nth, nvar, varp, var, sol4plot_name)
!
!    implicit none
!
!     integer, intent(in) 		:: nph, nth, nvar
!     integer, intent(in) 		:: varp(:, :)
!     complex(cwp), intent(in)	:: var(:)
!
!    character(len=*) :: sol4plot_name
!
!      OPEN(10,status='unknown', file=sol4plot_name, form='unformatted', action='write', access='stream')
!
!      write(10) nph, nth, nvar
!      write(10) varp(:,1)
!      write(10) varp(:,2)
!      write(10) abs(var)
!      CLOSE(10)
!
!
!end subroutine save_sol4plot_cmplx
!====================================================================
!------- BILINEAR INTERPOLATION COMPLEX (GRID TO GRID) --------------
!====================================================================
subroutine bilinear_grid2grid_cmplx(nph, nth, ph, th, f, nph1, nth1, ph1, th1, f1, lib)!, dir_grid)
! (use simple bilinear interpolation)

!% given a function f on a {phi, theta} grid,
!% with -pi/2 <= th <= pi/2 and 0 <= ph < 2*pi,
!% interpolates to the grid {phi1, theta1}. The
!% output grid is specified by the two vectors
!% (th1,ph1).
!% Assumes that th, ph, th1, ph1 are all increasing,
!% with possibly irregular spacing.
    implicit none
!character(len=*) :: dir_grid

integer, intent(in)  :: nph, nth, nph1, nth1
complex(cwp),intent(in)  :: f(nph, nth) ! old data
real(wp),intent(in)  :: ph(nph), th(nth), ph1(nph1), th1(nth1)
character, intent(in)	:: lib

complex(wp),allocatable :: f1(:,:) ! interpolated data

type (triplet)	     :: im_ph, im_th ! interpolation matrix (common for x and y)
complex(cwp), allocatable:: vec_tmp(:), mat_tmp(:, :)

integer:: c, loc, istatus

! Prepare interpolation matrices
!print *, nph, nth, nph1, nth1

!**********************************************************************************************************
!% given a grid th, generates a matrix suitable
!% for interpolating to the grid th1.
call init_sparse_0(im_th, nth1, nth, 2*nth1)
do c = 1,nth1
  ! find neighbors of th1(c) in th grid
  call locate(nth, th, th1(c), loc)
!  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=1
	im_th%vals(2*c-1)=1
  elseif (loc == nth) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=nth
	im_th%vals(2*c-1)=1
  else
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=loc
	im_th%vals(2*c-1)=(th(loc+1)-th1(c))/(th(loc+1)-th(loc))

	im_th%indi(2*c)=c
	im_th%indj(2*c)=loc+1
	im_th%vals(2*c)=(th1(c)-th(loc))/(th(loc+1)-th(loc))
  end if
!  print *, im_th%vals(2*c-1),im_th%vals(2*c)
end do

!print *, ""

!% given a grid ph with 0 <= ph < 2*pi,
!% generates a matrix suitable for interpolating
!% to the grid ph1, with 0 <= ph1 < 2*pi.
! Takes into account periodicity in phi
call init_sparse_0(im_ph, nph1, nph, 2*nph1)
do c = 1,nph1

  ! find neighbors of ph1(c) in ph grid
  call locate(nph, ph, ph1(c), loc)
!  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=1

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=nph

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph1(c)-(ph(nph)-2*pi))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph(1)-ph1(c))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  elseif (loc == nph) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=nph

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=1

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph(1)+2*pi-ph1(c))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph1(c)-ph(nph))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  else
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=loc
	im_ph%vals(2*c-1)=(ph(loc+1)-ph1(c))/(ph(loc+1)-ph(loc))

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=loc+1
	im_ph%vals(2*c)=(ph1(c)-ph(loc))/(ph(loc+1)-ph(loc))
  end if
!  print *, im_ph%vals(2*c-1),im_ph%vals(2*c)
end do
! first average in ph
    allocate(vec_tmp(nph1), stat = istatus)
    allocate(mat_tmp(nth, nph1), stat = istatus)

do c = 1, nth
	call coo_vec_mul(im_ph,f(:, c),lib, vec_tmp)
	mat_tmp(c, :) = vec_tmp
end do

! then average in th
    if (allocated(vec_tmp)) deallocate(vec_tmp)
    allocate(vec_tmp(nth1), stat = istatus)

    if (allocated(f1)) deallocate(f1)
    allocate(f1(nph1, nth1), stat = istatus)

do c = 1, nph1
	call coo_vec_mul(im_th,mat_tmp(:, c),lib, vec_tmp)
	f1(c, :) = vec_tmp
end do

!call save_cmat(im_ph, dir_grid // 'im_ph.dat')
!call save_cmat(im_th, dir_grid // 'im_th.dat')

end subroutine bilinear_grid2grid_cmplx

!====================================================================
!---------- BILINEAR INTERPOLATION (GRID TO GRID) ------------------
!====================================================================
subroutine bilinear_grid2grid_4(nph, nth, ph, th, f, nph1, nth1, ph1, th1, f1,lib)!, dir_grid)
! (use simple bilinear interpolation)

!% given a function f on a {phi, theta} grid,
!% with -pi/2 <= th <= pi/2 and 0 <= ph < 2*pi,
!% interpolates to the grid {phi1, theta1}. The
!% output grid is specified by the two vectors
!% (th1,ph1).
!% Assumes that th, ph, th1, ph1 are all increasing,
!% with possibly irregular spacing.
    implicit none
!character(len=*) :: dir_grid

integer, intent(in)  :: nph, nth, nph1, nth1
real(sp),intent(in)  :: f(nph, nth) ! old data
real(sp),intent(in)  :: ph(nph), th(nth), ph1(nph1), th1(nth1)
character, intent(in)	:: lib

real(sp),allocatable :: f1(:,:) ! interpolated data

type (triplet_sp)	     :: im_ph, im_th ! interpolation matrix (common for x and y)
real(sp), allocatable:: vec_tmp(:), mat_tmp(:, :)

integer:: c, loc, istatus

! Prepare interpolation matrices
!print *, nph, nth, nph1, nth1
!**********************************************************************************************************
!% given a grid th, generates a matrix suitable
!% for interpolating to the grid th1.

call init_sparse_0(im_th, nth1, nth, 2*nth1)

do c = 1,nth1
  ! find neighbors of th1(c) in th grid
  call locate_4(nth, th, th1(c), loc)
!  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=1
	im_th%vals(2*c-1)=1
  elseif (loc == nth) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=nth
	im_th%vals(2*c-1)=1
  else
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=loc
	im_th%vals(2*c-1)=(th(loc+1)-th1(c))/(th(loc+1)-th(loc))

	im_th%indi(2*c)=c
	im_th%indj(2*c)=loc+1
	im_th%vals(2*c)=(th1(c)-th(loc))/(th(loc+1)-th(loc))
  end if
!  print *, im_th%vals(2*c-1),im_th%vals(2*c)
end do

!% given a grid ph with 0 <= ph < 2*pi,
!% generates a matrix suitable for interpolating
!% to the grid ph1, with 0 <= ph1 < 2*pi.
! Takes into account periodicity in phi
call init_sparse_0(im_ph, nph1, nph, 2*nph1)
do c = 1,nph1

  ! find neighbors of ph1(c) in ph grid
  call locate_4(nph, ph, ph1(c), loc)
!  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=1

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=nph

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph1(c)-(ph(nph)-2*pi))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph(1)-ph1(c))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  elseif (loc == nph) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=nph

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=1

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph(1)+2*pi-ph1(c))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph1(c)-ph(nph))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  else
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=loc
	im_ph%vals(2*c-1)=(ph(loc+1)-ph1(c))/(ph(loc+1)-ph(loc))

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=loc+1
	im_ph%vals(2*c)=(ph1(c)-ph(loc))/(ph(loc+1)-ph(loc))
  end if
!  print *, im_ph%vals(2*c-1),im_ph%vals(2*c)
end do
! first average in ph
    allocate(vec_tmp(nph1), stat = istatus)
    allocate(mat_tmp(nth, nph1), stat = istatus)

do c = 1, nth
	call coo_vec_mul(im_ph,f(:, c),lib, vec_tmp)
	mat_tmp(c, :) = vec_tmp
end do

! then average in th
    if (allocated(vec_tmp)) deallocate(vec_tmp)
    allocate(vec_tmp(nth1), stat = istatus)

    if (allocated(f1)) deallocate(f1)
    allocate(f1(nph1, nth1), stat = istatus)

do c = 1, nph1
	call coo_vec_mul(im_th,mat_tmp(:, c),lib, vec_tmp)
	f1(c, :) = vec_tmp
end do

!call save_cmat(im_ph, dir_grid // 'im_ph.dat')
!call save_cmat(im_th, dir_grid // 'im_th.dat')

end subroutine bilinear_grid2grid_4

subroutine bilinear_grid2grid_8(nph, nth, ph, th, f, nph1, nth1, ph1, th1, f1, lib)!, dir_grid)
! (use simple bilinear interpolation)

!% given a function f on a {phi, theta} grid,
!% with -pi/2 <= th <= pi/2 and 0 <= ph < 2*pi,
!% interpolates to the grid {phi1, theta1}. The
!% output grid is specified by the two vectors
!% (th1,ph1).
!% Assumes that th, ph, th1, ph1 are all increasing,
!% with possibly irregular spacing.
    implicit none
!character(len=*) :: dir_grid

integer, intent(in)  :: nph, nth, nph1, nth1
real(dp),intent(in)  :: f(nph, nth) ! old data
real(dp),intent(in)  :: ph(nph), th(nth), ph1(nph1), th1(nth1)
character, intent(in)	:: lib

real(dp),allocatable :: f1(:,:) ! interpolated data

type (triplet_dp)	     :: im_ph, im_th ! interpolation matrix (common for x and y)
real(dp), allocatable:: vec_tmp(:), mat_tmp(:, :)

integer:: c, loc, istatus

! Prepare interpolation matrices
!print *, nph, nth, nph1, nth1
!**********************************************************************************************************
!% given a grid th, generates a matrix suitable
!% for interpolating to the grid th1.
call init_sparse_0(im_th, nth1, nth, 2*nth1)

do c = 1,nth1
  ! find neighbors of th1(c) in th grid
  call locate_8(nth, th, th1(c), loc)
!  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=1
	im_th%vals(2*c-1)=1
  elseif (loc == nth) then
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=nth
	im_th%vals(2*c-1)=1
  else
	im_th%indi(2*c-1)=c
	im_th%indj(2*c-1)=loc
	im_th%vals(2*c-1)=(th(loc+1)-th1(c))/(th(loc+1)-th(loc))

	im_th%indi(2*c)=c
	im_th%indj(2*c)=loc+1
	im_th%vals(2*c)=(th1(c)-th(loc))/(th(loc+1)-th(loc))
  end if
!  print *, im_th%vals(2*c-1),im_th%vals(2*c)
end do

!% given a grid ph with 0 <= ph < 2*pi,
!% generates a matrix suitable for interpolating
!% to the grid ph1, with 0 <= ph1 < 2*pi.
! Takes into account periodicity in phi
call init_sparse_0(im_ph, nph1, nph, 2*nph1)
do c = 1,nph1

  ! find neighbors of ph1(c) in ph grid
  call locate_8(nph, ph, ph1(c), loc)
!	  write(*, '(i3)', advance = 'no') loc
  ! if not between two points use just the closest neighbour
  if (loc == 0) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=1

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=nph

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph1(c)-(ph(nph)-2*pi))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph(1)-ph1(c))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  elseif (loc == nph) then
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=nph

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=1

	if (2*pi-(ph(nph)-ph(1)) > 0) then
		im_ph%vals(2*c-1)=(ph(1)+2*pi-ph1(c))/(2*pi-(ph(nph)-ph(1)))
		im_ph%vals(2*c)=(ph1(c)-ph(nph))/(2*pi-(ph(nph)-ph(1)))
	else
		im_ph%vals(2*c-1)=1
		im_ph%vals(2*c)=0
	endif
  else
	im_ph%indi(2*c-1)=c
	im_ph%indj(2*c-1)=loc
	im_ph%vals(2*c-1)=(ph(loc+1)-ph1(c))/(ph(loc+1)-ph(loc))

	im_ph%indi(2*c)=c
	im_ph%indj(2*c)=loc+1
	im_ph%vals(2*c)=(ph1(c)-ph(loc))/(ph(loc+1)-ph(loc))
  end if
!  print *, im_ph%vals(2*c-1),im_ph%vals(2*c)
end do
! first average in ph
    allocate(vec_tmp(nph1), stat = istatus)
    allocate(mat_tmp(nth, nph1), stat = istatus)

do c = 1, nth
	call coo_vec_mul(im_ph,f(:, c),lib, vec_tmp)
	mat_tmp(c, :) = vec_tmp
end do

! then average in th
    if (allocated(vec_tmp)) deallocate(vec_tmp)
    allocate(vec_tmp(nth1), stat = istatus)

    if (allocated(f1)) deallocate(f1)
    allocate(f1(nph1, nth1), stat = istatus)

do c = 1, nph1
	call coo_vec_mul(im_th,mat_tmp(:, c),lib, vec_tmp)
	f1(c, :) = vec_tmp
end do

!call save_cmat(im_ph, dir_grid // 'im_ph.dat')
!call save_cmat(im_th, dir_grid // 'im_th.dat')

end subroutine bilinear_grid2grid_8
!**********************************************************************************************************
subroutine locate(n, xx, x, j)
!% given an array xx(1:n), and a value x, returns a value
!% j such that x is between xx(j) and xx(j+1). xx(1:n) must
!% be monotonic, either increasing or decreasing. j=0 or j=n
!% is returned to indicate that x is out of range.
    implicit none

integer, intent(in)  :: n
real(wp),intent(in)  :: xx(n), x

integer :: j, jl, ju, jm

! First consider the out of range cases
if ( (xx(n)>=xx(1)).and.(x <= xx(1)) .or. (xx(n)<=xx(1)).and.(x >= xx(1)) ) then
  j=0
elseif ( (xx(n)>=xx(1)).and.(x >= xx(n)) .or. (xx(n)<=xx(1)).and.(x <= xx(n)) ) then
  j=n
else

	jl=0
	ju=n+1

	do while (ju-jl > 1)
	  jm =(ju+jl)/2 ! implicit floor rounding
	  if ( (xx(n)>=xx(1)) .eqv. (x>=xx(jm)) ) then
	     jl=jm
	  else
	    ju=jm
	  end if
	end do

	j = jl

end if
!
!if (j == 0) then
!print *, x
!endif

end subroutine locate

subroutine locate_4(n, xx, x, j)

    implicit none

integer, intent(in)  :: n
real(sp),intent(in)  :: xx(n), x

integer :: j, jl, ju, jm

! First consider the out of range cases
if ( (xx(n)>=xx(1)).and.(x <= xx(1)) .or. (xx(n)<=xx(1)).and.(x >= xx(1)) ) then
  j=0
elseif ( (xx(n)>=xx(1)).and.(x >= xx(n)) .or. (xx(n)<=xx(1)).and.(x <= xx(n)) ) then
  j=n
else

	jl=0
	ju=n+1

	do while (ju-jl > 1)
	  jm =(ju+jl)/2 ! implicit floor rounding
	  if ( (xx(n)>=xx(1)) .eqv. (x>=xx(jm)) ) then
	     jl=jm
	  else
	    ju=jm
	  end if
	end do

	j = jl

end if
!
!if (j == 0) then
!print *, x
!endif

end subroutine locate_4

subroutine locate_8(n, xx, x, j)

    implicit none

integer, intent(in)  :: n
real(dp),intent(in)  :: xx(n), x

integer :: j, jl, ju, jm

! First consider the out of range cases
if ( (xx(n)>=xx(1)).and.(x <= xx(1)) .or. (xx(n)<=xx(1)).and.(x >= xx(1)) ) then
  j=0
elseif ( (xx(n)>=xx(1)).and.(x >= xx(n)) .or. (xx(n)<=xx(1)).and.(x <= xx(n)) ) then
  j=n
else

	jl=0
	ju=n+1

	do while (ju-jl > 1)
	  jm =(ju+jl)/2 ! implicit floor rounding
	  if ( (xx(n)>=xx(1)) .eqv. (x>=xx(jm)) ) then
	     jl=jm
	  else
	    ju=jm
	  end if
	end do

	j = jl

end if
!
!if (j == 0) then
!print *, x
!endif

end subroutine locate_8
!**********************************************************************************************************
!**********************************************************************************************************

end module interpolate
