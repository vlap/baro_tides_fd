! --------------------------------------------------------------------
! MODULE  MyTrigonometricFunctions:
!    This module provides the following functions and constants
!    (1) r2d()     - converts its argument in radian to
!                               degree
!    (2) d2r()     - converts its argument in degree to
!                               radian
!    (3) MySIN()              - compute the sine of its argument in
!                               degree
!    (4) MyCOS()              - compute the cosine of its argument
!                               in degree
! --------------------------------------------------------------------

MODULE  my_trigs

   use precisions, only: sp, dp, wp, cwp

!==========================================================================================
!**************************** INTERFACES **************************************************
!==========================================================================================
!---------------------------- LOGICAL TO INTEGER --------------------------------------
     interface log2int
          module procedure log2int_4
          module procedure log2int_8
          module procedure log2int_vec_4
          module procedure log2int_vec_8
     end interface
!---------------------------- RADIAN TO DEGREE ---------------------------------------
     interface r2d
          module procedure RadianToDegree
          module procedure RadianToDegreeDP
          module procedure RadianToDegree_vec
          module procedure RadianToDegreeDP_vec
     end interface
!---------------------------- DEGREE TO RADIAN ---------------------------------------
     interface d2r
          module procedure DegreeToRadian
          module procedure DegreeToRadianDP
          module procedure DegreeToRadian_vec
          module procedure DegreeToRadianDP_vec
     end interface
!---------------------------- THETA TO TAU ---------------------------------------
     interface th2ta
          module procedure th2ta_4
          module procedure th2ta_8
          module procedure th2ta_vec_4
          module procedure th2ta_vec_8
     end interface
!---------------------------- TAU TO THETA ---------------------------------------
     interface ta2th
          module procedure ta2th_4
          module procedure ta2th_8
          module procedure ta2th_vec_4
          module procedure ta2th_vec_8
     end interface
!---------------------------- IMAGINARY PART --------------------------------------
     interface imag
          module procedure imag_my_4
          module procedure imag_my_8
          module procedure imag_my_vec_4
          module procedure imag_my_vec_8
     end interface
!---------------------------- SINC FUNCTION ---------------------------------------
     interface sinc
          module procedure sinc_4
          module procedure sinc_8
          module procedure sinc_vec_4
          module procedure sinc_vec_8
     end interface
!==========================================================================================

   complex(cwp), parameter :: i = dcmplx(0, 1)
   real(dp), PARAMETER :: PI        = 2*ACOS(0.)      ! some constants
   real(dp), PARAMETER :: Degree180 = 180.0
   real(dp), PARAMETER :: R_to_D    = Degree180/PI
   real(dp), PARAMETER :: D_to_R    = PI/Degree180

CONTAINS

! --------------------------------------------------------------------
! FUNCTIONS  log2int():
!    This function does the natural conversion of LOGICAL to INTEGER
! --------------------------------------------------------------------
    integer(4) FUNCTION log2int_4(log4)
      IMPLICIT  NONE
      logical(4), INTENT(IN):: log4

      if (log4) then
      	log2int_4 = 1
      else
      	log2int_4 = 0
      endif
   END FUNCTION log2int_4

    integer(8) FUNCTION log2int_8(log8)
      IMPLICIT  NONE
      logical(8), INTENT(IN):: log8

      if (log8) then
      	log2int_8 = 1
      else
      	log2int_8 = 0
      endif
   END FUNCTION log2int_8

   FUNCTION  log2int_vec_4(log4_vec)
      IMPLICIT  NONE
      logical(4), DIMENSION(:), INTENT(IN) :: log4_vec
      INTEGER ::n, istat
      integer(4), DIMENSION(:), allocatable :: log2int_vec_4
      n = size(log4_vec)

	allocate(log2int_vec_4(n), stat = istat)
	log2int_vec_4 = 0

	where (log4_vec) log2int_vec_4 = 1

   END FUNCTION  log2int_vec_4

   FUNCTION  log2int_vec_8(log8_vec)
      IMPLICIT  NONE
      logical(8), DIMENSION(:), INTENT(IN) :: log8_vec
      INTEGER ::n, istat
      integer(8), DIMENSION(:), allocatable :: log2int_vec_8
      n = size(log8_vec)

	allocate(log2int_vec_8(n), stat = istat)
	log2int_vec_8 = 0

	where (log8_vec) log2int_vec_8 = 1

   END FUNCTION  log2int_vec_8

! --------------------------------------------------------------------
! FUNCTIONS  imag_my():
!    This function takes a complex(4) or complex(8) argument and returns it's imaginary part in real(4) or real(8).
! --------------------------------------------------------------------
    real(4) FUNCTION imag_my_4(z4)
      IMPLICIT  NONE
      complex(4), INTENT(IN):: z4

      imag_my_4 = aimag(z4)

   END FUNCTION imag_my_4

    real(8) FUNCTION imag_my_8(z8)
      IMPLICIT  NONE
      complex(8), INTENT(IN):: z8

      imag_my_8 = dimag(z8)

   END FUNCTION imag_my_8

   FUNCTION  imag_my_vec_4(vec_4)
      IMPLICIT  NONE
      complex(4), DIMENSION(:), INTENT(IN) :: vec_4
      INTEGER ::n, istat
      real(4), DIMENSION(:), allocatable :: imag_my_vec_4
      n = size(vec_4)

	allocate(imag_my_vec_4(n), stat = istat)

      imag_my_vec_4 = aimag(vec_4)

   END FUNCTION  imag_my_vec_4

   FUNCTION  imag_my_vec_8(vec_8)
      IMPLICIT  NONE
      complex(8), DIMENSION(:), INTENT(IN) :: vec_8
      INTEGER ::n, istat
      real(8), DIMENSION(:), allocatable :: imag_my_vec_8
      n = size(vec_8)

	allocate(imag_my_vec_8(n), stat = istat)

      imag_my_vec_8 = dimag(vec_8)

   END FUNCTION  imag_my_vec_8
! --------------------------------------------------------------------
! FUNCTION  r2d():
!    This function takes a real(wp) argument in radian and converts it to
! the equivalent degree.
! --------------------------------------------------------------------
    real(4) FUNCTION RadianToDegree(Radian)
      IMPLICIT  NONE
      real(4), INTENT(IN) :: Radian

      RadianToDegree = Radian * R_to_D

   END FUNCTION RadianToDegree

    real(8) FUNCTION RadianToDegreeDP(Radian)
      IMPLICIT  NONE
      real(8), INTENT(IN) :: Radian

      RadianToDegreeDP = Radian * R_to_D

   END FUNCTION RadianToDegreeDP

    FUNCTION RadianToDegree_vec(Radian)
      IMPLICIT  NONE
      real(4), DIMENSION(:), INTENT(IN) :: Radian
      INTEGER ::n, istat
      real(4), DIMENSION(:), allocatable :: RadianToDegree_vec
      n = size(Radian)

	allocate(RadianToDegree_vec(n), stat = istat)

      RadianToDegree_vec = Radian * R_to_D

   END FUNCTION RadianToDegree_vec

    FUNCTION RadianToDegreeDP_vec(Radian)
      IMPLICIT  NONE
      real(8), DIMENSION(:), INTENT(IN) :: Radian
      INTEGER ::n, istat
      real(8), DIMENSION(:), allocatable :: RadianToDegreeDP_vec
      n = size(Radian)

	allocate(RadianToDegreeDP_vec(n), stat = istat)

      RadianToDegreeDP_vec = Radian * R_to_D

   END FUNCTION RadianToDegreeDP_vec
! --------------------------------------------------------------------
! FUNCTION  d2r():
!    This function takes a real(wp) argument in degree and converts it to
! the equivalent radian.
! --------------------------------------------------------------------
   real(4) FUNCTION  DegreeToRadian(Degree)
      IMPLICIT  NONE
      real(4), INTENT(IN) :: Degree

      DegreeToRadian = Degree * D_to_R
   END FUNCTION  DegreeToRadian

   real(8) FUNCTION  DegreeToRadianDP(Degree)
      IMPLICIT  NONE
      real(8), INTENT(IN) :: Degree

      DegreeToRadianDP = Degree * D_to_R
   END FUNCTION  DegreeToRadianDP

   FUNCTION  DegreeToRadian_vec(Degree)
      IMPLICIT  NONE
      real(4), DIMENSION(:), INTENT(IN) :: Degree
      INTEGER ::n, istat
      real(4), DIMENSION(:), allocatable :: DegreeToRadian_vec
      n = size(Degree)

	allocate(DegreeToRadian_vec(n), stat = istat)

      DegreeToRadian_vec = Degree * D_to_R
   END FUNCTION  DegreeToRadian_vec

   FUNCTION  DegreeToRadianDP_vec(Degree)
      IMPLICIT  NONE
      real(8), DIMENSION(:), INTENT(IN) :: Degree
      INTEGER ::n, istat
      real(8), DIMENSION(:), allocatable :: DegreeToRadianDP_vec
      n = size(Degree)

	allocate(DegreeToRadianDP_vec(n), stat = istat)

      DegreeToRadianDP_vec = Degree * D_to_R
   END FUNCTION  DegreeToRadianDP_vec


! --------------------------------------------------------------------
! FUNCTION  th2ta():
!    converts a value of the latitude th (in radians),
!    to a generalized coordinate tau.
! --------------------------------------------------------------------

   real(sp) FUNCTION  th2ta_4(th,coor)
      IMPLICIT  NONE
      real(sp), INTENT(IN) :: th
      integer, intent(in)  :: coor

	if (coor == 1) then
      th2ta_4 = atanh(sin(th))
    else
      th2ta_4 = th
    endif

   END FUNCTION  th2ta_4
   real(dp) FUNCTION  th2ta_8(th,coor)
      IMPLICIT  NONE
      real(dp), INTENT(IN) :: th
      integer, intent(in)  :: coor

	if (coor == 1) then
      th2ta_8 = atanh(sin(th))
    else
      th2ta_8 = th
    endif

   END FUNCTION  th2ta_8

   FUNCTION  th2ta_vec_4(th,coor)
      IMPLICIT  NONE
      real(sp), DIMENSION(:), INTENT(IN) :: th
      integer, intent(in)  :: coor
      INTEGER ::n, istat
      real(sp), DIMENSION(:), allocatable :: th2ta_vec_4
      n = size(th)

	allocate(th2ta_vec_4(n), stat = istat)

	if (coor == 1) then
      th2ta_vec_4 = atanh(sin(th))
    else
      th2ta_vec_4 = th
    endif

   END FUNCTION  th2ta_vec_4
   FUNCTION  th2ta_vec_8(th,coor)
      IMPLICIT  NONE
      real(dp), DIMENSION(:), INTENT(IN) :: th
      integer, intent(in)  :: coor
      INTEGER ::n, istat
      real(dp), DIMENSION(:), allocatable :: th2ta_vec_8
      n = size(th)

	allocate(th2ta_vec_8(n), stat = istat)

	if (coor == 1) then
      th2ta_vec_8 = atanh(sin(th))
    else
      th2ta_vec_8 = th
    endif

   END FUNCTION  th2ta_vec_8

! --------------------------------------------------------------------
! FUNCTION  ta2th():
!    This function takes a real(wp) argument in degree and computes its
! cosine value.  It does the computation by converting its argument to
! radian and uses Fortran's cos().
! --------------------------------------------------------------------

   real(sp) FUNCTION  ta2th_4(ta,coor)
      IMPLICIT  NONE
      real(sp), INTENT(IN) :: ta
      integer, intent(in)  :: coor

	if (coor == 1) then
      ta2th_4 = asin(tanh(ta))
    else
      ta2th_4 = ta
    endif

   END FUNCTION  ta2th_4
   real(dp) FUNCTION  ta2th_8(ta,coor)
      IMPLICIT  NONE
      real(dp), INTENT(IN) :: ta
      integer, intent(in)  :: coor

	if (coor == 1) then
      ta2th_8 = asin(tanh(ta))
    else
      ta2th_8 = ta
    endif

   END FUNCTION  ta2th_8

   FUNCTION  ta2th_vec_4(ta,coor)
      IMPLICIT  NONE
      real(sp), DIMENSION(:), INTENT(IN) :: ta
      integer, intent(in)  :: coor
      INTEGER ::n, istat
      real(sp), DIMENSION(:), allocatable :: ta2th_vec_4
      n = size(ta)

	allocate(ta2th_vec_4(n), stat = istat)


    if (coor == 1) then
      ta2th_vec_4 = asin(tanh(ta))
    else
      ta2th_vec_4 = ta
    endif

   END FUNCTION  ta2th_vec_4
   FUNCTION  ta2th_vec_8(ta,coor)
      IMPLICIT  NONE
      real(dp), DIMENSION(:), INTENT(IN) :: ta
      integer, intent(in)  :: coor
      INTEGER ::n, istat
      real(dp), DIMENSION(:), allocatable :: ta2th_vec_8
      n = size(ta)

	allocate(ta2th_vec_8(n), stat = istat)

	if (coor == 1) then
      ta2th_vec_8 = asin(tanh(ta))
    else
      ta2th_vec_8 = ta
    endif

   END FUNCTION  ta2th_vec_8

! --------------------------------------------------------------------
! FUNCTION  sinc():
!    The normalized sinc function is the Fourier transform of the rectangular function with no scaling.
!		sinc(x) = sin(pi*x)/(pi*x)
! When coeffs of a partial Fourier series (n<N) are premultiplied by sinc(k/N) ost of the Gibbs phenomenon
! is eliminated. In coordinate space is equivalent to integral-averaging of the partial Fourier series over at point x
! over the interval [x-pi/n, x+pi/n]. sinc(k/N) is known as the Lanczos σ factor
! --------------------------------------------------------------------

   real(sp) FUNCTION  sinc_4(x)
      IMPLICIT  NONE
      real(sp), INTENT(IN) :: x


	if (x == 0) then
      sinc_4 = 1
    else
      sinc_4 = sin(pi*x)/(pi*x)
    endif

   END FUNCTION  sinc_4
   real(dp) FUNCTION  sinc_8(x)
      IMPLICIT  NONE
      real(dp), INTENT(IN) :: x


	if (x == 0) then
      sinc_8 = 1
    else
      sinc_8 = sin(pi*x)/(pi*x)
    endif

   END FUNCTION  sinc_8

   FUNCTION  sinc_vec_4(x)
      IMPLICIT  NONE
      real(sp), DIMENSION(:), INTENT(IN) :: x

      INTEGER ::n, istat, j
      real(sp), DIMENSION(:), allocatable :: sinc_vec_4
      n = size(x)

	allocate(sinc_vec_4(n), stat = istat)

	do j = 1,n
		sinc_vec_4(j) = sinc_4(x(j))
	enddo

   END FUNCTION  sinc_vec_4
   FUNCTION  sinc_vec_8(x)
      IMPLICIT  NONE
      real(dp), DIMENSION(:), INTENT(IN) :: x

      INTEGER ::n, istat, j
      real(dp), DIMENSION(:), allocatable :: sinc_vec_8
      n = size(x)

	allocate(sinc_vec_8(n), stat = istat)

	do j = 1,n
		sinc_vec_8(j) = sinc_8(x(j))
	enddo

   END FUNCTION  sinc_vec_8

! --------------------------------------------------------------------
! Subroutine calc_latlon_nonrot
!% transforming from a computational lat-lon grid to a standard
!% lat-lon grid. The computational grid is shifted by ph0 in lon,
!% and dropped down th0 from the pole (in the direction of ph0).
!
!% In the new grid, thc runs from -pi/2 to pi/2
!%                  phc runs from 0 to 2*pi.
!% Inverse formulas yield th from -pi/2 to pi/2
!%                    and ph from 0 to 2*pi.
! --------------------------------------------------------------------
subroutine calc_latlon_nonrot(phc,thc,ph0,th0, ph_nonrot, th_nonrot)

	implicit none

     real(wp), allocatable	:: thc(:), phc(:)
     real(wp), allocatable	:: ph_nonrot(:), th_nonrot(:)
     real(wp)	:: ph0, th0
     integer	:: istat, n

	if (size(thc) /= size(phc)) then
   		print *, "Error in calc_latlon_nonrot"
   		stop
   	end if

   	n = size(phc)

	if (.not. allocated(th_nonrot)) allocate(th_nonrot(n), stat = istat)
	if (.not. allocated(ph_nonrot)) allocate(ph_nonrot(n), stat = istat)

	th_nonrot=asin(-sin(th0)*cos(thc)*cos(phc)+cos(th0)*sin(thc))
	!ph=ph0-real(i*log((cos(th0)*cos(thc)*cos(phc)+sin(th0)*sin(thc)+i*cos(thc)*sin(phc))/cos(th)));
	ph_nonrot=ph0+angle((cos(th0)*cos(thc)*cos(phc)+sin(th0)*sin(thc)+i*cos(thc)*sin(phc))/cos(th_nonrot))
	ph_nonrot = modulo(ph_nonrot,2*pi)

end subroutine calc_latlon_nonrot

! --------------------------------------------------------------------
! FUNCTION  ANGLE():
!%   ANGLE(H) returns the phase angles, in radians, of a matrix with
!%   complex elements.
!%
!%   Class support for input X:
!%      float: double, single
!%
!%   See also ABS, UNWRAP.
!
!% Clever way:
!% p = imag(log(h));
! --------------------------------------------------------------------
   FUNCTION  angle(z)

      IMPLICIT  NONE

      complex(cwp), DIMENSION(:), INTENT(IN) :: z
      INTEGER ::n, istat
      real(wp), DIMENSION(:), allocatable	 :: angle
      n = size(z)

	allocate(angle(n), stat = istat)
      angle = imag(log(z))

   END FUNCTION  angle

END MODULE  my_trigs
