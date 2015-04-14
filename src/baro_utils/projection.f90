module baro_projection

     use precisions, only: wp, cwp
     use my_trigs
     use save_load
     use my_sparse
     use my_sparse_aggregate
     use control
     use sal
     use my_blas

     integer, parameter :: sind = 86400 ! number of seconds in a day
     integer, parameter :: ppp = 8 ! the default value for P%nppp when P%nppp is not specified

contains

!**********************************************************************************************************
!**********************************************************************************************************

subroutine calc_dhat(cpts, nu, nv, P, dir_cols, dir_mats, dir_sols)

!% calculates, and writes to file, quantities related to the
!% projection of the bottom drag onto the component frequencies.

implicit none

     integer, intent(in)      :: nu, nv
     type(params)		:: P
     character(len=*)   :: cpts
     character(len=2)   :: cpt
     integer            :: ppc, coor, cooruv
     real(wp)			:: cdg, tmp
     real    ::      T1, T2 ! for measuring CPU (NOT REAL TIME!)
     integer :: wall_t1, wall_t2, clock_rate, clock_max

     integer            :: istat, ccpt, ct, ncpts, m,n

     type (triplet)       :: u2v, v2u
!     type (tide_params)   :: pars_tmp
     real(wp),allocatable :: omega0(:)
     real(wp)             :: t, dt, tmaxd
     integer, allocatable :: ncycles(:)

     real(wp), allocatable:: ta_u(:),ta_v(:)
     real(wp), allocatable :: H_u(:), H_v(:)

     character(len = *) :: dir_cols, dir_mats, dir_sols

     complex(cwp), allocatable :: u(:), v(:)!, h(:)
     complex(cwp), allocatable :: ustore_u(:,:), ustore_v(:,:), vstore_u(:,:), vstore_v(:,:)
     complex(cwp), allocatable :: u_u(:), u_v(:), v_u(:), v_v(:)
     complex(cwp), allocatable :: Iu(:,:), Iv(:,:)
     complex(cwp), allocatable :: Du(:), Dv(:)

!%==============
!% short-cuts
!%==============
cpts = P%cpts
coor = P%coor
cdg = P%cd

!%====================================================================
!% differentiate between lat-lon/MERC, velocity/transport formulations
!%====================================================================
if (coor == 1) then ! MERCATOR
!  if (P%uvform == 1) then
!    cooruv=11 ! VELOCITY
!  elseif (Puvform == 2)
    cooruv=12 ! TRANSPORT
!  endif
elseif (coor == 2) then ! LAT/LON
!  if (P%uvform == 1) then
!    cooruv=21 ! VELOCITY
!  elseif (Puvform == 2)
    cooruv=22 ! TRANSPORT
!  endif
endif
!==========================================================================================
write(*, '("  Welcome to calc_dhat: ")', advance='no')

ncpts=len(cpts)/2

!% number of time discretization points per tidal period
if (P%nppp == 0) then
    ppc=ppp ! default is 8
else
    ppc = P%nppp
endif

!==========================================================================================
!%=============================================
!% set up arrays (tides, projection integrals)
!%=============================================

    allocate(u(nu), v(nv), stat = istat)
    allocate(ustore_u(nu,ncpts), ustore_v(nv,ncpts), vstore_u(nu,ncpts), vstore_v(nv,ncpts), stat = istat)
    allocate(u_v(nv), v_u(nu), stat = istat)

      call load_alloc_sparse(u2v, dir_mats // 'u2v.dat')
      call load_alloc_sparse(v2u, dir_mats // 'v2u.dat')

do ccpt = 1, ncpts

     cpt=cpts(2*ccpt-1:2*ccpt)

          call load_vector(u, dir_sols // cpt // '_u' // '_m0' // '.dat')
          call load_vector(v, dir_sols // cpt // '_v' // '_m0' // '.dat')
  ustore_u(:,ccpt)=u
  vstore_v(:,ccpt)=v

     call coo_vec_mul(u2v, u, P%lib, u_v)
     call coo_vec_mul(v2u, v, P%lib, v_u)
  ustore_v(:,ccpt)=u_v
  vstore_u(:,ccpt)=v_u

enddo
    call dealloc_sparse(v2u)
    call dealloc_sparse(u2v)
    deallocate(u,v, u_v, v_u)

!==========================================================================================
!%==============================
!% set up the time integration:
!%==============================
    allocate(omega0(ncpts), ncycles(ncpts), stat=istat)
    call get_tidal_frequencies(cpts, omega0)

if (ncpts == 1) then
  tmaxd=(2*pi/omega0(1))/sind
  ncycles(1)=1
else
  ! time (in days) over which analysis is performed
!  print *,"P%ndays", P%ndays
  if (P%ndays == 0) then
    call calc_projection_period(cpts, omega0, tmaxd)
    P%ndays = tmaxd
  else
    tmaxd = P%ndays
  endif

  ncycles = floor(sind*tmaxd*omega0/(2*pi))
  omega0 = (2*pi/(sind*tmaxd))*ncycles
endif

dt = sind*tmaxd/(ppc*maxval(ncycles))

!==========================================================================================
!%==================================================
!% calculate the projection integrals (nh x ncpts):
!%==================================================

write(*, '("Projecting: ")', advance='no')
call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

    allocate(u_u(nu), u_v(nv), v_u(nu), v_v(nv), stat = istat)

    allocate(Iu(nu, ncpts), Iv(nv, ncpts), stat = istat)
    Iu=0.
    Iv=0.

do ct=1,ppc*maxval(ncycles)

  t=ct*dt

!  u_u = real( matmul(ustore_u,exp(-i*omega0*t)) )
!  v_u = real( matmul(vstore_u,exp(-i*omega0*t)) )
!  v_v = real( matmul(vstore_v,exp(-i*omega0*t)) )
!  u_v = real( matmul(ustore_v,exp(-i*omega0*t)) )

  	call mat_mat_mul_d(ustore_u,exp(-i*omega0*t), P%blas_num_threads, u_u)
	u_u = real(u_u)
  	call mat_mat_mul_d(vstore_u,exp(-i*omega0*t), P%blas_num_threads, v_u)
	v_u = real(v_u)
	call mat_mat_mul_d(ustore_v,exp(-i*omega0*t), P%blas_num_threads, u_v)
	u_v = real(u_v)
  	call mat_mat_mul_d(vstore_v,exp(-i*omega0*t), P%blas_num_threads, v_v)
	v_v = real(v_v)

  do m = 1,nu
  	tmp = sqrt(u_u(m)**2+v_u(m)**2)*u_u(m)
  	Iu(m,:) = Iu(m,:) + dt*(tmp * exp(i*omega0(:)*t) )
!    do n = 1,ncpts
!       Iu(m,n) = Iu(m,n) + dt*( tmp * exp(i*omega0(n)*t) )
!    enddo
  enddo
  do m = 1,nv
  	tmp = sqrt(u_v(m)**2+v_v(m)**2)*v_v(m)
  	Iv(m,:) = Iv(m,:) + dt*( tmp * exp(i*omega0(:)*t) )
!    do n = 1,ncpts
!       Iv(m,n) = Iv(m,n) + dt*( tmp * exp(i*omega0(n)*t) )
!    enddo
  enddo

  if (modulo(ct,maxval(ncycles)) == 0) then
      write(*, '("*")', advance='no')
  endif

enddo
deallocate(u_u, u_v, v_u, v_v)

Iu=(2/t)*Iu
Iv=(2/t)*Iv

  call CPU_Time(T2)
  call system_clock ( wall_t2, clock_rate, clock_max )
	if (P%messages >= 1 ) then
		write(*, '(a, a, a, a, a)') ' done (CPU: ', trim(ADJUSTL(conv_secs(T2-T1))), ', Wall: ',&
						trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) ))), ')'
	endif

!==========================================================================================
!%=======================
!% calculate the Q and D
!%=======================
call CPU_Time(T1)
call system_clock ( wall_t1, clock_rate, clock_max )

write(*, '("  Storing: ")', advance='no')

     call load_alloc_vector(H_u, dir_cols // 'H_u.dat')
     call load_alloc_vector(H_v, dir_cols // 'H_v.dat')
     call load_alloc_vector(ta_u, dir_cols // 'ta_u.dat')
     call load_alloc_vector(ta_v, dir_cols // 'ta_v.dat')

    allocate(Du(nu), Dv(nv), stat = istat)
do ccpt = 1, ncpts

  cpt=cpts(2*ccpt-1:2*ccpt)

	select case (cooruv)
		case (12) ! MERCATOR, TRANSPORT
			Du=cdg*Iu(:,ccpt)*(cosh(ta_u)/H_u)**2
			Dv=cdg*Iv(:,ccpt)*(cosh(ta_v)/H_v)**2
		case (22) ! LAT-LON, TRANSPORT
			Du=cdg*Iu(:,ccpt)/(H_u)**2
			Dv=cdg*Iv(:,ccpt)/(H_v)**2
	end select

     call save_vector(Du, dir_sols // cpt // '_Du' // '_m0_drag' // '.dat')
     call save_vector(Dv, dir_sols // cpt // '_Dv' // '_m0_drag' // '.dat')

     write(*, '("*")', advance='no')

enddo

deallocate(H_u, H_v, ta_u, ta_v, Du, Dv)

  call CPU_Time(T2)
  call system_clock ( wall_t2, clock_rate, clock_max )
	if (P%messages >= 1 ) then
		write(*, '(a, a, a, a, a)') ' done (CPU: ', trim(ADJUSTL(conv_secs(T2-T1))), ', Wall: ',&
						trim(ADJUSTL(conv_secs( real(wall_t2-wall_t1)/real(clock_rate) ))), ')'
	elseif (P%messages == 0 ) then
		write(*, '(" ")')
	endif

end subroutine calc_dhat

!**********************************************************************************************************
!**********************************************************************************************************

subroutine get_tidal_frequencies(cpts, omega0)

implicit none
!% for the string cpts='m2k1...', calculates an
!% array of frequencies omega=[omega1 omega2 ...]
     character(len=*)   :: cpts
     character(len=2)   :: cpt
     type (tide_params) :: pars
     real(wp)           :: omega0(*)
     integer            :: ccpt, ncpts!, istat



ncpts=len(cpts)/2

do ccpt = 1, ncpts

     cpt=cpts(2*ccpt-1:2*ccpt)
     pars = get_pars(cpt)
     omega0(ccpt) = pars%omega0

enddo

end subroutine get_tidal_frequencies

!**********************************************************************************************************
!**********************************************************************************************************
subroutine calc_projection_period(cpts, omega0, ndays)

implicit none
!% Given a list of components of the form 'm2k1...',
!% determine the shortest number of days over which the
!% number of complete tidal cycles of all constituents
!% will differ by at least one.
     character(len=*)   :: cpts
!     type (tide_params) :: pars_tmp
     real(wp)           :: omega0(:), ndays
     integer            :: ncpts, n_unique, cvgd, istat
     real(wp),allocatable:: period(:)
     integer,allocatable :: ncycles(:)


ncpts=len(cpts)/2
!%================
!% start the loop
!%================
ndays=0
cvgd=1
allocate(period(ncpts), ncycles(ncpts), stat=istat)
do while (cvgd > 0)
  ndays=ndays+1
  period=2*pi/omega0
  ncycles=floor(sind*ndays/period) !rounding here: complete cycles
  call unique(ncycles, n_unique)
  cvgd=(abs(size(ncycles)-n_unique))
!  print *, size(ncycles), n_unique
end do

end subroutine calc_projection_period
!**********************************************************************************************************

!**********************************************************************************************************
Subroutine unique(XDONT, NUNI)
!   unique = (Stable unique) Removes duplicates from an array,
!            leaving unique entries in the order of their first
!            appearance in the initial set.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (InOut)  :: XDONT
      Integer, Intent (Out) :: NUNI
! __________________________________________________________
!
      Integer, Dimension (Size(XDONT)) :: IWRKT
      Logical, Dimension (Size(XDONT)) :: IFMPTYT
      Integer :: ICRS
! __________________________________________________________
      Call UNIINV (XDONT, IWRKT)
      IFMPTYT = .True.
      NUNI = 0
      Do ICRS = 1, Size(XDONT)
         If (IFMPTYT(IWRKT(ICRS))) Then
            IFMPTYT(IWRKT(ICRS)) = .False.
            NUNI = NUNI + 1
            XDONT (NUNI) = XDONT (ICRS)
         End If
      End Do
      Return
!
End Subroutine unique

!**********************************************************************************************************
Subroutine uniinv (XDONT, IGOEST)
! __________________________________________________________
!   UNIINV = Merge-sort inverse ranking of an array, with removal of
!   duplicate entries.
!   The routine is similar to pure merge-sort ranking, but on
!   the last pass, it sets indices in IGOEST to the rank
!   of the value in the ordered set with duplicates removed.
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XDONT
      Integer, Dimension (:), Intent (Out) :: IGOEST
! __________________________________________________________
      Integer :: XTST, XDONA, XDONB
!
! __________________________________________________________
      Integer, Dimension (SIZE(IGOEST)) :: JWRKT, IRNGT
      Integer :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2, NUNI
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IGOEST))
!
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IGOEST (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) < XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo (NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 4) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (2*LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!  One steps in the C subset, that we create in the final rank array
!
!  Make a copy of the rank array for the iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
            XDONA = XDONT (JWRKT(IINDA))
            XDONB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XDONA > XDONB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XDONB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XDONA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
!   Last merge of A and B into C, with removal of duplicates.
!
      IINDA = 1
      IINDB = LMTNA + 1
      NUNI = 0
!
!  One steps in the C subset, that we create in the final rank array
!
      JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
      If (IINDB <= NVAL) Then
        XTST = NEARLESS (Min(XDONT(JWRKT(1)), XDONT(IRNGT(IINDB))))
      Else
        XTST = NEARLESS (XDONT(JWRKT(1)))
      Endif
      Do IWRK = 1, NVAL
!
!  We still have unprocessed values in both A and B
!
         If (IINDA <= LMTNA) Then
            If (IINDB <= NVAL) Then
               If (XDONT(JWRKT(IINDA)) > XDONT(IRNGT(IINDB))) Then
                  IRNG = IRNGT (IINDB)
                  IINDB = IINDB + 1
               Else
                  IRNG = JWRKT (IINDA)
                  IINDA = IINDA + 1
               End If
            Else
!
!  Only A still with unprocessed values
!
               IRNG = JWRKT (IINDA)
               IINDA = IINDA + 1
            End If
         Else
!
!  Only B still with unprocessed values
!
            IRNG = IRNGT (IWRK)
         End If
         If (XDONT(IRNG) > XTST) Then
            XTST = XDONT (IRNG)
            NUNI = NUNI + 1
         End If
         IGOEST (IRNG) = NUNI
!
      End Do
!
      Return
!
End Subroutine uniinv

Function nearless (XVAL) result (I_nl)
!  Nearest value less than given value
! __________________________________________________________
      Integer, Intent (In) :: XVAL
      Integer :: I_nl
! __________________________________________________________
      I_nl = XVAL - 1
      return
!
End Function nearless


end module baro_projection
