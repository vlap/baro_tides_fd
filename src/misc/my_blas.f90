module my_blas

     use precisions, only: sp, dp, wp, cwp
	 use dispmodule

	integer, parameter :: auto_num_threads = 0 ! choose num_threads automatically
!==========================================================================================
!**************************** INTERFACES **************************************************
!==========================================================================================
!---------------------------- DOT PRODUCT ---------------------------------------
     interface vec_vec_dot
          module procedure vec_vec_dotproduct_sp
          module procedure vec_vec_dotproduct_dp
     end interface  
!---------------------------- INITIATE WITH VALS -------------------------------------
     interface mat_mat_mul_d
          module procedure dense_vec_mult
          module procedure dense_vec_mult_cmplx_cmplx
     end interface  
!==========================================================================================

contains
!==========================================================================================
real(sp) function vec_vec_dotproduct_sp(x, y, n_threads)
! c := mat*vec
!or
! vec := mat*vec
	implicit none
	include 'mkl_blas.fi'


	 integer :: nx, ny, n
     real(sp):: x(:), y(:)
     integer, optional :: n_threads
     integer :: num_threads, mkl_get_max_threads

     integer   :: incx, incy

	num_threads = auto_num_threads
	nx = size(x)
	ny = size(y)
	if (nx/= ny) then
		print *, "vec_vec_dotproduct: dimensions must agree nx \= ny:", nx, ny
		stop
	end if
	n = nx

	if (present(n_threads)) then
		num_threads = n_threads
	endif

if (num_threads < 0) then
	! built in Fortran function
	vec_vec_dotproduct_sp = dot_product(x,y)
else
	!mkl
	!mkl
	if ( (num_threads >0).and.( num_threads /= mkl_get_max_threads()) ) then
		call mkl_set_num_threads(num_threads) ! # threads
	endif


	incx = 1 ! increment for the elements of x.
	incy = 1 ! increment for the elements of y.

	    vec_vec_dotproduct_sp = sdot(n, x, incx, y, incy)

endif

end function vec_vec_dotproduct_sp
!==========================================================================================
real(dp) function vec_vec_dotproduct_dp(x, y, n_threads)
! c := mat*vec
!or
! vec := mat*vec
	implicit none
	include 'mkl_blas.fi'


	 integer :: nx, ny, n
     real(dp):: x(:), y(:)
     integer, optional :: n_threads
     integer :: num_threads, mkl_get_max_threads

     integer   :: incx, incy

	num_threads = auto_num_threads
	nx = size(x)
	ny = size(y)
	if (nx/= ny) then
		print *, "vec_vec_dotproduct: dimensions must agree nx \= ny:", nx, ny
		stop
	end if
	n = nx

	if (present(n_threads)) then
		num_threads = n_threads
	endif

if (num_threads < 0) then
	! built in Fortran function
	vec_vec_dotproduct_dp = dot_product(x,y)
else
	!mkl
	!mkl
	if ( (num_threads >0).and.( num_threads /= mkl_get_max_threads()) ) then
		call mkl_set_num_threads(num_threads) ! # threads
	endif


	incx = 1 ! increment for the elements of x.
	incy = 1 ! increment for the elements of y.

		vec_vec_dotproduct_dp = ddot(n, x, incx, y, incy)

endif

end function vec_vec_dotproduct_dp
!==========================================================================================
!==========================================================================================
subroutine dense_vec_mult(mat, vec, n_threads, c)
! c := mat*vec
!or
! vec := mat*vec
implicit none

	 integer :: m, n, nvec, nc
     real(wp)                         :: mat(:,:), vec(:)
     real(wp), allocatable            :: tmp(:)
     real(wp), optional               :: c(:)
     integer :: n_threads, mkl_get_max_threads

     integer   :: istatus, lda, incx, incy
     character :: transa
     real(wp)  :: alpha, beta

	m = size(mat,1)
	n = size(mat,2)
	nvec = size(vec)

	if (n /= nvec) then
		print *, "dense_vec_mult: dimensions must agree n \= nvec:", n, nvec
		stop
	elseif (present(c)) then
	    nc = size(c)
	    if (m /= nc)  then
	    	print *, "dense_vec_mult: dimensions must agree m \= nc:", m, nc
		stop
		endif
	end if

if (n_threads < 0) then
	! built in Fortran function
    if (present(c)) then
          c = matmul(mat, vec)
    else
         vec = matmul(mat, vec)
    endif
else
	!mkl
	if ( (n_threads >0).and.( n_threads /= mkl_get_max_threads()) ) then
		call mkl_set_num_threads(n_threads) ! # threads
	endif

	transa = 'N'
	alpha = 1.
	beta = 0.
	lda = m
	incx = 1 ! increment for the elements of x.
	incy = 1 ! increment for the elements of y.

	if (wp==dp) then
	    if (present(c)) then
	          call dgemv(transa,m,n,alpha,mat,lda,vec,incx,beta,c,incy)
	    else
	          allocate(tmp(m),  stat = istatus)
	          call dgemv(transa,m,n,alpha,mat,lda,vec,incx,beta,tmp,incy)
	          ! copy result back to vec
	          vec = tmp
	    endif
	elseif (wp == sp) then
	    if (present(c)) then
	          call sgemv(transa,m,n,alpha,mat,lda,vec,incx,beta,c,incy)
	    else
	          allocate(tmp(m),  stat = istatus)
	          call sgemv(transa,m,n,alpha,mat,lda,vec,incx,beta,tmp,incy)
	          ! copy result back to vec
	          vec = tmp
	    endif
	endif
endif

end subroutine dense_vec_mult
!==========================================================================================
subroutine dense_vec_mult_cmplx_cmplx(mat, vec, n_threads, c)
! c := mat*vec
!or
! vec := mat*vec
implicit none

	 integer :: m, n, nvec, nc
     complex(cwp)                         :: mat(:,:), vec(:)
     complex(cwp), allocatable            :: tmp(:)
     complex(cwp), optional               :: c(:)
     integer :: n_threads, mkl_get_max_threads

     integer   :: istatus, lda, incx, incy
     character :: transa
     complex(cwp)  :: alpha, beta

	m = size(mat,1)
	n = size(mat,2)
	nvec = size(vec)

	if (n /= nvec) then
		print *, "dense_vec_mult: dimensions must agree n \= nvec:", n, nvec
		stop
	elseif (present(c)) then
	    nc = size(c)
	    if (m /= nc)  then
	    	print *, "dense_vec_mult: dimensions must agree m \= nc:", m, nc
		stop
		endif
	end if

if (n_threads < 0) then
	! built in Fortran function
    if (present(c)) then
          c = matmul(mat, vec)
    else
         vec = matmul(mat, vec)
    endif
else
	!mkl
	if ( (n_threads >0).and.( n_threads /= mkl_get_max_threads()) ) then
		call mkl_set_num_threads(n_threads) ! # threads
	endif

	transa = 'N'
	alpha = 1.
	beta = 0.
	lda = m
	incx = 1 ! increment for the elements of x.
	incy = 1 ! increment for the elements of y.

	if (cwp==dp) then
	    if (present(c)) then
	          call zgemv(transa,m,n,alpha,mat,lda,vec,incx,beta,c,incy)
	    else
	          allocate(tmp(m),  stat = istatus)
	          call zgemv(transa,m,n,alpha,mat,lda,vec,incx,beta,tmp,incy)
	          ! copy result back to vec
	          vec = tmp
	    endif
	elseif (cwp == sp) then
	    if (present(c)) then
	          call cgemv(transa,m,n,alpha,mat,lda,vec,incx,beta,c,incy)
	    else
	          allocate(tmp(m),  stat = istatus)
	          call cgemv(transa,m,n,alpha,mat,lda,vec,incx,beta,tmp,incy)
	          ! copy result back to vec
	          vec = tmp
	    endif
	endif
endif

end subroutine dense_vec_mult_cmplx_cmplx
!==========================================================================================
end module my_blas
