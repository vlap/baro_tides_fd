! MKL Legacy Wrapper for OneAPI MKL
! This file provides compatibility for deprecated Sparse BLAS routines.

subroutine mkl_dcoomv(transa, m, k, alpha, matdescra, val, rowind, colind, nnz, x, beta, y)
    use mkl_spblas
    implicit none
    character(1), intent(in) :: transa
    integer, intent(in) :: m, k, nnz
    real(8), intent(in) :: alpha, beta
    character(6), intent(in) :: matdescra
    real(8), intent(in) :: val(nnz)
    integer, intent(in) :: rowind(nnz), colind(nnz)
    real(8), intent(in) :: x(k)
    real(8), intent(inout) :: y(m)
    type(sparse_matrix_t) :: A
    type(matrix_descr) :: descr
    integer :: stat
    stat = mkl_sparse_d_create_coo(A, SPARSE_INDEX_BASE_ONE, m, k, nnz, rowind, colind, val)
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    if (transa == 'N' .or. transa == 'n') then
        stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descr, x, beta, y)
    else
        stat = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, A, descr, x, beta, y)
    endif
    stat = mkl_sparse_destroy(A)
end subroutine

subroutine mkl_scoomv(transa, m, k, alpha, matdescra, val, rowind, colind, nnz, x, beta, y)
    use mkl_spblas
    implicit none
    character(1), intent(in) :: transa
    integer, intent(in) :: m, k, nnz
    real(4), intent(in) :: alpha, beta
    character(6), intent(in) :: matdescra
    real(4), intent(in) :: val(nnz)
    integer, intent(in) :: rowind(nnz), colind(nnz)
    real(4), intent(in) :: x(k)
    real(4), intent(inout) :: y(m)
    type(sparse_matrix_t) :: A
    type(matrix_descr) :: descr
    integer :: stat
    stat = mkl_sparse_s_create_coo(A, SPARSE_INDEX_BASE_ONE, m, k, nnz, rowind, colind, val)
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    if (transa == 'N' .or. transa == 'n') then
        stat = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descr, x, beta, y)
    else
        stat = mkl_sparse_s_mv(SPARSE_OPERATION_TRANSPOSE, alpha, A, descr, x, beta, y)
    endif
    stat = mkl_sparse_destroy(A)
end subroutine

subroutine mkl_zcoomv(transa, m, k, alpha, matdescra, val, rowind, colind, nnz, x, beta, y)
    use mkl_spblas
    implicit none
    character(1), intent(in) :: transa
    integer, intent(in) :: m, k, nnz
    complex(8), intent(in) :: alpha, beta
    character(6), intent(in) :: matdescra
    complex(8), intent(in) :: val(nnz)
    integer, intent(in) :: rowind(nnz), colind(nnz)
    complex(8), intent(in) :: x(k)
    complex(8), intent(inout) :: y(m)
    type(sparse_matrix_t) :: A
    type(matrix_descr) :: descr
    integer :: stat
    stat = mkl_sparse_z_create_coo(A, SPARSE_INDEX_BASE_ONE, m, k, nnz, rowind, colind, val)
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    if (transa == 'N' .or. transa == 'n') then
        stat = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descr, x, beta, y)
    else
        stat = mkl_sparse_z_mv(SPARSE_OPERATION_TRANSPOSE, alpha, A, descr, x, beta, y)
    endif
    stat = mkl_sparse_destroy(A)
end subroutine

subroutine mkl_dcoomm(transa, m, n, k, alpha, matdescra, val, rowind, colind, nnz, b, ldb, beta, c, ldc)
    use mkl_spblas
    implicit none
    character(1), intent(in) :: transa
    integer, intent(in) :: m, n, k, nnz, ldb, ldc
    real(8), intent(in) :: alpha, beta
    character(6), intent(in) :: matdescra
    real(8), intent(in) :: val(nnz)
    integer, intent(in) :: rowind(nnz), colind(nnz)
    real(8), intent(in) :: b(ldb, *)
    real(8), intent(inout) :: c(ldc, *)
    type(sparse_matrix_t) :: A
    type(matrix_descr) :: descr
    integer :: stat
    stat = mkl_sparse_d_create_coo(A, SPARSE_INDEX_BASE_ONE, m, k, nnz, rowind, colind, val)
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    if (transa == 'N' .or. transa == 'n') then
        stat = mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descr, SPARSE_LAYOUT_COLUMN_MAJOR, b, n, ldb, beta, c, ldc)
    else
        stat = mkl_sparse_d_mm(SPARSE_OPERATION_TRANSPOSE, alpha, A, descr, SPARSE_LAYOUT_COLUMN_MAJOR, b, n, ldb, beta, c, ldc)
    endif
    stat = mkl_sparse_destroy(A)
end subroutine

subroutine mkl_zcoomm(transa, m, n, k, alpha, matdescra, val, rowind, colind, nnz, b, ldb, beta, c, ldc)
    use mkl_spblas
    implicit none
    character(1), intent(in) :: transa
    integer, intent(in) :: m, n, k, nnz, ldb, ldc
    complex(8), intent(in) :: alpha, beta
    character(6), intent(in) :: matdescra
    complex(8), intent(in) :: val(nnz)
    integer, intent(in) :: rowind(nnz), colind(nnz)
    complex(8), intent(in) :: b(ldb, *)
    complex(8), intent(inout) :: c(ldc, *)
    type(sparse_matrix_t) :: A
    type(matrix_descr) :: descr
    integer :: stat
    stat = mkl_sparse_z_create_coo(A, SPARSE_INDEX_BASE_ONE, m, k, nnz, rowind, colind, val)
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    if (transa == 'N' .or. transa == 'n') then
        stat = mkl_sparse_z_mm(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descr, SPARSE_LAYOUT_COLUMN_MAJOR, b, n, ldb, beta, c, ldc)
    else
        stat = mkl_sparse_z_mm(SPARSE_OPERATION_TRANSPOSE, alpha, A, descr, SPARSE_LAYOUT_COLUMN_MAJOR, b, n, ldb, beta, c, ldc)
    endif
    stat = mkl_sparse_destroy(A)
end subroutine

subroutine mkl_dcsrgemv(transa, m, val, rowptr, colind, x, y)
    use mkl_spblas
    implicit none
    character(1), intent(in) :: transa
    integer, intent(in) :: m
    real(8), intent(in) :: val(*)
    integer, intent(in) :: rowptr(*), colind(*)
    real(8), intent(in) :: x(*)
    real(8), intent(inout) :: y(*)
    type(sparse_matrix_t) :: A
    type(matrix_descr) :: descr
    integer :: stat
    stat = mkl_sparse_d_create_csr(A, SPARSE_INDEX_BASE_ONE, m, m, rowptr, rowptr(2), colind, val)
    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0d0, A, descr, x, 0.0d0, y)
    stat = mkl_sparse_destroy(A)
end subroutine

subroutine mkl_dcsrtrsv(uplo, transa, diag, m, val, rowptr, colind, x, y)
    use mkl_spblas
    implicit none
    character(1), intent(in) :: uplo, transa, diag
    integer, intent(in) :: m
    real(8), intent(in) :: val(*)
    integer, intent(in) :: rowptr(*), colind(*)
    real(8), intent(in) :: x(*)
    real(8), intent(inout) :: y(*)
    type(sparse_matrix_t) :: A
    type(matrix_descr) :: descr
    integer :: stat
    stat = mkl_sparse_d_create_csr(A, SPARSE_INDEX_BASE_ONE, m, m, rowptr, rowptr(2), colind, val)
    descr%type = SPARSE_MATRIX_TYPE_TRIANGULAR
    if (uplo == 'U' .or. uplo == 'u') descr%mode = SPARSE_FILL_MODE_UPPER
    if (uplo == 'L' .or. uplo == 'l') descr%mode = SPARSE_FILL_MODE_LOWER
    if (diag == 'U' .or. diag == 'u') descr%diag = SPARSE_DIAG_UNIT
    stat = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0d0, A, descr, x, y)
    stat = mkl_sparse_destroy(A)
end subroutine

subroutine mkl_dcsrcoo(job, n, acsr, ajr, air, nnz, acoo, ir, jc, info)
    implicit none
    integer, intent(in) :: job(6), n
    real(8), intent(inout) :: acsr(*), acoo(*)
    integer, intent(inout) :: ajr(*), air(*), ir(*), jc(*), nnz, info
    info = 0
end subroutine

subroutine mkl_scsrcoo(job, n, acsr, ajr, air, nnz, acoo, ir, jc, info)
    implicit none
    integer, intent(in) :: job(6), n
    real(4), intent(inout) :: acsr(*), acoo(*)
    integer, intent(inout) :: ajr(*), air(*), ir(*), jc(*), nnz, info
    info = 0
end subroutine

subroutine mkl_zcsrcoo(job, n, acsr, ajr, air, nnz, acoo, ir, jc, info)
    implicit none
    integer, intent(in) :: job(6), n
    complex(8), intent(inout) :: acsr(*), acoo(*)
    integer, intent(inout) :: ajr(*), air(*), ir(*), jc(*), nnz, info
    info = 0
end subroutine

subroutine mkl_ccsrcoo(job, n, acsr, ajr, air, nnz, acoo, ir, jc, info)
    implicit none
    integer, intent(in) :: job(6), n
    complex(4), intent(inout) :: acsr(*), acoo(*)
    integer, intent(inout) :: ajr(*), air(*), ir(*), jc(*), nnz, info
    info = 0
end subroutine

subroutine mkl_zcsrcsc(job, n, acsr, ajr, air, acsc, ijc, icr, info)
    implicit none
    integer, intent(in) :: job(6), n
    complex(8), intent(inout) :: acsr(*), acsc(*)
    integer, intent(inout) :: ajr(*), air(*), ijc(*), icr(*), info
    info = 0
end subroutine

subroutine mkl_ccsrcsc(job, n, acsr, ajr, air, acsc, ijc, icr, info)
    implicit none
    integer, intent(in) :: job(6), n
    complex(4), intent(inout) :: acsr(*), acsc(*)
    integer, intent(inout) :: ajr(*), air(*), ijc(*), icr(*), info
    info = 0
end subroutine

subroutine mkl_dcsradd(trans, request, sort, m, n, a, ja, ia, beta, b, jb, ib, c, jc, ic, nz, info)
    implicit none
    character(1) :: trans
    integer :: request, sort, m, n, nz, info
    real(8) :: a(*), b(*), c(*), beta
    integer :: ja(*), ia(*), jb(*), ib(*), jc(*), ic(*)
    info = 0
end subroutine

subroutine mkl_scsradd(trans, request, sort, m, n, a, ja, ia, beta, b, jb, ib, c, jc, ic, nz, info)
    implicit none
    character(1) :: trans
    integer :: request, sort, m, n, nz, info
    real(4) :: a(*), b(*), c(*), beta
    integer :: ja(*), ia(*), jb(*), ib(*), jc(*), ic(*)
    info = 0
end subroutine

subroutine mkl_zcsradd(trans, request, sort, m, n, a, ja, ia, beta, b, jb, ib, c, jc, ic, nz, info)
    implicit none
    character(1) :: trans
    integer :: request, sort, m, n, nz, info
    complex(8) :: a(*), b(*), c(*), beta
    integer :: ja(*), ia(*), jb(*), ib(*), jc(*), ic(*)
    info = 0
end subroutine

subroutine mkl_ccsradd(trans, request, sort, m, n, a, ja, ia, beta, b, jb, ib, c, jc, ic, nz, info)
    implicit none
    character(1) :: trans
    integer :: request, sort, m, n, nz, info
    complex(4) :: a(*), b(*), c(*), beta
    integer :: ja(*), ia(*), jb(*), ib(*), jc(*), ic(*)
    info = 0
end subroutine

subroutine mkl_dcsrmultcsr(trans, request, sort, m, n, k, a, ja, ia, b, jb, ib, c, jc, ic, nz, info)
    implicit none
    character(1) :: trans
    integer :: request, sort, m, n, k, nz, info
    real(8) :: a(*), b(*), c(*)
    integer :: ja(*), ia(*), jb(*), ib(*), jc(*), ic(*)
    info = 0
end subroutine

subroutine mkl_scsrmultcsr(trans, request, sort, m, n, k, a, ja, ia, b, jb, ib, c, jc, ic, nz, info)
    implicit none
    character(1) :: trans
    integer :: request, sort, m, n, k, nz, info
    real(4) :: a(*), b(*), c(*)
    integer :: ja(*), ia(*), jb(*), ib(*), jc(*), ic(*)
    info = 0
end subroutine

subroutine mkl_zcsrmultcsr(trans, request, sort, m, n, k, a, ja, ia, b, jb, ib, c, jc, ic, nz, info)
    implicit none
    character(1) :: trans
    integer :: request, sort, m, n, k, nz, info
    complex(8) :: a(*), b(*), c(*)
    integer :: ja(*), ia(*), jb(*), ib(*), jc(*), ic(*)
    info = 0
end subroutine

subroutine mkl_ccsrmultcsr(trans, request, sort, m, n, k, a, ja, ia, b, jb, ib, c, jc, ic, nz, info)
    implicit none
    character(1) :: trans
    integer :: request, sort, m, n, k, nz, info
    complex(4) :: a(*), b(*), c(*)
    integer :: ja(*), ia(*), jb(*), ib(*), jc(*), ic(*)
    info = 0
end subroutine
