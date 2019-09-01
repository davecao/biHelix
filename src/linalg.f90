! ******************************************************************************
!
! file: linalg.f90
!
!
! author: Cao Wei
! Timestamp: Sun Jul  7 19:40:49 2019
!
! Copyright (C) 2019 Cao Wei. All rights reserved.
!
!
! The following statement of license applies *only* to this header file,
! and *not* to the other files distributed with FFTW or derived therefrom:
!
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
! OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
! GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! ******************************************************************************

! gfortran -cpp -DQUAD_PRECISION -g -c linalg.f90
!   option -cpp: used to process preprocessing option

module linalg

  use precision
  use utilities
  use ieee_arithmetic

  real(DP), public, parameter :: TOL_SVD = 1e-13_dp

  private :: lapack_dgesvd_interface
  private :: L2norm, L1norm, Lmax

  public:: svd
  interface svd
     module procedure lapack_dgesvd_interface
  end interface svd

  public :: outer, cross, vecDegree
  public :: matinv, M33INV, inverse, inv
  public :: prj

  !public:: norm
  !interface norm
  !  module procedure L2norm ! interface for a module;procedure is implicit
  !  module procedure L1norm
  !  module procedure Lmax
  !end interface norm

contains
  !*****************************************************************************
  ! Description:
  !   clip values of a vector between a given range
  !
  ! Arguments:
  !    a(array of rank1): the first vector
  !    lowest : the lowest value of the given range
  !    highest: the highest value of the given range
  !
  ! Returns:
  !     clipped vector
  ! status:
  !     not implemented yet! Is it necessary?
  !*****************************************************************************
  !*****************************************************************************
  ! Description:
  !   Project the vector 'u' to the vector 'v'
  !                      dot(U, V)       V
  !       Prj(u to v) = ----------- * --------
  !                       norm(V)      norm(V)
  ! Arguments:
  !    u(array of rank1): the first vector
  !    v(array of rank1): the second vector
  !
  ! Returns:
  !    c(array of rank1): a component on vector b
  ! status:
  !
  !*****************************************************************************
  subroutine prj(u, v, p)
    real(DP), dimension(3), intent(in) :: u, v
    real(DP), dimension(3), intent(out) :: p
    real(DP), dimension(3) :: norm_v
    real(DP) :: magnitude_v, dot_uv

    magnitude_v = norm2(v)
    norm_v = v / magnitude_v
    dot_uv = dot_product(u, v)
    p = dot_uv / magnitude_v * norm_v

  end subroutine prj
  !*****************************************************************************
  !  matinv  -  Compute the inverse of a 3x3 matrix.
  !   from helanal.f
  !    Note: Could not process ill-conditioned matrix
  ! Arguments:
  !  h       = input 3x3 matrix to be inverted
  ! Return:
  !  s       = output 3x3 inverse of matrix A
  !
  !*****************************************************************************
  subroutine matinv(h, s)
    real(DP), dimension(3, 3), intent(in):: h
    real(DP), dimension(3, 3), intent(out):: s
    real(DP) :: c11,c12,c13,c21,c22,c23,c31,c32,c33
    real(DP) :: deth1,deth2,deth3,deth
    integer:: j, k
    do j = 1, 3
      do k=1, 3
        s(j,k) = 0.0000000
      end do
    end do
    c11=h(2,2)*h(3,3)-h(3,2)*h(2,3)
    c12=h(2,1)*h(3,3)-h(3,1)*h(2,3)
    c13=h(2,1)*h(3,2)-h(3,1)*h(2,2)
    c21=h(1,2)*h(3,3)-h(1,3)*h(3,2)
    c22=h(1,1)*h(3,3)-h(1,3)*h(3,1)
    c23=h(1,1)*h(3,2)-h(3,1)*h(1,2)
    c31=h(1,2)*h(2,3)-h(1,3)*h(2,2)
    c32=h(1,1)*h(2,3)-h(1,3)*h(2,1)
    c33=h(1,1)*h(2,2)-h(1,2)*h(2,1) 
    deth1=h(1,1)*c11
    deth2=-h(1,2)*c12
    deth3=h(1,3)*c13
    deth=deth1+deth2+deth3
    write(*, *) 'det=', deth
    if (deth.ne.0.0) then
      s(1,1)=c11/deth+s(1,1)
      s(1,2)=-c21/deth+s(1,2)
      s(1,3)=c31/deth+s(1,3)
      s(2,1)=-c12/deth+s(2,1)
      s(2,2)=c22/deth+s(2,2)
      s(2,3)=-c32/deth+s(2,3)
      s(3,1)=c13/deth+s(3,1)
      s(3,2)=-c23/deth+s(3,2)
      s(3,3)=c33/deth+s(3,3)
    else
      write(*,*) 'Matrix is singular'
    end if
    return
  end subroutine matinv
  !*****************************************************************************
  !  M33INV  -  Compute the inverse of a 3x3 matrix.
  !
  !  A       = input 3x3 matrix to be inverted
  !  AINV    = output 3x3 inverse of matrix A
  !  OK_FLAG = (output) .TRUE. if the input matrix could be inverted,
  !            and .FALSE. if the input matrix is singular.
  !*****************************************************************************
  subroutine M33INV (A, Ainv, OK_FLAG)
    implicit none

    real(DP), dimension(3, 3), intent(in):: A
    real(DP), dimension(3, 3), intent(out):: Ainv
    logical, intent(out):: OK_FLAG

    real(DP), parameter :: eps = 3.0d-10
    real(DP) :: det
    real(DP), dimension(3, 3) :: cofactor
    det =   A(1,1)*A(2,2)*A(3,3)  &
          - A(1,1)*A(2,3)*A(3,2)  &
          - A(1,2)*A(2,1)*A(3,3)  &
          + A(1,2)*A(2,3)*A(3,1)  &
          + A(1,3)*A(2,1)*A(3,2)  &
          - A(1,3)*A(2,2)*A(3,1)
    write(*, *) 'det=', det
    if (abs(det) .le. eps) then
      Ainv = 0.0d0
      OK_FLAG = .FALSE.
      return
    end if

    cofactor(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
    cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    cofactor(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
    cofactor(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
    cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
    cofactor(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
    cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
    cofactor(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

    Ainv = transpose(cofactor) / det

    OK_FLAG = .TRUE.

    return
  end subroutine M33INV

  !=============================================================================
  ! Inverse matrix
  ! Method: Based on Doolittle LU factorization for Ax=b
  ! Alex G. December 2009
  !-----------------------------------------------------------------------------
  ! input ...
  ! a(n,n) - array of coefficients for matrix A
  ! n      - dimension
  ! output ...
  ! c(n,n) - inverse matrix of A
  ! comments ...
  ! the original matrix a(n,n) will be destroyed
  ! during the calculation
  !=============================================================================
  subroutine inverse(a,c,n)
    implicit none
    integer n
    double precision a(n,n), c(n,n)
    double precision L(n,n), U(n,n), b(n), d(n), x(n)
    double precision coeff
    integer i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    ! step 1: forward elimination
    do k=1, n-1
       do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          do j=k+1,n
             a(i,j) = a(i,j)-coeff*a(k,j)
          end do
       end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
      L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
      do i=1,j
        U(i,j) = a(i,j)
      end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
      b(k)=1.0
      d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
    ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
    ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
        c(i,k) = x(i)
      end do
      b(k)=0.0
    end do
  end subroutine inverse

  !*****************************************************************************
  ! Description:
  !   Find inverse of a matrix using dgetrf and dgetri of LAPACK
  !
  ! Standard format:
  !    f90 and later
  !
  ! Status:
  !    OK
  ! Note: To ill-conditioned matrix, severe problem with dgetrf and dgetri
  !
  ! Arguments:
  !   A (array): a vector of MxN array
  !
  ! Returns:
  !   Ainv(array) : the inverse of A
  !*****************************************************************************
  function inv(A) result(Ainv)
    real(DP), dimension(:,:), intent(in) :: A
    real(DP), dimension(size(A,1),size(A,2)) :: Ainv

    real(DP), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A

    n = size(A,1)
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)
    if (info /= 0) then
      stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)
    !write(*, *) 'Check inverse'
    !call mat2d_print(A)
    !M = matmul(Ainv,A)
    if (info /= 0) then
      stop 'Matrix inversion failed!'
    end if
  end function inv

  !*****************************************************************************
  ! Description:
  !   Calculate the outer product of two vectors, a and b
  !   (Not used now)
  !
  ! Standard format:
  !    f90 and later
  !
  ! Status:
  !    OK
  ! Arguments:
  !   a (array): a vector of N elements
  !   b (array): a vector of N elements
  !
  ! Returns:
  !    N x N outer product
  !*****************************************************************************
  function outer(a, b)
    integer, parameter :: SP = kind(1.0)
    integer, parameter :: DP = kind(1.0d0)
    ! description of input and output
    real(DP), dimension(:), intent(in) :: a, b
    real(DP), dimension(size(a), size(b)) :: outer
    !f2py real(DP) intent(hide), depend(a, b):: outer

    outer = spread(a, dim=2, ncopies=size(b)) * &
         spread(b, dim=1, ncopies=size(a))
  end function outer


  !*****************************************************************************
  ! Description:
  !    Calculate cross product of two vectors, a and b
  !
  ! Standard format:
  !    f90 and later
  !
  ! Status:
  !    OK
  !
  ! Arguments:
  !   a (array): a vector of 3 elements
  !   b (array): a vector of 3 elements
  !
  ! Returns:
  !   a unit vector of 3 elements, which is perpendicular
  !   to the ab-plane
  !*****************************************************************************
  function cross(a, b)
    integer, parameter :: SP = kind(1.0)
    integer, parameter :: DP = kind(1.0d0)
    real(DP), dimension(3), intent(in) :: a, b
    real(DP), dimension(3) :: cross
    real(DP) :: length

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

    length = sqrt(sum(cross * cross))
    ! normalization
    cross = cross * 1.0 / length
  end function cross

  !*****************************************************************************
  ! Description:
  !   Calculate angle between any two vectors
  !   The two vectors will be normalized by L2 form
  ! Arguments:
  !    a(array of rank1): the first vector
  !    b(array of rank1): the second vector
  !    inRad (boolean): default is True, the returned angle is in radian
  !                    if False, the returned angle is in degree
  ! Returns:
  !     angle (real(DP)): the acuted angle between A and B
  !*****************************************************************************

  function vecDegree(a, b, radians) result(angle)
    real(DP), intent(in):: a(:)
    real(DP), intent(in):: b(:)
    logical, optional, intent(in) :: radians
    real(DP), dimension(size(a)) :: v1
    real(DP), dimension(size(b)) :: v2
    logical :: rrad = .TRUE.
    real(DP) :: angle_rad, angle
    integer :: v1_size, v2_size
    ! initialize variables
    angle = -10.0
    if (present(radians)) then
      rrad = radians
    end if
    v1_size = size(a)
    v2_size = size(b)
    if (v1_size.ne.v2_size) then
      stop 'The size of two vectors should be same.'
    end if
    v1 = a / norm2(a)
    v2 = b / norm2(b)
    angle_rad = dot_product(v1, v2)
    ! *-- clip [-1.0, 1.0]
    if (angle_rad <= -1.0d0) then 
      angle_rad = -1.0d0
    else if(angle_rad >= 1.0d0) then
      angle_rad = 1.0d0
    else
    end if
    angle = acos(angle_rad)  ! in radians
    if (.not.rrad) then 
      angle = angle * 180.0/ pi
      if (angle > 90) then
        angle = 180.0 - angle
      end if
    end if
  end function vecDegree

  !*****************************************************************************
  ! Description:
  !    normalize a vector with the maximum
  !
  ! Arguments:
  !    A(array of rank1): a d-dimensional vector
  !
  ! Returns:
  !   A normalized vector
  !*****************************************************************************
  function Lmax(A) result(res)
    real(DP), intent(in):: A(:)
    real(DP), dimension(size(A)) :: res
    real(DP) :: n_val
    n_val = maxval(A)
    res = A / n_val
  end function Lmax

  !*****************************************************************************
  ! Description:
  !    L1 normalize a vector
  !
  ! Arguments:
  !    A(array of rank1): a d-dimensional vector
  !
  ! Returns:
  !   A L1 normalized vector
  !*****************************************************************************
  function L1norm(A) result(res)
    real(DP), intent(in):: A(:)
    real(DP), dimension(size(A)) :: res
    real(DP) :: n_val
    if (ieee_support_inf(n_val)) then
      n_val = ieee_value(n_val,  ieee_negative_inf)
    end if
    n_val = sum(abs(A))
    res = A / n_val
  end function L1norm

  !*****************************************************************************
  ! Description:
  !   L2 form, normalize the vector/matrix
  !
  ! Arguments:
  !     A (array of rank 1)
  ! Returns
  !     res(arrray of rank 1): L2 normalized of A
  !*****************************************************************************
  function L2norm(A) result(res)
    real(DP), intent(in):: A(:)
    real(DP), dimension(size(A)) :: res
    real(DP) :: n_val
    ! n_val = sqrt(dot_product(A, A))
    n_val = norm2(A)
    res = A / n_val
  end function L2norm

  !*****************************************************************************
  ! Description:
  !   Compute the matrix singular value decomposition with DGESVD
  !   of Lapack library
  !
  !      U, S, V = svd(A)
  ! Arguments:
  !   A (array): M x N . M>>N  
  !   M (integer): M of A, the number of rows of A
  !   N (integer): N of A, the number of columns of A
  !
  ! Returns:
  !   U (array): M x M
  !   S (array): N x N, diagonal
  !   V (array): N x N
  !*****************************************************************************
  subroutine lapack_dgesvd_interface(A, M, N, U, S, VT)

    real(DP), dimension(M, N), intent(in):: A
    integer, intent(in):: M, N
    real(DP), dimension(M, M), intent(out):: U ! dummy
    real(DP), dimension(N), intent(out):: S
    real(DP), dimension(N, N), intent(out):: VT
    !real(DP), dimension(N, N):: VT
    real(DP), dimension(:), allocatable:: work
    integer :: lda, ldu, lwork, lwkopt, ldvt, info
    integer :: AllocateStatus, DeAllocateStatus
    character JOBU, JOBVT
    real(DP) :: eps, serrbd
    !real(DP), dimension(N) :: rcondu, rcondv, uerrbd, verrbd

    ! *-- external subroutine --*
    !external ddisna
    external dgesvd

    ! rcondu(nmax), rcondv(nmax)
    ! ----------------------------------------------------------------
    !  Compute the singular values and left and right singular vectors
    !        of A (A = U*S*(V**T), m.ge.n)
    ! ----------------------------------------------------------------
    JOBU = 'A'
    JOBVT = 'A'
    lda = M
    ldu = M
    ldvt = N

    eps = dlamch('Eps')
    serrbd = eps * S(1)

    lwork = max(1, 3*min(M, N) + max(M,N), 5*min(M, N))

    allocate(work(lwork), stat=AllocateStatus)
    if ( AllocateStatus /= 0 ) stop 'Failed to allocate memory for work' 

    call dgesvd(JOBU, JOBVT, M, N, A, lda, S, U, ldu, VT, ldvt, &
         work, lwork, info)
    lwkopt = work(1)

!    if (info == 0) then
!       ! *-- print singular values, implicit loop --*
!       print *, 'Singular values'
!       print *, (S(j), j=1, N)
!       !
!       !  Call DDISNA (F08FLF) to estimate reciprocal condition
!       !  numbers for the singular vectors
!       !
!       call ddisna('Left', M, N, S, rcondu, info)
!       call ddisna('Right', N, N, S, rcondv, info)
!       !
!       !  Compute the error estimates for the singular vectors
!       !
!       do i = 1, n
!          uerrbd(i) = serrbd/rcondu(i)
!          verrbd(i) = serrbd/rcondv(i)
!       end do
!       !
!       ! Print the approximate error bounds for the singular values
!       ! and vectors
!       !
!       write (nout, *) ''
!       write (nout, *) 'Error estimate for the singular values'
!       write (nout, *) serrbd
!       write (nout, *) ''
!       write (nout, *) 'Error estimates for the left singular vectors'
!       write (nout, *)(uerrbd(i), i=1, n)
!       write (nout, *) ''
!       write (nout, *) 'Error estimates for the right singular vectors'
!       write (nout, *)(verrbd(i), i=1, n)
!    else
!       write (nout, *) 'Failure in DGESVD. INFO =', info
!    end if
!    ! *-- info --*
!    if (lwork < lwkopt) then
!       write (nout, *) 'Optimum workspace required = ', lwkopt, &
!            'provided         = ', lwork
!    end if
    !V = transpose(VT)
    !V = VT
    deallocate(work, stat=DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop 'Failed to release memory for work'
  end subroutine lapack_dgesvd_interface


end module linalg
