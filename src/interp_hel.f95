!*******************************************************************************
! Description:
!
! Catmull-Rom interpolation
! 
!
! Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!
!
!  Created:
!    11 Janurary 2016
!
!  Author:
!
!    Wei Cao
!*******************************************************************************

module interp_hel
  use precision
  use linalg
  ! ****** 3d basis ******
  !spline (u [0,1])
  !
  !                            |2, -2, 1,   1| | x1,  y1,  z1|
  !(x, y, z) =|u^3, u^2, u, 1| |-3, 3, -2, -1| | x2,  y2,  z2|
  !  point       parameter     |0, 0,    1, 0| |x1', y1', z1'|
  !   on          vector       |1, 0,    0, 0| |x2', y2', z2'|
  !  the spline                    basis     control matrix
  ! *--- b splines basis ---*
  real(DP) :: bsplinesBasis(4, 4) = reshape( &
    (/-1.0,  3.0, -3.0, 1.0, &
       3.0, -6.0,  3.0, 0.0, -3.0,  0.0,  3.0, 0.0, 1.0,  4.0,  1.0, 0.0/), (/4&
       &, 4/)) 
  ! *--- cubic hermite basis
  real(DP) :: cubicHermiteBasis(4, 4) = reshape((/ 2.0, -2.0, -1.0,  1.0, -3.0,&
       &  3.0, -2.0, -1.0, 0.0,  0.0,  1.0,  0.0, 1.0,  0.0,  0.0,  0.0/),(/4,&
       & 4/)) 
  ! *--- bezier basis ---*
  real(DP) :: bezierBasis(4, 4) = reshape((/-1.0,  3.0, -3.0,  1.0, 3.0, -6.0, &
       & 3.0,  0.0, -3.0,  3.0,  0.0,  0.0, 1.0,  0.0,  0.0,  0.0/),(/4, 4/))
  ! *--- catmull rom basis ---*
  ! s : is so-called tension parameter
  ! For centripetal Catmull-Rom spline, the value of s(\alpha) is 0.5
  ! When s(\alpha) = 0, the resulting curve is the standard Catmull-Rom spline 
  !     or uniform Catmull-Rom spline.
  ! When s(\alpha) = 1, the product is a chordal Catmull-Rom spline.
  !| 2, -2, 1,  1| | 0, 1, 0, 0|   | -s,  2-s,   s-2,  s|
  !|-3,  3,-2, -1| | 0, 0, 1, 0| = |2*s,  s-3, 3-2*s, -s|
  !| 0,  0, 1,  0| |-s, 0, s, 0|   | -s,    0,     s,  0|
  !| 1,  0, 0,  0| | 0,-s, 0, s|   |  0,    1,     0,  0|
  ! ---------------------------------------------------------------------------
  real(DP) :: cRomBasis(4, 4) = reshape((/ 2.0, -2.0,  1.0,  1.0, -3.0,  3.0, &
       &-2.0, -1.0, 0.0,  0.0,  1.0,  0.0, 1.0,  0.0,  0.0,  0.0/),(/4, 4/))
  public :: cRomBasis

contains
  function catmullRomBasis(s) result(cr)
    real(DP), intent(in):: s
    real(DP) :: cr(4, 4)

    cr(1,1) = -s
    cr(1,2) = 2.0 - s
    cr(1,3) = s - 2.0
    cr(1,4) = s

    cr(2,1) = 2.0*s
    cr(2,2) = s - 3.0
    cr(2,3) = 3.0 - 2.0*s
    cr(2,4) = -s

    cr(3,1) = -s
    cr(3,2) = 0.0
    cr(3,3) = s
    cr(3,4) = 0.0

    cr(4,1) = 0.0
    cr(4,2) = 1.0
    cr(4,3) = 0.0
    cr(4,4) = 0.0
    !real(DP) :: catmullRomBasis(4, 4) = reshape( &
    !(/   -s,  2.0-s,     s-2.0,  s, &
    !  2.0*s,  s-3.0, 3.0-2.0*s, -s, &
    !     -s,    0.0,         s,0.0, &
    !    0.0,    1.0,       0.0,0.0/), (/4, 4/))
    !real(DP) :: basis(4, 4) = reshape( &
    !(/   0.0, 1.0, 0.0, 0.0, &
    !     0.0, 0.0, 1.0, 0.0d0, &
    !      -s, 0.0,  -s, 0.0, &
    !     0.0,  -s, 0.0,   s/), (/4, 4/))
    !catmullRomBasis = matmul(cRomBasis, basis)

  end function catmullRomBasis
end module interp_hel
