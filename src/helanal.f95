!*******************************************************************************
! Description:
!
!
! Reference:
!  Kumar, S. and Bansal, M. (1996). Structural and Sequence Characteristics
! of Long Alpha Helices in Globular Proteins. Biophysical J.,71, 1574-1586.
!
!  Sugeta H, Miyazawa T. General method for calculating helical parameters of
! polymer chains from bond lengths, bond angles, and internal-rotation angles.
! Biopolymers. 1967;5: 673–679. 
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
!
! Compilation: (DEBUG)
!  gfortran:
!     gfortran -g -cpp -DQUAD_PRECISION -c helanal.f95
!     gfortran -g -DQUAD_PRECISION -m64 -fbounds-check -Wall -fbacktrace -Wtabs\
!            -finit-real=nan \
!            test_f95.f95 helanal.o -o test_f95
!  f2py: manually generate python extension or see setup.py 
!    f2py-2.7 helanal.f95 -m _helanal -h helanal.pyf --overwrite-signature
!    f2py-2.7 -c helanal.pyf helanal.f95
!*******************************************************************************
!****************************************************************************
!                                                                           *
! PROGRAM TO CHARACTERISE THE GEOMETRIES OF ALPHA HELICES IN PROTEINS.      *
!                                                                           *
! AUTHORS: Sandeep Kumar and Prof. Manju Bansal,                            *
!          MBU, Indian Institute of Science, Bangalore 560012, India        *
!                                                                           *
!  e-mail address: mb@mbu.iisc.ernet.in                                     *
!                                                                           *
! Method outlined in the following paper:                                   *
!                                                                           *
! Kumar, S. and Bansal, M. (1996). Structural and Sequence Characteristics  *
! of Long Alpha Helices in Globular Proteins. Biophysical J.,71, 1574-1586. *
!                                                                           *
! If any Bug is found, please report it to the authors.                     *
!                                                                           *
! SUMMARY OF THE ALGORITHM                                                  *
!                                                                           *
! Geometry of an alpha helix is characterised in terms of the angles between*
! local helix axes and the path traced by the local helix origins,          *
! calculated using the procedure of Sugeta and Miyazawa (1967), for every   *
! set of four contiguous C-alpha atoms, and sliding this window over the    *
! length of the helix in steps of one C-alpha atom.                         *
! Matrix M(I, J) contains the bending angles between local helix axes I & J *
! which are used to characterize the overall geometry of the  helix.        *
! The local helix origins trace out the path described by the helix in 3-D  *
! space. These origins are reoriented in X-Y plane and the reoriented       *
! points are used to fit a circle and a line by least squares method.       *
! Unit twist and unit height of the alpha helix are also calculated.        *
!                                                                           *
! A maximum of 5000 helices, each with 100 or less number of residues can   *
! be analysed at one time.                                                  *
!                                                                           *
! In order to analyse a larger number of helices, increase the value of 'i' *
! in the  following statement:                                              *
!                                                                           *
!      parameter (i=5000) in the main program.                              *
!                                                                           *
! In order to analyse helices with lengths greater than 100 residues,       *
! increase the dimensions of the appropriate variables in the dimension     *
! statements of the main program and subroutines.                           *
!                                                                           *
! This program has several options for input and is fully interactive.      *
!                                                                           *
! INPUT FILES :                                                             *
!                                                                           *
! Input files to this program can be different depending upon the options   *
! chosen. This program can directly read HELIX records in one or more PDB   *
! files. In order to analyse helices found in helix records of a PDB file,  *
! give the PDB file name as input to the program. In order to analyse the   *
! helices found in the HELIX records of more than one PDB  files, give the  *
! PDB file names sequentially when prompted by the program or input the     *
! name of the file containing the PDB files names in format (5x,a11) as the *
! input to the program.                                                     *
!                                                                           *
! HELIX records can be read from PDB files or can be read from a file       *
!      containing information about the helix start/end residues written in *
!      the same format as PDB HELIX records or in a different format, in    *
!      which case the format has to be keyed in when running the program.   *
! In the last case, files other than the PDB files or files containing      *
!      C-alpha coordinates in the PDB format, can also be given as input to *
!      the program, upon specifying their format.                           *
!                                                                           *
!                                                                           *
! All PDB files and other input files (if any) should be in the same        *
! directory.                                                                *
!                                                                           *
! OUTPUT FILES :                                                            *
!                                                                           *
! For each  run of HELANAL on file(s) containing alpha helices with length  *
! greater than or equal to 9 residues, the following output files are       *
! created.                                                                  *
!                                                                           *
! RUN.ANS contains the questions and their answers during a run of HELANAL. *
!                                                                           *
! HELINFO.OUT file created only when the HELIX records in the PDB files are *
! used for information on helix start/end residues, contains these helix    *
! records.                                                                  *
!                                                                           *
! HELCA.OUT contains Coordinates of the C-alpha atoms constituting the      *
! helices.                                                                  *
!                                                                           *
! AXES.OUT contains the local helix axes fitted to 4 consecutive C-alpha    *
! atoms along with a matrix  M(I, J) whose elements are the angles between  *
! local helix axes I and J.                                                 *
!                                                                           *
! ANGLE.OUT contains the angle between successive local helix axes. It also *
! lists mean bending angle and maximum bending angle for each helix.        *
!                                                                           *
! ORIGIN.OUT contains the local helix origins for the helix along with the  *
! statistics obtained by fitting least square plane, circle and line to the *
! local helix origins.                                                      *
!                                                                           *
! NH.OUT contains unit height and unit twist for every turn of the helix as *
! well as average unit height and average unit twist for the whole helix.   *
!                                                                           *
! ***.PRM contains summary of various parameters obtained by HELANAL.       *
!                                                                           *
! ***.TAB contains the parameters from ***.PRM file in a tabular form along *
! the overall geometry assignment.                                          *
!                                                                           *
!****************************************************************************
module helanal
  use precision
  use utilities
  use linalg
  ! *-- Global --*
  implicit none

  private
  ! *-- symbolic name: f2py could not recognize them if used in sub/func --*
  ! *-- f2py not recognize them --* 
!#ifdef QUAD_PRECISION
!  integer, private, parameter :: DP = 16 ! kind(1.0d0)
!#elif TEN_DIGIT_PRECISION
!  integer, private, parameter :: DP = selected_real_kind(10) !kind(1.0d0)
!#else
!  integer, private, parameter :: DP = 8 ! kind(1.0d0)
!#endif

  !integer, parameter :: SP = kind(1.0)
  !integer, parameter :: DP = kind(1.0d0)
  ! *-- public: default is private --*
  public :: fit
  !public :: local_helix
  !public :: fit_origins_lsq

contains

  !*****************************************************************************
  ! Description:
  !   For a given 1d array, compute average, standard deviation, mean absolute
  !   deviation
  !
  ! Standard format:
  !    f90 and later
  !
  ! Status:
  !    OK 
  !
  ! Arguments:
  !   p (array): a vector of 1xN array
  !
  ! Returns:
  !   ave (real) : average of p, sum(x)/npts
  !   var (real) : variance of p, E^2-E
  !   sdev (real) : 
  !   mean_abs_std(real)
  !*****************************************************************************
  subroutine stats(p, ave, var, sdev, mean_abs_std)
    real(DP), dimension(:), intent(in):: p
    real(DP), intent(OUT) :: ave, var, sdev, mean_abs_std
    real(DP) :: sum_
    integer :: n
    n = size(p)
    sum_ = sum(p)
    ave = sum(p) / n
    var = (dot_product(p, p) - sum_*sum_/n)/(n-1)
    if (var < 0.) then
      ! an error exist
      write(*, *) 'Negative variance: ', var
      sdev = -1
    else
      sdev = sqrt(var)
    endif
    mean_abs_std = sum(abs(p - ave))/n
  end subroutine stats


  !*****************************************************************************
  ! Description:
  !   For a given four continuous 3d points, find direction along the helix 
  !   axis, origins, twist angle, and height of the turn of a helix
  !  
  ! Standard format:
  !    f90 and later
  !
  ! Status:
  !    OK 
  !
  ! Arguments:
  !   p1 (array): a vector of 1x3 array
  !   p2 (array): a vector of 1x3 array
  !   p3 (array): a vector of 1x3 array
  !   p4 (array): a vector of 1x3 array
  !
  ! Returns:
  !   direct
  !   origin
  !   twist
  !   height
  !*****************************************************************************
  subroutine local_helix(p1, p2, p3, p4, direct, origin, twist, height, info)
    !integer, parameter :: SP = kind(1.0)
    !integer, parameter :: DP = kind(1.0d0)
    implicit none
    real(DP), dimension(:), intent(in):: p1, p2, p3, p4
    logical, optional :: info
    real(DP), dimension(:), intent(out):: direct
    real(DP), dimension(:, :), intent(out):: origin
    real(DP), intent(out):: twist, height

    ! local variables
    real(DP), dimension(3) :: v12, v23, v34, dv13, dv24
    real(DP) :: length1, length2, r, cos_theta
    ! intrinsic :: norm2
    real(DP) :: norm2
    logical :: quiet

10  format(1x,'twist:',1X,F8.3,1X,'radius:',F8.3,1X,'height:',F8.3)

    if (present(info)) then
      quiet = info
    else
      quiet = .FALSE.
    endif

    v12 = p2 - p1
    v23 = p3 - p2
    v34 = p4 - p3

    dv13 = v12 - v23
    dv24 = v23 - v34

    ! normal vector
    direct = cross(dv13, dv24)

    ! angle
    length1 = sqrt(dot_product(dv13, dv13))
    length2 = sqrt(dot_product(dv24, dv24))
    cos_theta = dot_product(dv13, dv24) / (length1 * length2)
    !write(*, *) "length1", length1, "length2", length2, "cos_theta", cos_theta
    ! twist: in degree
    twist = acos(cos_theta) * 180.0 / pi

    ! radius of local helix cylinder
    r = sqrt(length1 * length2) / (2.0 * (1 - cos_theta))

    ! height of local helix cylinder
    height = dot_product(v23, direct)

    ! unit vectors
    dv13 = dv13 / norm2(dv13, 1)
    dv24 = dv24 / norm2(dv24, 1)

    ! origins
    origin(1, :) = p2 - dv13 * r
    origin(2, :) = p3 - dv24 * r
    if (.not.(quiet)) then
      write(6, *) 'cross product:'
      write(6, '(*(F8.3,2X))') direct
      write(6, 10) twist, r, height
      write(6, *) 'origins:'
      call report(origin, nout=6)
    endif
    ! rotation vector
    ! rot_vectors(1, :) = dv13
    ! rot_vectors(2, :) = dv24

  end subroutine local_helix

  !*****************************************************************************
  ! Description:
  !   Fit least-squared plane to local helix origins (Not finished )
  !   (not used now)
  ! Standard format:
  !   f90 and later
  !
  ! Arguments:
  !   origins (array): a vector of Nx3 array
  !   reference_axis (array, optional): a vector of 1x3 array default is [0, 0, 1]
  !
  ! Returns:
  !   tilt_angle (dp): angle  between helix axis and reference axis
  !   axvec (array of 3): normal vector from svd
  !   radc (dp): 
  !   rmsdc (dp):
  !   rmsdl (dp):
  !   r2 (dp):
  !*****************************************************************************
  subroutine fit_circle_lsq(origins, tilt_angle, axvec, reference_axis)
    !                      radc, rmsdc, rmsdl, r2, reference_axis)
    !integer, parameter :: SP = kind(1.0)
    !integer, parameter :: DP = kind(1.0d0)
    real(DP), dimension(:,:), intent(in):: origins
    real(DP), dimension(3), optional, intent(in) :: reference_axis
    real(DP), intent(out) :: tilt_angle
    !real(DP), intent(out) :: radc, rmsdc, rmsdl, r2
    real(DP), dimension(:), intent(out) :: axvec

    ! *-- local variables --*
    integer :: nrows, ncols
    real(DP), dimension(size(origins, 2)):: centroid ! dim=2: get columns
    real(DP), dimension(size(origins, 1), size(origins, 2)):: A
    real(DP), dimension(size(origins, 2), size(origins, 2)):: U
    real(DP), dimension(size(origins, 2)):: S
    real(DP), dimension(size(origins, 2), size(origins, 2)):: VT, covMat
    real(DP), dimension(3) :: ref_axis = (/0.0, 0.0, 1.0/)
    real(DP) :: ax(3) = 0.0
    real(DP) :: agreement = 0.0

    ! *-- Optional argument --*
    if ( present(reference_axis) ) then
      ref_axis = reference_axis
    endif

    nrows = size(origins, 1)
    ncols = size(origins, 2)
    ! get arithmetic mean of each cols
    centroid = sum(origins, 1) / nrows

    ! centralize
    ! helanl.f did not do this step
    A = origins - spread(centroid, dim=1, ncopies=nrows)
    ! covariance matrix
    ! A: N x d
    ! A^T * A = (d x N) (N x d) = (d, d)
    covMat = matmul(transpose(A), A)
    nrows = size(covMat, 1)
    ncols = size(covMat, 2)
    ! * -------------------------------------------------
    ! Singular Value Decomposition:
    ! svd: covMat will be changed in svd
    ! U is column-wised singular vectors
    ! V is row-wised singular vectors
    ! Note:  pay attention to the ill-conditioned matrix
    ! *--------------------------------------------------
    call svd(covMat, nrows, ncols, U, S, VT)

    ! Point to the first point
    agreement = vecDegree(A(1,:), VT(1,:))
    ! test agreement in (-90, 90)
    if ((agreement < pi/2.0).and.(agreement > -pi/2.0)) then
      ax = VT(1, :)
    else
      ax = -1*VT(1, :)
    end if
    tilt_angle = vecDegree(ax, reference_axis, radians=.FALSE.)
    axvec = ax

    ! *-- SVD -> pseudo inv--*

    !bp = sum(origins, 2)
    !write(*, *) 'bp:'
    !write(*, '(*(F8.3, 2X))') bp
    !matp = matp + outer(bp, bp)
    !write(*, *) 'matp'
    !call mat2d_print(matp)

    !matp_inv = inv(matp)
    !call inverse(matp, matp_inv, 3) 
    !call M33INV (matp, matp_inv, OK_FLAG)
    !if (.not. OK_FLAG) then
    !  write(*, *) 'The input matrix is singular'
    !  stop 
    !end if
    !call matinv(matp, matp_inv)
    !write(*, *) 'matp inv:'
    !call mat2d_print(matp_inv)
    !write(*, *) 'matp after:'
    !call mat2d_print(matp)

    !ap = sum(matp_inv * spread(bp, dim=1, ncopies=nrows), 1)
    !write(*, *) 'ap:'
    !write(*, '(*(F8.3, 2X))') ap
    !rmp = origins * spread(ap, dim=1, ncopies=nrows)

  end subroutine fit_circle_lsq


  !*****************************************************************************
  ! Description:
  !   For a set of CA coordinates, find bending angles along the helix axis
  !   using a sliding window of 9 points.
  !
  ! Standard format:
  !   f90 and later
  !
  ! Status:
  !   partially completed. fit_origins_lsq will be completed in future.
  !
  ! Arguments:
  !   points (array): Nx3 , row-wised matrix
  !
  ! Returns:
  !   bending_angles (array):
  !   directions (array): helix axis
  !   origins (array): origins along helix axis
  !   radc (real(DP)):
  !   rmsdc (real(DP)):
  !   rmsdl (real(DP)):
  !   r2 (real(DP)) : estimated radius
  !*****************************************************************************
  subroutine fit(points, bending_angles, directions, origins, &
                 tilt, reference_axis, info)

    !             tilt, radc, rmsdc, rmsdl, r2, &
    !             reference_axis, info)
    !integer, parameter :: SP = kind(1.0)
    !integer, parameter :: DP = kind(1.0d0)

    real(DP),  intent(in):: points(:, :)
    real(DP), intent(out):: bending_angles(:)
    real(DP), intent(out):: directions(:, :), origins(:, :)
    real(DP), intent(out):: tilt
    real(DP), optional, intent(in):: reference_axis(:)
    logical, optional , intent(in):: info
    !f2py real(DP), intent(out) bending_angles
    ! real(DP), intent(out)::  radc, rmsdc, rmsdl, r2
    ! real(DP), intent(out):: r2
    
    ! *-- local variables --*
    logical :: quiet = .FALSE.
    real(DP):: direct(3), axvec(3)
    real(DP):: origin(2, 3)
    real(DP) :: twist, height, angle, tmp, bend2
    real(DP) :: ave, var, sdev, mean_abs_std, max_angle
    integer :: nrows, ncols, i, j, dsize, ibend
    integer :: AllocateStatus, DeAllocateStatus
    real(DP), dimension(3) :: ref_axis = (/0.0, 0.0, 1.0/)

    ! *-- dynamic local variables --*
    real(DP), dimension(:), allocatable :: twists, heights
    integer(DP), allocatable:: bending_angles_matrix(:, :)

    ! *-- output format --*
!16  format (1x,'MATRIX M(I,J) FOR ANGLES BETWEEN LOCAL HELIX AXES I AND J')
17  format (i3,1x,40i3)
19  format(1x,'AVERAGE TWIST (DEG.)',6x,f8.2,2x,'S.D.',f8.2,5x,'ADV',f8.2) 
!20  format(1x,'AVERAGE N',17x,f8.2,2x,'S.D.',f8.2,5x,'ADV',f8.2) 
21  format(1x,'AVERAGE UNIT HEIGHT (ANG.)',f8.2,2x,'S.D.',f8.2,5x,'ADV',f8.2) 
27  format(1x,'MEAN BENDING ANGLE',8x,f8.3,2x,'S.D.',f8.3,5x,'ADV',f8.3) 
!28  format(1x,'SUMMARY OF PARAMETERS OBTAINED BY HELANAL'/) 
!29  format(1x,'MEAN BENDING ANGLE (DEG.)',15x,f8.3,2x,'S.D.',f8.3,2x,'ADV',f8.3) 
30  format(1x,'MAXIMUM BENDING ANGLE (DEG.)',12x,f8.3)
!37  format(1x,'AVERAGE UNIT HEIGHT (ANG.)',14x,f8.3,2x,'S.D.',f8.3,2x,'ADV',f8.3)
    ! *-- processing optional arguments --*
    if( present(info) ) quiet = info
    if ( present(reference_axis) ) ref_axis = reference_axis

    nrows = size(points, 1)
    ncols = size(points, 2)
    dsize = nrows - 3
    ! *-- allocate memory --*
    allocate(twists(dsize), stat=AllocateStatus)
    if (AllocateStatus /= 0) stop 'Failed to allocate memory for twists'
    allocate(heights(dsize), stat=AllocateStatus)
    if (AllocateStatus /= 0) stop 'Failed to allocate memory for heights'
    allocate(bending_angles_matrix(dsize, dsize), stat=AllocateStatus)
    if (AllocateStatus /= 0) stop 'Failed to allocate memory for bending_angles_matrix'
    ! *-- initialize --*
    twists = 0.0
    heights = 0.0

    ! *-- allocate(bending_angles(dsize - 3)) --*
    if (ncols.ne.3) then
      stop 'Wrong shape of the input points'
    end if

    do i=1, dsize
      ! *-- slide along the helix spiral curve --*
      call local_helix(points(i,:), points(i+1,:), &
                       points(i+2, :), points(i+3,:), &
                       direct, origin, twist, height, quiet)
      ! *-- save parameters of local helix --*
      !call report(direct, nout=6)
      !write(*,*) 'shape:', shape(directions)
      !write(*,*) 'shape:', shape(direct)
      directions(i, :) = direct
      origins(i, :) = origin(1, :)
      origins(i+1, :) = origin(2, :)
      twists(i) = twist
      heights(i) = height
    end do

    ! *-- bending angle --*
    ! *-- Local (k--k+3, k+3--k+6) bending angles
    do i=1, dsize - 3 
      angle = dot_product(directions(i,:), directions(i+3, :))
      if (abs(angle - 1.0).le.1.0E-06) angle = 1.0
      ! *-- angle in degree --*
      bending_angles(i) = acos(angle)*180/pi
    end do
    ! *-- Mean bending angle for the whole helix
    call stats(bending_angles, ave, var, sdev, mean_abs_std)
    ! *-- Maximum bending angle in the whole helix
    max_angle = maxval(bending_angles)
    if(.not.(quiet)) then
      write(6, 27) ave, sdev, mean_abs_std
      write(6, 30) max_angle
    end if
    ! *-- local bending angle matrix --*
    do i=1, dsize
      do j=1, dsize
        angle = dot_product(directions(i,:), directions(j, :))
        if (abs(angle - 1.0).le.1.0E-06) angle = 1.0
        ! *-- angle in degree --*
        tmp = acos(angle) * 180/ pi
        ibend = int(tmp)
        bend2 = tmp-float(ibend)
        if (bend2.lt.0.5) then
          bending_angles_matrix(i, j) = int(tmp)
        else
          bending_angles_matrix(i, j) = int(tmp + 1.0)
        endif
      enddo
    enddo
    if (.not.(quiet)) then
      do i=1, dsize
        write(6, 17) i,(bending_angles_matrix(i, j), j=1,dsize)
      end do
      ! *-- Average helical parameters --*
      ! 1. average twist
      call stats(twists, ave, var, sdev, mean_abs_std)
      write(6, 19) ave, sdev, mean_abs_std
      ! 2. average unit height
      call stats(heights, ave, var, sdev, mean_abs_std)
      write(6, 21) ave, sdev, mean_abs_std
    endif
    ! *-- fit to circle --*
    !call fit_circle_lsq(origins, tilt, axvec, radc, rmsdc, rmsdl, r2, ref_axis)
    call fit_circle_lsq(origins, tilt, axvec, ref_axis)
    !if (.not.(quiet)) then
    !  write(6, *) 'fit_circle_lsq fininished'
    !  write(6, *) 'points', shape(points)
    !  write(6, *) 'origins', shape(origins)
    !  write(6, *) 'directions', shape(directions)
    !  write(6, *) 'dsize=', dsize
    !  write(6, *)' twists', shape(twists)
    !  write(6, *)' heights', shape(heights)
    !  write(6, *)' bending_angles_matrix', shape(bending_angles_matrix)
    !  write(6, *) 'best tilt angle', tilt
    !  write(6, *) 'fit vector', axvec
    !  write(6,*) '----'
    !end if
    ! *-- deallocate --*
    if (allocated(twists)) then
      deallocate(twists, stat=DeAllocateStatus)
      if (DeAllocateStatus /= 0) stop 'Failed to release memory for heights'
!      write(*,*) 'Release twists'
!    else
!      write(*,*) 'Failed to release twists'
    end if
    if (allocated(heights)) then
      deallocate(heights, stat=DeAllocateStatus)
      if (DeAllocateStatus /= 0) stop 'Failed to release memory for heights'
!      write(*,*) 'Release heights'
!    else
!      write(*,*) 'Failed to release heights'
    end if

    if (allocated(bending_angles_matrix)) then
      deallocate(bending_angles_matrix, stat=DeAllocateStatus)
      if (DeAllocateStatus /= 0) stop 'Failed to release memory for bending_angles_matrix'
!      write(*,*) 'Release bending_angles_matrix'
!    else
!      write(*,*) 'Failed to release bending_angles_matrix'
    end if
  end subroutine fit

end module helanal
