!*******************************************************************************
! Description:
!
! Utitlity functions/ subroutines
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

module utilities
  use precision

  private:: report_r1, report_r2

  public:: report
  interface report
     module procedure report_r1
     module procedure report_r2
  end interface report

contains
  !*****************************************************************************
  ! Description:
  !   print 2d matrix for a given points (array)
  !
  ! Standard format:
  !    f90 and later
  !
  ! Status: 
  !    OK
  ! Arguments:
  !   mat (array): Nxd elements, N samples and d dimensions
  !
  ! Returns:
  !   None
  !*****************************************************************************
  subroutine report_r2(mat, nout, fmt, msg)
    !integer, parameter :: SP = kind(1.0)
    !integer, parameter :: DP = kind(1.0d0)
    real(DP), dimension(:, :), intent(in):: mat
    integer, optional, intent(in) :: nout
    character(*), optional, intent(in) :: msg
    character(*), optional, intent(in) :: fmt
    character(*), parameter :: lfmt = '(i3,1x,(*(F12.3, 2X)))'
    integer :: nrows, ncols, i
    integer :: unit_num = 6
    if (present(nout)) then
      unit_num = nout
    end if
!    if (present(fmt)) then
!      lfmt = fmt
!    end if
    nrows = size(mat, dim=1)
    ncols = size(mat, dim=2)
    if (present(msg)) then
      write(unit_num, *) msg
    end if
    do i=1, nrows
       write(unit_num, lfmt) i, mat(i, :)
    end do
  end subroutine report_r2
  !*****************************************************************************
  ! Description:
  !   print 1d matrix for a given points (array)
  !
  ! Standard format:
  !    f90 and later
  !
  ! Status: 
  !    OK
  ! Arguments:
  !   mat (array): 1xd elements, 1 samples and d dimensions
  !
  ! Returns:
  !   None
  !*****************************************************************************
  subroutine report_r1(mat, nout, fmt, msg)
    !integer, parameter :: SP = kind(1.0)
    !integer, parameter :: DP = kind(1.0d0)
    real(DP), dimension(:), intent(in):: mat
    integer, optional, intent(in) :: nout
    character(*), optional, intent(in) :: fmt
    character(*), optional, intent(in) :: msg
    character(*), parameter :: lfmt = '(*(F12.3, 2X))'
    integer :: unit_num = 6
    if (present(nout)) then
      unit_num = nout
    end if
!    if (present(fmt)) then
!      lfmt = fmt
!    end if
    if (present(msg)) then
      write(unit_num, *) msg
    end if
    write(unit_num, lfmt) mat
  end subroutine report_r1

  !*****************************************************************************
  ! Description:
  !   print 1d matrix for a given points (array)
  !
  ! Standard format:
  !    f90 and later
  !
  ! Status: 
  !    OK
  ! Arguments:
  !   mat (array): 1xd elements, 1 samples and d dimensions
  !
  ! Returns:
  !   None
  !*****************************************************************************
  subroutine report_r0(mat, nout, fmt, msg)
    real(DP), intent(in):: mat
    integer, optional, intent(in) :: nout
    character(*), optional, intent(in) :: fmt
    character(*), optional, intent(in) :: msg
    character(*), parameter :: lfmt = '(10A, F12.3)'
    integer :: unit_num = 6
    integer :: msg_len

    if (present(nout)) then
      unit_num = nout
    end if
    if (present(msg)) then
      msg_len = len(msg)
      write(unit_num, lfmt) msg, mat
    else
      write(unit_num, '(F12.3)') mat
    end if
  end subroutine report_r0
end module utilities
