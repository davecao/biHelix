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

  public:: report, reallocate, error

  interface report
    module procedure report_r1
    module procedure report_r2
  end interface report

  interface reallocate
    module procedure reallocate_rv
    module procedure reallocate_rm
    module procedure reallocate_iv
    module procedure reallocate_im
    module procedure reallocate_hv
  end interface reallocate

  interface array_copy
    module procedure array_copy_d
    module procedure array_copy_i
  end interface array_copy

contains

  subroutine warn(msg)
    character(len=*), intent(in) :: msg
    write (*,*) 'Warn: ', msg
  end subroutine warn

  subroutine info(msg)
    character(len=*), intent(in) :: msg
    write (*,*) 'Info: ', msg
  end subroutine info

  subroutine error(msg)
    character(len=*), intent(in) :: msg
    write (*,*) 'Error: ', msg
    stop 'program terminated.'
  end subroutine error
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
  !*****************************************************************************
  ! Description:
  !   reallocate functions (see nrutil.f90)
  !
  ! Standard format:
  !    f90 and later
  !
  ! Status: 
  !    OK
  ! Arguments:
  !
  ! Returns:
  !   None
  !*****************************************************************************

  function reallocate_rv(p, n)
    real(DP), dimension(:), pointer :: p, reallocate_rv
    integer, intent(in) :: n
    integer :: nold, ierr
    allocate(reallocate_rv(n), stat=ierr)
    if (ierr /= 0) then
      call error('reallocate_rv: failed to allocate memory')
    end if
    if ( .not. associated(p) ) return
    nold = size(p)
    reallocate_rv(1:min(nold,n)) = p(1:min(nold, n))
    deallocate(p)
  end function reallocate_rv

  function reallocate_rm(p, n, m)
    real(DP), dimension(:,:), pointer :: p, reallocate_rm
    integer, intent(in) :: n,m
    integer :: nold, mold, ierr
    allocate(reallocate_rm(n,m), stat=ierr)
    if (ierr /= 0) call &
      error('reallocate_rm: problem in attempt to allocate memory')
    if (.not. associated(p)) return
    nold=size(p,1)
    mold=size(p,2)
    reallocate_rm(1:min(nold,n),1:min(mold,m))=&
      p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  end function reallocate_rm
  
  function reallocate_iv(p, n)
    integer, dimension(:), pointer :: p, reallocate_iv
    integer, intent(in) :: n
    integer :: nold, mold, ierr
    allocate(reallocate_iv(n), stat=ierr)
    if (ierr /= 0) call &
      error('reallocate_iv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold = size(p)
    reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  end function reallocate_iv

  function reallocate_im(p, n, m)
    integer, dimension(:,:), pointer :: p, reallocate_im
    integer, intent(in) :: n,m
    integer :: nold,mold,ierr
    allocate(reallocate_im(n,m), stat=ierr)
    if (ierr /= 0) call &
      error('reallocate_im: problem in attempt to allocate memory')
    if (.not. associated(p)) return
    nold=size(p,1)
    mold=size(p,2)
    reallocate_im(1:min(nold, n),1:min(mold, m))=&
      p(1:min(nold, n),1:min(mold, m))
    deallocate(p)
  end function reallocate_im

  function reallocate_hv(p, n)
    character(1), dimension(:), pointer :: p, reallocate_hv
    integer, intent(in) :: n
    integer :: nold, ierr
    allocate(reallocate_hv(n),stat=ierr)
    if (ierr /= 0) call &
      error('reallocate_hv: problem in attempt to allocate memory')
    if (.not. associated(p)) return
    nold=size(p)
    reallocate_hv(1:min(nold, n))=p(1:min(nold, n))
    deallocate(p)
  end function reallocate_hv

  !*****************************************************************************
  ! Description:
  !   move data functions (see nrutil.f90)
  !
  ! Standard format:
  !    f90 and later
  !
  ! Status: 
  !    OK
  ! Arguments:
  !
  ! Returns:
  !   None
  !*****************************************************************************

  subroutine array_copy_d(src,dest,n_copied,n_not_copied)
    real(DP), dimension(:), intent(in) :: src
    real(DP), dimension(:), intent(out) :: dest
    integer, intent(out) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  end subroutine array_copy_d

  subroutine array_copy_i(src,dest,n_copied,n_not_copied)
    integer, dimension(:), intent(in) :: src
    integer, dimension(:), intent(out) :: dest
    integer, intent(out) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  end subroutine array_copy_i

end module utilities
