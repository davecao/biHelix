!******************************************************************************
! Description:
!
!  memory: 
!
!******************************************************************************

module memory

  use precision

  private:: hpAllocMemR1, hpAllocMemR2

  public:: hpAllocMem
!  public:: hpDeAllocMem

  interface hpAllocMem
     module procedure hpAllocMemR1
     module procedure hpAllocMemR2
  end interface hpAllocMem

!  interface hpDeAllocMem
!     module procedure hpDeAllocMem
!  end interface hpDeAllocMem

contains
  !**********************************************************
  ! Description:
  !    allocate a rank-1 array
  ! Return
  !    a pointer
  !**********************************************************
  function hpAllocMemR1(in, sizeofIn) result(out)
    real(DP), target, allocatable, dimension(:):: in
    integer, intent(in):: sizeofIn
    integer :: AllocateStatus
    real(DP), pointer, dimension(:):: out

    allocate(in(sizeofIn), stat=AllocateStatus)
    if (AllocateStatus /= 0) stop 'Failed to allocate memory'
    out => in
  end function hpAllocMemR1

  !**********************************************************
  ! Description:
  !    allocate a rank-2 array
  ! Return
  !    a pointer
  !**********************************************************
  function hpAllocMemR2(in, nrows, ncols) result(out)
    real(DP), target, allocatable, dimension(:, :)::in
    integer, intent(in):: nrows, ncols
    integer :: AllocateStatus
    real(DP), pointer, dimension(:, :):: out

    allocate(in(nrows, ncols), stat=AllocateStatus)
    if (AllocateStatus /= 0) stop 'Failed to allocate memory for rank(2) array'
    out => in
  end function hpAllocMemR2

  !**********************************************************
  ! Description:
  !    Deallocate array
  ! Return
  !    None
  !**********************************************************
  
end module memory
