module GROUP_MOD

  !****************
  ! other modules
  !****************
  use precision
  use atom_mod
  implicit none ! use strong type

  private  ! hide the type-bound procedure implementation procedures

  public :: group

  ! ---- Group: store an array of atoms ----
  type group
    private
    integer :: num_atoms
    integer :: full
    type(atom), dimension(:), allocatable :: atoms

  contains
    procedure :: add
    procedure :: getCoords
    procedure :: resize => reallocate
    procedure :: printf
  end type group

  interface group
    module procedure new_group
  end interface group

contains
    !**********************************************
  ! Group Constructor
  !**********************************************

  function new_group(natoms)
    integer, optional, intent(in) :: natoms
    integer :: init_natoms = 10
    integer :: AllocateStatus, DeAllocateStatus
    type(group) :: new_group ! an instance of group

    if (present(natoms) .and. natoms > 10) then
      init_natoms = natoms
    end if
    new_group%num_atoms = 0
    new_group%full = init_natoms
    allocate(new_group%atoms(init_natoms), stat=AllocateStatus)
    if (AllocateStatus /= 0) stop 'Failed to allocate memory for atoms'
    call new_group%printf()
  end function new_group

  !**********************************************
  ! Group - realloc the objects of Atom
  !**********************************************

  subroutine reallocate(this)
    class(group), intent(inout) :: this
    !integer, optional, intent(in) :: nsize
    type(atom), dimension(:), allocatable :: tmp
    !class(atom), pointer :: atom_pt
    integer :: AllocateStatus, DeAllocateStatus
    integer :: num_atoms
    integer :: fsize

    num_atoms = this%num_atoms
    fsize = num_atoms * 2

    call this%printf()

    allocate(tmp(fsize), stat=AllocateStatus)
    if (AllocateStatus /= 0) stop 'Atom: Failed to allocate memory'
    write(*,*) 'num_atoms:', num_atoms, 'size: ', size(this%atoms)
    tmp(1:num_atoms) = this%atoms
    !deallocate(this%atoms, stat=DeAllocateStatus)
    !if (DeAllocateStatus /= 0) stop 'Group: Failed to release memory'
    call move_alloc(from=tmp, to=this%atoms)
    deallocate(tmp)
    call this%printf()
  end subroutine reallocate
  !**********************************************
  ! Group - add
  !**********************************************

  subroutine add(this, atom_obj)
    class(group), intent(inout) :: this
    class(atom), intent(in) :: atom_obj
    integer :: natoms
    call this%printf()

    natoms = this%num_atoms
    if (natoms > this%full) then
      ! reallocate
      ! gp_pt=>reallocate(gp_pt, 2*natoms)
      call this%resize()
    end if
    this%atoms(natoms + 1) = atom_obj
    this%num_atoms = natoms + 1
    this%full = 2*natoms
  end subroutine
  !**********************************************
  ! Group - getCoords
  !**********************************************

  function getCoords(this) result(coords)
    class(group), intent(in) :: this
    real(DP), dimension(:, :), allocatable :: coords
    integer :: i, natoms

    natoms = this%num_atoms
    allocate(coords(natoms, 3))
    do i=1, natoms
      coords(i,:) = this%atoms(i)%getCoord()
    end do
  end function getCoords

  !**********************************************
  ! Group - getCoords
  !**********************************************

  subroutine printf(this)
    class(group), intent(in) :: this
    integer :: i
    integer :: iostat
    character(LEN=100) :: iomsg

    write(*,*) 'natoms:', this%num_atoms, 'full:', this%full
    do i=1, this%num_atoms
      call this%atoms(i)%writef(unit=6, iostat=iostat, iomsg=iomsg)
      if ( iostat /= 0 ) then
        write(*,*)'error:', iomsg
      end if
    end do
  end subroutine printf
end module GROUP_MOD