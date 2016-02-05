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
    procedure :: getNumAtoms
    procedure :: setAtomBendingAngleAt
    procedure :: printf_db ! for debug
    procedure :: printf_rl ! for both screen and file
    generic :: printf => printf_db, printf_rl
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
    integer :: AllocateStatus
    type(group) :: new_group ! an instance of group

    if (present(natoms) .and. natoms > 10) then
      init_natoms = natoms
    end if
    new_group%num_atoms = 0
    new_group%full = init_natoms
    allocate(new_group%atoms(init_natoms), stat=AllocateStatus)
    if (AllocateStatus /= 0) stop 'Failed to allocate memory for atoms'
  end function new_group

  !**********************************************
  ! Group - realloc the objects of Atom
  !**********************************************

  subroutine reallocate(this)
    class(group), intent(inout) :: this
    !integer, optional, intent(in) :: nsize
    type(atom), dimension(:), allocatable :: tmp
    !class(atom), pointer :: atom_pt
    integer :: AllocateStatus
    integer :: full
    integer :: fsize

    full = this%full
    fsize = full * 2

    allocate(tmp(fsize), stat=AllocateStatus)
    if (AllocateStatus /= 0) stop 'Atom: Failed to allocate memory'
    tmp(1:full) = this%atoms(1:full)
    !deallocate(this%atoms, stat=DeAllocateStatus)
    !if (DeAllocateStatus /= 0) stop 'Group: Failed to release memory'
    call move_alloc(from=tmp, to=this%atoms)
    this%full = fsize
    !deallocate(tmp)
  end subroutine reallocate
  !**********************************************
  ! Group - add
  !**********************************************

  subroutine add(this, atom_obj)
    class(group), intent(inout) :: this
    class(atom), intent(in) :: atom_obj
    integer :: natoms
    natoms = this%num_atoms
    
    if (natoms >= this%full) then
      ! reallocate
      call this%resize()
    end if
    this%atoms(natoms+1) = atom_obj
    this%num_atoms = natoms + 1
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
  ! Group - getNumAtoms
  !**********************************************

  function getNumAtoms(this) result(natoms)
    class(group), intent(in) :: this
    integer :: natoms
    natoms = this%num_atoms
  end function getNumAtoms

  !**********************************************
  ! Group - setAtomBendingAngleAt
  !**********************************************

  subroutine setAtomBendingAngleAt(this, pos, angle)
    class(group), intent(inout) :: this
    real(DP), intent(in) :: angle
    integer, intent(in) :: pos

    this%atoms(pos)%bending_angle = angle
  end subroutine setAtomBendingAngleAt
  !**********************************************
  ! Group - printf interface
  !**********************************************

  subroutine printf_db(this)
    class(group), intent(in) :: this
    integer :: i
    integer :: iostat
    character(LEN=100) :: iomsg

    !write(*,*) 'natoms:', this%num_atoms, 'full:', this%full
    do i=1, this%num_atoms
      call this%atoms(i)%writef(unit=6, iostat=iostat, iomsg=iomsg)
      if ( iostat /= 0 ) then
        write(*,*)'error:', iomsg
      end if
    end do
  end subroutine printf_db

  subroutine printf_rl(dtv, unit, iotype, v_list, iostat, iomsg)
    class(group), intent(in) :: dtv
    integer, intent(in)        :: unit      ! Internal unit to write to.
    character(*), optional, intent(in)   :: iotype    ! LISTDIRECTED or DTxxx
    integer, optional, intent(in)        :: v_list(:) ! parameters from fmt spec.
    integer, intent(out)       :: iostat    ! non zero on error, etc.
    character(*), intent(inout):: iomsg     ! define if iostat non zero.
    integer :: i

    do i=1, dtv%num_atoms
      call dtv%atoms(i)%writef(unit=unit, iostat=iostat, iomsg=iomsg)
      if ( iostat /= 0 ) then
        write(*,*)'error:', iomsg
      end if
    end do
  end subroutine printf_rl

end module GROUP_MOD
