! * -- output format -- *
! 1. 1-6(A6): 'ATOM  '
! 2. 7-11(I5): Atom serial number 
! 3. 13-16(A4):Atom name
! 4. 17(A1): Alternate location indicator
! 5. 18-20(A3): Residume name
! 6. 21 blank(1x): blank
! 7. 22(A1): chain ID
! 8. 23-26 (I4): Residue sequence number
! 9. 27 (A1): iCode
! 10. 31-38(F8.3): x of CA: x coordinate of CA
! 11. 39-46(F8.3): y of CA: y coordinate of CA
! 12. 47-54(F8.3): z of CA: z coordinate of CA
! 13. bending angle of the residue: f8.3 angle in degree
! 14. distance to upper layer: f8.3 angstroms to the upper layer
! 15. distance to lower layer: f8.3 angstroms to the lower layer
! ************************************************************************
! 20 format(A6,I5,1x,A4,A1,A3,1x,A1,I4,A1,3x,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3)
!

module GROUP_MOD

  !****************
  ! other modules
  !****************
  use precision
  use atom_mod
  use utilities, only:report
  implicit none ! use strong type

  private  ! hide the type-bound procedure implementation procedures

  public :: group

  ! ---- Group: store an array of atoms ----
  type group
    private
    integer :: num_atoms
    integer :: full
    real(DP), public, allocatable:: directions(:, :)
    real(DP), public, allocatable:: helix_origins(:, :)
    real(DP), public :: tilt
    real(DP), public, dimension(3) :: reference_axis
    real(DP), public, dimension(3) :: upper
    real(DP), public, dimension(3) :: lower
    real(DP), public, dimension(3) :: mem_normal
    character(4), public :: idCode                    ! HEADER
    character(40), public :: classification           ! HEADER
    character(9), public :: depDate                   ! HEADER
    type(atom), dimension(:), allocatable :: atoms    ! ATOM

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
    new_group%reference_axis = (/0.0, 0.0, 1.0/)
    new_group%upper = (/0.0, 0.0, 0.0/)
    new_group%lower = (/0.0, 0.0, 0.0/)
    new_group%mem_normal = (/0.0, 0.0, 0.0/)
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
    integer :: iostat, unit
    character(LEN=100) :: iomsg
    unit = 6
    write(unit,fmt='(A6,4x,A40,A9,3x,A4)') "HEADER", this%classification, &
                                     this%depDate, this%idCode
    call report(this%directions, nout=unit, msg="Directions:")
    call report(this%helix_origins, nout=unit, msg="Helix origins:")
    call report(this%reference_axis, nout=unit, msg="Reference axis:")
    call report(this%upper, nout=unit, msg="The center of upper layer:")
    call report(this%lower, nout=unit, msg="The center of lower layer:")
    write(unit, '(A, F12.3)') "Tilt angle w.r.t Reference axis:", this%tilt
    !write(*,*) 'natoms:', this%num_atoms, 'full:', this%full
    do i=1, this%num_atoms
      call this%atoms(i)%writef(unit=unit, iostat=iostat, iomsg=iomsg)
      if ( iostat /= 0 ) then
        write(unit,*)'error:', iomsg
      end if
    end do
  end subroutine printf_db

  subroutine printf_rl(this, fp)
    class(group), intent(in) :: this
    integer, intent(in) :: fp
    integer :: iostat    ! non zero on error, etc.
    character(LEN=100) :: iomsg     ! define if iostat non zero.
    integer :: i

    write(fp,'(A6,4x,A40,A9,3x,A4)') "HEADER", this%classification, &
                                     this%depDate, this%idCode
    ! *-- write parameters --*
    call report(this%directions, fp, msg="Directions:")
    call report(this%helix_origins, fp, msg="Helix origins:")
    call report(this%reference_axis, fp, msg="Reference axis:")
    call report(this%upper, fp, msg="The center of upper layer:")
    call report(this%lower, fp, msg="The center of lower layer:")
    write(fp, '(A, F12.3)') "Tilt angle w.r.t Reference axis:", this%tilt

    ! *-- write data --*
    do i=1, this%num_atoms
      call this%atoms(i)%writef(fp, iostat=iostat, iomsg=iomsg)
      if ( iostat /= 0 ) then
        write(*,*)'error:', iomsg
      end if
    end do
    ! *-- close file ---
  end subroutine printf_rl

end module GROUP_MOD
