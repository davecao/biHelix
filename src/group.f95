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
  use mem_mod
  use desp_mod
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
    real(DP), public, dimension(3) :: reference_axis, helix_axis
    type(membraneInfo), public :: mem_info

    !real(DP), public, dimension(3) :: upper
    !real(DP), public, dimension(3) :: lower
    !real(DP), public, dimension(3) :: mem_normal
    character(4), public :: idCode                    ! HEADER
    character(40), public :: classification           ! HEADER
    character(9), public :: depDate                   ! HEADER
    type(atom), dimension(:), allocatable :: atoms    ! ATOM
  contains
    procedure :: add
    procedure :: getCoords
    procedure :: resize => reallocate
    procedure :: getNumAtoms
    procedure :: getReferenceAxis
    procedure :: findDistToMem
    procedure :: setAtomBendingAngleAt
    procedure :: setMemLayerNum
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
    new_group%mem_info%upcenter = (/999.0, 0.0, 0.0/)
    new_group%mem_info%lowcenter = (/999.0, 0.0, 0.0/)
    new_group%mem_info%upNormVec = (/999.0, 0.0, 0.0/)
    new_group%mem_info%lowNormVec = (/999.0, 0.0, 0.0/)
    new_group%helix_axis = (/999.0, 0.0, 0.0/)
  end function new_group

  !**********************************************
  ! Group - getReferenceAxis
  !**********************************************

  function getReferenceAxis(this) result(ref_ax)
    class(group), intent(in) :: this
    real(DP), dimension(3) :: ref_ax
    ref_ax = this%mem_info%getRefAxis()
  end function getReferenceAxis

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
  ! Group - setMemLayerNum
  !**********************************************

  subroutine setMemLayerNum(this, hasUpper, hasLower)
    class(group), intent(inout) :: this
    logical, intent(in) :: hasUpper, hasLower
    call this%mem_info%setLayerNum(hasUpper=hasUpper, hasLower=hasLower)
  end subroutine setMemLayerNum
  !**********************************************
  ! Group - findDistToMem
  !**********************************************

  subroutine findDistToMem(this)
    class(group), intent(inout) :: this
    integer :: natoms, i
    natoms = this%num_atoms

    if (this%mem_info%doubleLayers) then
      caloop: do i=1, natoms
        call this%atoms(i)%setDistToMem(this%mem_info%upcenter, &
                                   this%mem_info%upNormVec, &
                                   this%mem_info%lowcenter, &
                                   this%mem_info%lowNormVec)
      end do caloop
    else if (this%mem_info%upperOnly) then
      do i=1, natoms
        call this%atoms(i)%setDistToMem(this%mem_info%upcenter, &
                                   this%mem_info%upNormVec, 0)
      end do
    else if (this%mem_info%lowerOnly) then
      do i=1, natoms
        call this%atoms(i)%setDistToMem(this%mem_info%lowcenter, &
                                   this%mem_info%lowNormVec, 1)
      end do
    else
    end if
  end subroutine findDistToMem

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

    call report(this%mem_info%upcenter, nout=unit, msg="The center of upper layer:")
    call report(this%mem_info%lowcenter, nout=unit, msg="The center of lower layer:")
    call report(this%mem_info%upNormVec, nout=unit, msg="The normal of upper layer:")
    call report(this%mem_info%lowNormVec, nout=unit, msg="The normal of lower layer:")

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
    character(len=80) :: st1, ed1
    integer :: i
    call glue(fp, "REM   ")
    write(fp,'(A6,4x,A40,A9,3x,A4)') "HEADER", this%classification, &
                                     this%depDate, this%idCode
    ! *-- write parameters --*
    write(fp, '(A6,1x,3F8.3)') 'REFAXS', this%reference_axis
    write(fp, '(A6,1x,3F8.3)') 'MEMUPP', this%mem_info%upcenter
    write(fp, '(A6,1x,3F8.3)') 'MEMUNM', this%mem_info%upNormVec
    write(fp, '(A6,1x,3F8.3)') 'MEMLOW', this%mem_info%lowcenter
    write(fp, '(A6,1x,3F8.3)') 'MEMLNM', this%mem_info%lowNormVec

    do i=1, size(this%directions, 1)
      st1 = this%atoms(i)%getIdentifier()
      ed1 = this%atoms(i+3)%getIdentifier()
      write(fp, '(A6,1x,A,A3,A,3X,3F8.3)') 'DIRECT', trim(st1), " - ", &
                trim(ed1), this%directions(i,:)
    end do
    do i=1,size(this%helix_origins, 1)
      st1 = this%atoms(i)%getIdentifier()
      ed1 = this%atoms(i+1)%getIdentifier()
      write(fp, '(A6,1x,A,A3,A,3X,3F8.3)') 'HELORG', trim(st1), " - ", &
                trim(ed1), this%helix_origins(i,:)
    end do
    write(fp, '(A6, 1x, 3F12.3)') "HELAXS", this%helix_axis
    write(fp, '(A6, 1x, F12.3)') "TILANG", this%tilt
    !call report(this%directions, nout=fp, msg="Directions:")
    !call report(this%helix_origins, nout=fp, msg="Helix origins:")
    !call report(this%reference_axis, nout=fp, msg="Reference axis:")

    !call report(this%mem_info%upcenter, nout=fp, msg="The center of upper layer:")
    !call report(this%mem_info%lowcenter, nout=fp, msg="The center of lower layer:")
    !call report(this%mem_info%upNormVec, nout=fp, msg="The normal of upper layer:")
    !call report(this%mem_info%lowNormVec, nout=fp, msg="The normal of lower layer:")
    

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
