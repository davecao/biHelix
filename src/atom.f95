module ATOM_MOD
  !****************
  ! other modules
  !****************
  use precision
  implicit none ! use strong type

  private  ! hide the type-bound procedure implementation procedures

  public :: atom

  ! user-defined type
  type atom
    !private ! hide the underlying details
    integer:: resnum, atNum
    character(6):: recname
    character(4):: name
    character(3):: resname
    character:: chId, iCode, altLoc
    real(DP):: x, y, z
    real(DP):: bending_angle = 0.0
    real(DP):: d2top = 0.0
    real(DP):: d2bot = 0.0
  contains
    ! operator
    !generic, public :: assignment(=) => assign
    !procedure, private :: assign
    ! type-bound procedures
    procedure :: getCoord

    !procedure :: readline
    !generic, public :: read => readline

    ! operators
    !procedure, pass:: assign ! it is not necessary for f2003
    !generic, public:: assignment(=) => assign
    !generic:: operator(+) => add

    ! procedure, pass:: writef
    ! generic, public:: write => writef
    procedure :: writef
    generic, public :: write => writef 

    !generic, public :: write(formatted) => writef  ! available in fortran 2008
  end type atom

  ! user-defined constructor
  interface atom
    module procedure new_atom
  end interface atom

contains

  !**********************************************
  ! Atom Constructor
  !**********************************************

  function new_atom(recname, name, atNum, resname, resnum, chId, &
                    iCode, altLoc, x, y, z)
    real(DP), intent(in):: x, y, z
    character(*), intent(in):: recname, name
    character(3), intent(in):: resname
    character, intent(in):: chId, iCode, altLoc
    integer, intent(in):: atNum, resnum
    type(atom) new_atom

    new_atom%recname = recname
    new_atom%iCode = iCode
    new_atom%altLoc = altLoc
    new_atom%name = name
    new_atom%atNum = atNum
    new_atom%resname = resname
    new_atom%resnum = resnum
    new_atom%chId = chId
    new_atom%x = x
    new_atom%y = y
    new_atom%z = z
    new_atom%bending_angle = 0.0
    new_atom%d2top = 0.0
    new_atom%d2bot = 0.0
  end function new_atom

  !***********************************************
  ! Atom operator '=''
  !***********************************************
  !subroutine assign(this, from)
  !  class(atom), intent(out):: this
  !  class(atom), intent(in):: from
  !  ! Nonintrinsic assignment
  !  this = from
  !end subroutine assign
  !***********************************************
  ! Atom readline
  !***********************************************
!  subroutine readline(dtv, unit, iotype, v_list, iostat, iomsg)
!    class(atom), intent(inout) :: dtv
!    integer, intent(in) :: unit
!    character(*), intent(in) :: iotype
!    integer, intent(in):: v_list(:)
!    integer, intent(out):: iostat
!    character(*), intent(inout):: iomsg
!    ! * -- input format -- *
!    ! 1. Atom name: A4
!    ! 2. Atom number: I5
!    ! 3. blank: 1X
!    ! 4. Residue name: A3
!    ! 5. Residue number: I4
!    ! 6. Chain Id: A1
!    ! 7. x of CA: x coordinate of CA
!    ! 8. y of CA: y coordinate of CA
!    ! 9. z of CA: z coordinate of CA
!10 format(A4,I5,1x,A3,I4,A1,f8.3,f8.3,f8.3)
!    read(unit, fmt=10, IOSTAT=iostat, IOMSG=iomsg) dtv%name, dtv%atNum, &
!                  dtv%resname, dtv%resnum, dtv%chId, &
!                  dtv%x, dtv%y, dtv%z
!
!  end subroutine readline

  !***********************************************
  ! Atom writef, user-defined derived-type I/O
  !***********************************************

  subroutine writef(dtv, unit, iotype, v_list, iostat, iomsg)
    class(atom), intent(in)    :: dtv       ! Object to write.
    integer, intent(in)        :: unit      ! Internal unit to write to.
    character(*), optional, intent(in)   :: iotype    ! LISTDIRECTED or DTxxx
    integer, optional, intent(in)        :: v_list(:) ! parameters from fmt spec.
    integer, intent(out)       :: iostat    ! non zero on error, etc.
    character(*), intent(inout):: iomsg     ! define if iostat non zero.

20 format(A6,I5,1x,A4,A1,A3,1x,A1,I4,A1,3x,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3)
    write(unit=unit, fmt=20, IOSTAT=iostat,IOMSG=iomsg) dtv%recname, dtv%atNum, &
                    dtv%name, dtv%altLoc, dtv%resname, dtv%chId, dtv%resnum, &
                    dtv%iCode, dtv%x, dtv%y, dtv%z, &
                    dtv%bending_angle, dtv%d2top, dtv%d2bot

  end subroutine writef

  !*********************************************
  ! Atom - getCoord
  ! Return an array of coordinates
  !*********************************************

  function getCoord(this)
    class(atom), intent(in):: this
    real(DP), dimension(3) :: getCoord
    getCoord(1) = this%x
    getCoord(2) = this%y
    getCoord(3) = this%z
  end function getCoord


end module ATOM_MOD
