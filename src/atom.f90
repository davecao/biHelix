! ******************************************************************************
!
! file: atom.f90
!
!
! author: Cao Wei
! Timestamp: Sat Jul  6 07:52:01 2019
!
! Copyright (C) 2019 Cao Wei. All rights reserved.
!
!
! The following statement of license applies *only* to this header file,
! and *not* to the other files distributed with FFTW or derived therefrom:
!
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
! OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
! GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! ******************************************************************************

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
    integer:: resnum
    integer:: atNum
    character(6):: recname
    character(4):: name
    character(3):: resname
    character(2):: ss
    character:: chId
    character:: iCode
    character:: altLoc
    real(DP):: x
    real(DP):: y
    real(DP):: z
    real(DP):: bending_angle = 0.0
    real(DP):: d2top = 0.0
    real(DP):: d2bot = 0.0
    ! *-- stride --*
    ! keywaords:
    ! extended --> secondary E
    ! helix    --> secondary H
    ! helix310 --> secondary G
    ! helixpi  --> secondary I
    ! turn     --> secondary T
    ! bridge   --> secondary B
    ! bend     --> secondary S
    ! coil     --> secondary C
    real(DP):: phi
    real(DP):: psi
    real(DP):: area

  contains
    ! operator
    !generic, public :: assignment(=) => assign
    !procedure, private :: assign
    ! type-bound procedures
    procedure :: getCoord
    procedure :: getIdentifier
    procedure :: distToDoubleMem
    procedure :: distToSingleMem
    generic, public:: setDistToMem=>distToDoubleMem, distToSingleMem
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
                    iCode, altLoc, x, y, z, phi, psi, area, ss)
    real(DP), intent(in):: x, y, z
    real(DP), intent(in):: phi, psi, area
    character(*), intent(in) :: ss
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
    new_atom%phi = phi
    new_atom%psi = psi
    new_atom%area = area
    new_atom%ss = ss
    new_atom%bending_angle = 0.0
    new_atom%d2top = 999.0
    new_atom%d2bot = 999.0
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

20 format(A6,I5,1x,A4,A1,A3,1x,A1,I4,A1,3x,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3,A2)
    write(unit=unit, fmt=20, IOSTAT=iostat,IOMSG=iomsg) dtv%recname, dtv%atNum, &
                    dtv%name, dtv%altLoc, dtv%resname, dtv%chId, dtv%resnum, &
                    dtv%iCode, dtv%x, dtv%y, dtv%z, &
                    dtv%bending_angle, dtv%d2top, dtv%d2bot, &
                    dtv%phi, dtv%psi, dtv%area, dtv%ss

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


  !*********************************************
  ! Atom - getIdentifier
  ! Return a string, resname resnum chId
  !*********************************************

  function getIdentifier(this) result(atId)
    class(atom), intent(in):: this
    character(len=80), allocatable :: atId
    character(len=80) :: identifier
    write(identifier,'(I5,1X,A3,1X,A1)') this%resnum, this%resname, this%chId
    atId = adjustl(trim(identifier))
  end function getIdentifier

  !*********************************************
  ! Atom - distToSingleMem
  !*********************************************

  subroutine distToSingleMem(this, center, normal, label)
    class(atom), intent(inout):: this
    real(DP), dimension(3), intent(in) :: center
    real(DP), dimension(3), intent(in) :: normal
    integer, intent(in) :: label
    real(DP), dimension(3) :: coords
    real(DP), dimension(3) :: normal_vec
    real(DP) :: dist
    ! label 0 --> upper
    ! label 1 --> lower
    normal_vec = normal/norm2(normal)

    coords = this%getCoord()
    dist = dot_product(coords - center, normal_vec)
    if (label == 0) then
      this%d2top = dist
    else if (label == 1) then
      this%d2bot = dist
    else
      stop "setDistToMem: label should be 0 or 1"
    end if
  end subroutine distToSingleMem
  !*********************************************
  ! Atom - distToDoubleMem
  !*********************************************

  subroutine distToDoubleMem(this, up_center, up_norm, low_center, low_norm)
    class(atom), intent(inout):: this
    real(DP), dimension(3), intent(in) :: up_center
    real(DP), dimension(3), intent(in) :: up_norm
    real(DP), dimension(3), intent(in) :: low_center
    real(DP), dimension(3), intent(in) :: low_norm
    real(DP), dimension(3) :: normal_vec
    normal_vec = up_norm/norm2(up_norm)

    call this%distToSingleMem(center=up_center, normal=normal_vec, label=0)
    normal_vec = low_norm/norm2(low_norm)
    call this%distToSingleMem(center=low_center, normal=normal_vec, label=1)
  end subroutine distToDoubleMem
end module ATOM_MOD
