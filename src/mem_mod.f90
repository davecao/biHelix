! ******************************************************************************
!
! file: mem_mod.f90
!
!
! author: Cao Wei
! Timestamp: Sun Jul  7 19:42:02 2019
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

! **************************************************
!  Membrane infomation
!   1. upper center and its normal vector
!   2. lower center and its norma vector
! **************************************************
module mem_mod
  !****************
  ! other modules
  !****************
  use precision
  implicit none ! use strong type

  private  ! hide the type-bound procedure implementation procedures

  public :: membraneInfo

  type membraneInfo
    real(DP), dimension(3) :: upcenter = (/999.0, 0.0, 0.0/)
    real(DP), dimension(3) :: lowcenter = (/999.0, 0.0, 0.0/)
    real(DP), dimension(3) :: upNormVec = (/999.0, 0.0, 0.0/)
    real(DP), dimension(3) :: lowNormVec = (/999.0, 0.0, 0.0/)
    logical :: doubleLayers = .FALSE.
    logical :: upperOnly = .FALSE.
    logical :: lowerOnly = .FALSE.
    contains
      procedure :: getRefAxis
      procedure :: setLayerNum
      procedure :: showMemInfo
  end type membraneInfo

contains
  !*********************************************
  ! membraneInfo - getRefAxis
  !*********************************************

  function getRefAxis(this) result(ref)
    class(membraneInfo), intent(in) :: this
    real(DP), dimension(3) :: ref
    real(DP) :: tolerance = 0.1e-7
    if (all(abs(this%upNormVec - this%lowNormVec) < tolerance))then
      ! double layers
       if (this%upNormVec(1) > 999.0 .or. this%upNormVec(1) < 999.0) then
        ref = this%upNormVec
      else
        stop "membrane info is not initialized correctly."
      end if
    else
      ! a single layer only
       if ((this%lowNormVec(1) > 999.0 .or. this%lowNormVec(1) < 999.0) .and. &
           (this%upNormVec(1) > 999.0 .or. this%upNormVec(1) < 999.0)) then
        ! the normal is not consistent
        stop "The two normal vectors of the membrane is not consistent"
      else
        ! select a vector normal
         if (this%upNormVec(1) > 999.0 .or. this%upNormVec(1) < 999.0) then
          ref = this%upNormVec
        else
          ref = this%lowNormVec
        end if
      end if
    end if
  end function getRefAxis

  !*********************************************
  ! membraneInfo - setLayerNum
  !*********************************************

  subroutine setLayerNum(this, hasUpper, hasLower)
    class(membraneInfo), intent(inout) :: this
    logical, intent(in):: hasUpper, hasLower

    if( hasUpper .and. hasLower) then
      this%doubleLayers = .TRUE.
      this%upperOnly = .TRUE.
      this%lowerOnly = .TRUE.
    else
      if (hasUpper) then
        this%upperOnly = .TRUE.
      else
        this%lowerOnly = .TRUE.
      end if
    end if

  end subroutine setLayerNum
  !*********************************************
  ! membraneInfo - showMemInfo
  !*********************************************

  subroutine showMemInfo(this)
    class(membraneInfo), intent(inout) :: this
    write(*, *) "Double: ", this%doubleLayers
    write(*, *) "upper: ", this%upperOnly
    write(*, *) "lower: ", this%lowerOnly
  end subroutine

end module mem_mod
