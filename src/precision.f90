! ******************************************************************************
!
! file: precision.f90
!
!
! author: Cao Wei
! Timestamp: Sun Jul  7 19:42:21 2019
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

module precision

#ifdef QUAD_PRECISION
  integer, parameter :: DP = 16 ! kind(1.0d0)
#elif TEN_DIGIT_PRECISION
  integer, parameter :: DP = selected_real_kind(10) !kind(1.0d0)
#else
  integer, parameter :: DP = 8 ! kind(1.0d0)
#endif
  integer, parameter :: kr4 = selected_real_kind(6, 37)       ! single precision real
  integer, parameter :: kr8 = selected_real_kind(15, 307)     ! double precision real

  ! Integer kinds
  integer, parameter :: ki4 = selected_int_kind(9)           ! single precision integer
  integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer

  !Complex kinds

  integer, parameter :: kc4 = kr4                            ! single precision complex
  integer, parameter :: kc8 = kr8                            ! double precision complex

  ! *-- const parameters --*
  real(DP), parameter :: pi = 3.141592653589793238462643383279502884197

  ! *-- stdin/ stdout --*
  integer, parameter :: stdin = 5   ! read from cli
  integer, parameter :: stdout = 6  ! write to screen


end module precision
