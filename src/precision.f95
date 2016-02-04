module precision

#ifdef QUAD_PRECISION
  integer, parameter :: DP = 16 ! kind(1.0d0)
#elif TEN_DIGIT_PRECISION
  integer, parameter :: DP = selected_real_kind(10) !kind(1.0d0)
#else
  integer, parameter :: DP = 8 ! kind(1.0d0)
#endif
  integer, parameter :: kr4 = selected_real_kind(6,37)       ! single precision real
  integer, parameter :: kr8 = selected_real_kind(15,307)     ! double precision real
  
  ! Integer kinds
  
  integer, parameter :: ki4 = selected_int_kind(9)           ! single precision integer
  integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer
  
  !Complex kinds
  
  integer, parameter :: kc4 = kr4                            ! single precision complex
  integer, parameter :: kc8 = kr8                            ! double precision complex

  ! *-- const parameters --*
  real(DP), parameter :: pi = 3.141592653589793238462643383279502884197

end module precision
