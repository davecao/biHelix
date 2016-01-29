module precision

#ifdef QUAD_PRECISION
  integer, parameter :: DP = 16 ! kind(1.0d0)
#elif TEN_DIGIT_PRECISION
  integer, parameter :: DP = selected_real_kind(10) !kind(1.0d0)
#else
  integer, parameter :: DP = 8 ! kind(1.0d0)
#endif
  ! *-- const parameters --*
  real(DP), parameter :: pi = 3.141592653589793238462643383279502884197
end module precision
