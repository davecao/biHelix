########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(HPSRC   ${SRC}/cla.f95
            ${SRC}/kinds.f95
            ${SRC}/stringmod.f90
            ${SRC}/precision.f95
            ${SRC}/interp_hel.f95
            ${SRC}/io.f95
            ${SRC}/desp.f95
            ${SRC}/mem_mod.f95
            ${SRC}/atom.f95
            ${SRC}/group.f95
            ${SRC}/utilities.f95
            ${SRC}/linalg.f95
            ${SRC}/helanal.f95
            ${SRC}/main.f95
)

#if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#  # gfortran
#  set (CMAKE_Fortran_FLAGS_RELEASE "-DQUAD_PRECISION -cpp -funroll-all-loops -fno-f2c -O3")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-DQUAD_PRECISION -cpp -fno-f2c -O0 -g")
#elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
#  # ifort (untested)
#  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
#elseif (Fortran_COMPILER_NAME MATCHES "g77")
#  # g77
#  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
#else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
#  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
#  message ("No optimized Fortran compiler flags are known, we just try -O2...")
#  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
#endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

#set_compile_flag(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#    Fortran "-DQUAD_PRECISION"
#)
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-std=gnu"        # GNU
                 )
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE " ${CMAKE_Fortran_FLAGS_RELEASE})
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG " ${CMAKE_Fortran_FLAGS_DEBUG})
#set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#    Fortran "-cpp"
#            "-DQUAD_PRECISION"
#            "-c"
#            "-fdefault-real-8"
#            "-fbacktrace" 
#            "-fno-align-commons" 
#            "-fbounds-check"
#)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${PROJECT_NAME} ${HPSRC})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links foo to the bar library
#TARGET_LINK_LIBRARIES(${EXECUTABLE} ${BARLIB})

# Uncomment if you need to link to BLAS and LAPACK
TARGET_LINK_LIBRARIES(${PROJECT_NAME} 
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT})

# Uncomment if you have parallization
#IF(USE_OPENMP)
#    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
#                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
#                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
#ELSEIF(USE_MPI)
#    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
#                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
#                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
#    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
#    TARGET_LINK_LIBRARIES(${FOOEXE} ${MPI_Fortran_LIBRARIES})
#ENDIF(USE_OPENMP)

#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
INSTALL(TARGETS ${EXE} RUNTIME DESTINATION bin)
