# Define the include directories
include_directories(include)

# MPI libraries are not required.
set(MPI FALSE)

 if(USE_SPARSE)
   add_definitions(-D_SPARSE_)
   if(NOT USE_MKL)
     message(FATAL_ERROR "Sparse support is only available through the MKL libraries")
   endif()
 endif()


# Define the source files required by the algorithm.
file(GLOB SPARSE_SALGORITHM_SRCS src/*.c)

# Rules for building the main part of the algorithm
add_executable(SAlgorithm ${SPARSE_SALGORITHM_SRCS} ${COMMON_SRCS})
target_link_libraries(SAlgorithm ${MATH_LIBS} ${ADWIN_LIBS} -lm -lhdf5 -lhdf5_hl)

if( BUILD_REALTIME )
  add_definitions(-DREAL_TIME_)
  add_executable(SAlgorithm-rt ${COMMON_SRCS} ${SPARSE_SALGORITHM_SRCS} )
  target_link_libraries(SAlgorithm ${MATH_LIBS} ${ADWIN_LIBS} -lm)
endif()