# Define the include directories
include_directories(include)

# MPI libraries are not required.
set(MPI FALSE)

# Define the source files required by the algorithm.
file(GLOB SPARSE_SALGORITHM_SRCS src/*.c)

# Rules for building the main part of the algorithm
add_executable(Sparse_Sequential ${SPARSE_SALGORITHM_SRCS} ${COMMON_SRCS})
target_link_libraries(Sparse_Sequential ${MATH_LIBS} ${ADWIN_LIBS})

if( BUILD_REALTIME )
  add_definitions(-DREAL_TIME_)
  add_executable(Sparse_Sequential-rt ${SPARSE_SALGORITHM_SRCS} ${COMMON_SRCS})
  target_link_libraries(Sparse_Sequential ${MATH_LIBS} ${ADWIN_LIBS})
endif()