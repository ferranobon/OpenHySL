# Define the include directories
include_directories(include)

# MPI libraries are not required.
set(MPI FALSE)

# Define the source files required by the algorithm.
file(GLOB SALGORITHM_SRCS src/*.c)

# Rules for building the main part of the algorithm
add_executable(Sequential ${SALGORITHM_SRCS} ${COMMON_SRCS})
target_link_libraries(Sequential ${MATH_LIBS} ${ADWIN_LIBS} -lm)

if( BUILD_REALTIME )
  add_definitions(-DREAL_TIME_)
  add_executable(Sequential-rt ${SALGORITHM_SRCS} ${COMMON_SRCS})
  target_link_libraries(Sequential-rt ${MATH_LIBS} ${ADWIN_LIBS} -lm)
endif()