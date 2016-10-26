# Define the include directories
include_directories(include)

# MPI libraries are required.
set(MPI TRUE)

# Define the source files required by the algorithm.
file(GLOB PALGORITHM_SRCS src/*.c)
add_definitions(-D_MPI_)
# Rules for building the main part of the algorithm
add_executable(Parallel ${PALGORITHM_SRCS} ${COMMON_SRCS})
target_link_libraries(Parallel ${MATH_LIBS} ${ADWIN_LIBS})