file( GLOB MPI_SRCS src/*.c )

add_library( core_mpi ${MPI_SRCS} )
target_link_libraries(core_mpi ${MATH_LIBS} -lm core)