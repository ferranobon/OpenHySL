file( GLOB SPARSE_SRCS src/*.c )

add_library(core_sp ${SPARSE_SRCS})
target_link_libraries(core_sp ${MATH_LIBS} -lm core)