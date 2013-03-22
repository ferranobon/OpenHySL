file( GLOB CORE_SRCS src/*.c )

if( WITH_MM )
  list(APPEND CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/../${MMPATH}/mmio.c)
endif()

add_library( core ${CORE_SRCS} )
target_link_libraries(core ${MATH_LIBS} -lm)