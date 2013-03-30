file( GLOB SUBS_SRCS src/*.c )

if( NOT USE_MPI )
  list(REMOVE_ITEM SUBS_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/Substructure_Auxiliary_MPI.c)
  list(REMOVE_ITEM SUBS_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/Substructure_CouplingNodes_MPI.c)
endif()

if( NOT WITH_ADWIN )
  list(REMOVE_ITEM SUBS_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/Substructure_ADwinRoutines.c)
endif()

add_library( substructure ${SUBS_SRCS} )
target_link_libraries(substructure ${MATH_LIBS} -lm core)