# Define the include directories
include_directories(include)

# Define the common sources for all the programs dealing with the substructure
#set(COMMON_SUB_SRCS src/RoutinesADwin.c)
#list(APPEND COMMON_SUB_SRCS src/Substructure.c)

set(COMMON_SUB_SRCS src/Substructure.c)

# Define whether the Substructure will be simulated or using ADwin
if (TEST)
  add_definitions(-DSIMULATE_SUB_)
endif()
# Rules for building the Custom Substructure server
add_executable(Substructure-Custom ${COMMON_SRCS} ${COMMON_SUB_SRCS} src/Custom_Server.c)
target_link_libraries(Substructure-Custom ${MATH_LIBS} ${ADWIN_LIBS} ${WIN_LIBS} m)

# Rules for building a NSEP compliant FCM
add_executable(Substructure-NSEP ${COMMON_SRCS} ${COMMON_SUB_SRCS} src/NSEP_Dorka_FCM.c)
target_link_libraries(Substructure-NSEP ${MATH_LIBS} ${ADWIN_LIBS} ${WIN_LIBS} m)

# Rules for building an OpenFresco compliant code
add_executable(Substructure-OP ${COMMON_SRCS} ${COMMON_SUB_SRCS} src/OpenFresco_Dorka.c)
target_link_libraries(Substructure-OP ${MATH_LIBS} ${ADWIN_LIBS} ${WIN_LIBS} m)