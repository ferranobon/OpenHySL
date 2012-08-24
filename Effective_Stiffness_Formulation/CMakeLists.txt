# Name of the project
cmake_minimum_required(VERSION 2.8)
project(SAlgorithm)

set(SAlgorithm_VERSION_MAJOR 0)
set(SAlgorithm_VERSION_MINOR 5)

# There are four build modes
# MPIDebug, MPIRelease, Debug and Release
if( NOT (CMAKE_BUILD_TYPE MATCHES Debug ) AND
    NOT (CMAKE_BUILD_TYPE MATCHES Release ))
  message( FATAL_ERROR
    "CMAKE_BUILD_TYPE must be one of: Debug or Release.")
endif()

# Ask the user for build input
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
  "Choose the build type. Accepted options are: Debug and Release"
  FORCE)

# Make will display verbose output
set(CMAKE_VERBOSE_MAKEFILE ON)

option(PROFILING      "Enable profiling using gprof"            OFF)
option(USE_THREADS    "Use thread implementation when possible" OFF)
option(USE_MKL        "Use Intel MKL libraries"                 OFF)
option(USE_ATLAS      "Use ATLAS implementation"                OFF)
option(USE_GOTOBLAS   "Use GOTOBLAS implementation"             OFF)
option(USE_MPI        "Build the parallel implementation"       OFF)
option(USE_SPARSE     "Build with sparse matrix support"        OFF)
option(BUILD_REALTIME "Build also realtime application"         OFF)
option(WITH_ADWIN     "Build with ADwin support"                OFF)

# Add support for Debian OS
if(NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  include_directories(/usr/include/x86_64-linux-gnu/)
endif()

if(WIN32)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  set(WIN_LIBS "-lws2_32 -lwsock32" )
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(DEBUG TRUE)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(RELEASE TRUE)
endif()

if(DEBUG)
  # -D_POSIX_SOURCE is added in order to make gcc be POSIX compliant
  if( "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS_DEBUG "-std=c99 -pedantic -Wall -W -Wextra -Wmissing-prototypes -Wstrict-prototypes -Wbad-function-cast -Wundef -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -Wstrict-overflow=5 -fshort-enums -fno-common -g3 -D_POSIX_SOURCE -D_BSD_SOURCE")
  endif()
  if( "${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_C_FLAGS_DEBUG "-std=c99 -pedantic -Wall -Wextra -Wmissing-prototypes -Wstrict-prototypes -Wbad-function-cast -Wundef -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -Wstrict-overflow=5 -fshort-enums -fno-common -g -D_POSIX_SOURCE -D_BSD_SOURCE")
  endif()
  if( "${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS_DEBUG "-std=c99 -pedantic -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wwrite-strings -D_POSIX_SOURCE -D_BSD_SOURCE -debug all" )
  endif()
  if("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    set(CMAKE_C_FLAGS_DEBUG "-Minform=inform")
  endif()
  # Add profiling flags.
  if(PROFILING)
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pg")
  endif() 
 list(APPEND COMPILE_FLAGS ${CMAKE_C_FLAGS_DEBUG})
endif()

if(RELEASE)
  set(CMAKE_C_FLAGS_RELEASE "-march=native -O2")
  # Add profiling flags.
  if(PROFILING)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -pg")
  endif()
  list(APPEND COMPILE_FLAGS ${CMAKE_C_FLAGS_RELEASE})
endif()

# Check for different implementations of BLAS and LAPACK and tag them as required.
include("CMakeCustomFiles/CMakeMathLibs.cmake")
set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})

# Common source and header files for all the programs
file(GLOB COMMON_SRCS Common_Files/src/*.c)
include_directories(Common_Files/include)

# Find the ADwin libraries and tag it as required
if(WITH_ADWIN)
  add_definitions(-DADWIN_)
  include("CMakeCustomFiles/CMakeADwin.cmake")
  set(CMAKE_REQUIRED_LIBRARIES ${ADWIN_LIBS})
else(NOT WITH_ADWIN)
  list(REMOVE_ITEM COMMON_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Common_Files/src/RoutinesADwin.c)
  MESSAGE( STATUS "Cmake binary dir ${CMAKE_CURRENT_SOURCE_DIR}")
  MESSAGE( STATUS "Compile Flags = ${COMMON_SRCS}" )
endif()

if(USE_MPI)
  # Rules to build the programs related to the parallel version of the algorithm
  add_subdirectory(Parallel)
else(NOT USE_MPI)
  # Rules to build the programs related to the sequential version of the algorithm.
  add_subdirectory(Sequential)
  add_subdirectory(Sparse)
  if(USE_SPARSE)
    add_definitions(-D_SPARSE_)
    if(NOT USE_MKL)
      message(FATAL_ERROR "Sparse support is only available through the MKL libraries")
    endif()
  endif()
endif()

# Rules to build the programs related to the substructure. Each program will
# be compliant with one of the different protocols: Custom, NSEP and OpenFresco
add_subdirectory(Substructure)

Find_Package(Doxygen)
if (DOXYGEN_FOUND)

  set(SOURCE_DOXYGEN_PATH ${CMAKE_BINARY_DIR}/..)
  set(SOURCE_DOXYGEN "${SOURCE_DOXYGEN_PATH}/Common_Files ${SOURCE_DOXYGEN_PATH}/Substructure ${SOURCE_DOXYGEN_PATH}/Sequential")

  configure_file( ${CMAKE_BINARY_DIR}/../Documentation/SAlgorithm.doxyfile ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)

  set(DOXYGEN_INPUT ${CMAKE_BINARY_DIR}/Doxyfile)
  set(DOC_DIR ${CMAKE_BINARY_DIR}/Doc CACHE PATH "Documentation directory")
  set(DOXYGEN_OUTPUT ${DOC_DIR}/html/index.html)

  file(MAKE_DIRECTORY ${DOC_DIR})

  ADD_CUSTOM_COMMAND(
    OUTPUT ${DOXYGEN_OUTPUT}
    COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
    COMMAND ${CMAKE_COMMAND} -E echo "Done."
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Doc
    DEPENDS ${DOXYGEN_INPUT}
    )

  ADD_CUSTOM_TARGET(doc ALL DEPENDS ${DOXYGEN_OUTPUT})
endif (DOXYGEN_FOUND)

MESSAGE( STATUS )
MESSAGE( STATUS "----------------------------------------------------------------------------")
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "Compile Flags = ${COMPILE_FLAGS}" )
MESSAGE( STATUS "Enable profiling using gprof: ${PROFILING}" )
MESSAGE( STATUS "Use thread implementation when possible: ${USE_THREADS}")
MESSAGE( STATUS "Use Reference implementation (Netlib): ${REFERENCE_FOUND}")
MESSAGE( STATUS "Use Intel MKL libraries: ${USE_MKL}")
MESSAGE( STATUS "Use ATLAS implementation: ${USE_ATLAS}")
MESSAGE( STATUS "Use GOTOBLAS implementation: ${USE_GOTOBLAS}")
if(REFERENCE_FOUND)
  MESSAGE( STATUS "WARNING. Using Reference implementation (Netlib). Possible performance drop.")
endif()
MESSAGE( STATUS "Use Parallel implementation: ${USE_MPI}")
MESSAGE( STATUS "Building also realtime applicaton: ${BUILD_REALTIME}" )
MESSAGE( STATUS "Building with ADwin support: ${WITH_ADWIN}" )
MESSAGE( STATUS "System Name: ${CMAKE_SYSTEM}")
MESSAGE( STATUS "----------------------------------------------------------------------------")
MESSAGE( STATUS )