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
option(WITH_MM        "Build with MatrixMarket suppoert"        OFF)

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
    set(CMAKE_C_FLAGS_DEBUG "-m64 -std=c99 -pedantic -Wall -W -Wextra -Wmissing-prototypes -Wstrict-prototypes -Wbad-function-cast -Wundef -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -Wstrict-overflow=5 -fshort-enums -fno-common -g3 -D_POSIX_SOURCE -D_BSD_SOURCE")
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
  set(CMAKE_C_FLAGS_RELEASE "-march=native -O2 -D_POSIX_SOURCE -D_BSD_SOURCE")
  # Add profiling flags.
  if(PROFILING)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -pg")
  endif()
  list(APPEND COMPILE_FLAGS ${CMAKE_C_FLAGS_RELEASE})
endif()

# Check for different implementations of BLAS and LAPACK and tag them as required.
include("CMakeCustomFiles/CMakeMathLibs.cmake")
set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})

# Check for HDF5
Find_Package(HDF5)

# Build the core library
file(GLOB CORE_SRCS src/Common_Formulation.c src/Conf_Parser.c src/EffK_Formulation.c src/Err_PID.c src/Init_GainMatrix.c src/Init_Rayleigh.c src/Input_Load.c src/MatrixVector.c src/MatrixVector_IOoperations.c src/Print_Messages.c src/Auxiliary_Math.c)

set( MMPATH "Extern/MatrixMarket" )
find_file( MMIOSRC mmio.c ${MMPATH})
if( NOT MMIOSRC )
  message( STATUS "Could not find mmio.c in ${MMPATH}")
endif()
find_file( MMIOHEA mmio.h ${MMPATH})
if( NOT MMIOHEA )
  message( STATUS "Could not find mmio.h in ${MMPATH}")
endif()

if( MMIOSRC AND MMIOHEA)
  MESSAGE( STATUS "Found MatrixMarket routines and header files in ${MMPATH}." )
  add_definitions(-D_MATRIXMARKET_)
  set( WITH_MM ON )
  list(APPEND CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/${MMPATH}/mmio.c)
  include_directories( ${MMPATH} )
else()
  MESSAGE( STATUS "WARNING. MatrixMarket routines not found in ${MMPATH}. Please, specify the MatrixMarket path manually!" )
endif()

# Find the ADwin libraries and tag it as required
if(WITH_ADWIN)
  add_definitions(-D_ADWIN_)
  include("CMakeCustomFiles/CMakeADwin.cmake")
  set(CMAKE_REQUIRED_LIBRARIES ${ADWIN_LIBS})
else(NOT WITH_ADWIN)
  list(REMOVE_ITEM CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Common_Files/src/RoutinesADwin.c)
  MESSAGE( STATUS "Cmake binary dir ${CMAKE_CURRENT_SOURCE_DIR}")
  MESSAGE( STATUS "Compile Flags = ${COMMON_SRCS}" )
endif()

if(USE_MKL)
  add_definitions(-D_MKL_)
endif()

include_directories(include)
add_library(core ${CORE_SRCS})
target_link_libraries(core ${MATH_LIBS} -lm)

if(USE_MPI)
  # Rules to build the programs related to the parallel version of the algorithm
  add_definitions(-D_MPI_)
else(NOT USE_MPI)
  # Rules to build the programs related to the sequential version of the algorithm.
  if(USE_SPARSE AND USE_MKL)
    add_definitions(-D_SPARSE_)
    file( GLOB SPARSE_SRCS src/*_Sp.c )
    add_library( core_sp ${SPARSE_SRCS})
    target_link_libraries(core_sp ${MATH_LIBS} -lm core )
  endif()
endif()

# Build Substructure Library
file( GLOB SUBS_SRCS src/Substructure.c src/Substructure_Exact.c src/Substructure_UHYDEfbr.c src/Substructure_SimMeasured.c src/Substructure_Experimental.c)
add_library( substructure ${SUBS_SRCS})
target_link_libraries( substructure -lm core )

add_executable( SAlgorithm SAlgorithm.c src/Algorithm_Aux.c)# ${CORE_SRCS} ${SPARSE_SRCS} ${SUBS_SRCS})
target_link_libraries(SAlgorithm core core_sp substructure)
# Rules to build the programs related to the substructure. Each program will
# be compliant with one of the different protocols: Custom, NSEP and OpenFresco


Find_Package(Doxygen)
if (DOXYGEN_FOUND)

  set(SOURCE_DOXYGEN_PATH ${CMAKE_BINARY_DIR}/..)
  set(SOURCE_DOXYGEN "${SOURCE_DOXYGEN_PATH}/include ${SOURCE_DOXYGEN_PATH}/src")

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
MESSAGE( STATUS "Building with MatrixMarket support: ${WITH_MM}" )
MESSAGE( STATUS "System Name: ${CMAKE_SYSTEM}")
MESSAGE( STATUS "----------------------------------------------------------------------------")
MESSAGE( STATUS )