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
option(USE_ATLAS      "Use ATLAS implementation"                OFF)
option(USE_GOTOBLAS   "Use GOTOBLAS implementation"             OFF)
option(USE_MKL        "Use Intel MKL libraries"                 OFF)
option(USE_MPI        "Build the parallel implementation"       OFF)
option(USE_OPENMP     "Use OpenMP to speed-up for loops"        OFF)
option(USE_SPARSE     "Build with sparse matrix support"        OFF)
option(USE_THREADS    "Use thread implementation when possible" OFF)
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
if (HDF5_FOUND)
    add_definitions(-D_HDF5_)
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${HDF5_C_FLAGS}")
endif()

if( USE_OPENMP )
Find_Package(OpenMP)
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

# Find the ADwin libraries and tag it as required
if(WITH_ADWIN)
  add_definitions(-D_ADWIN_)
  include("CMakeCustomFiles/CMakeADwin.cmake")
  set(CMAKE_REQUIRED_LIBRARIES ${ADWIN_LIBS})
endif()

# Build the core library
file(GLOB CORE_SRCS src/New_State.c src/Conf_Parser.c src/EffK_Formulation.c src/Err_PID.c src/GainMatrix.c src/Rayleigh.c src/Input_Load.c src/MatrixVector.c src/MatrixVector_IOoperations.c src/Print_Messages.c src/Auxiliary_Math.c src/MatrixVector_PS.c src/MatrixVector_IOoperations_PS.c)

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
  include_directories( ${MMPATH} )
else()
  MESSAGE( STATUS "WARNING. MatrixMarket routines not found in ${MMPATH}. Please, specify the MatrixMarket path manually!" )
endif()

include_directories( include )

if(USE_MKL)
  add_definitions(-D_MKL_)
  if(USE_SPARSE)
    add_definitions(-D_SPARSE_)
  endif()
endif()

add_subdirectory(libcore)

if( USE_MPI )
  add_definitions(-D_MPI_)
  add_subdirectory(libcore_mpi)
endif()

if( USE_SPARSE )
  add_subdirectory(libcore_sp)
endif()

add_subdirectory(libsubstructure)

add_executable( SAlgorithm src/SAlgorithm.c src/Algorithm_Aux.c src/HDF5_Operations.c)# ${CORE_SRCS} ${SPARSE_SRCS} ${SUBS_SRCS})
target_link_libraries(SAlgorithm core core_sp substructure hdf5 hdf5_hl)

add_executable( SAlgorithm_EffM src/SAlgorithm_EffM.c src/Algorithm_Aux.c src/HDF5_Operations.c)# ${CORE_SRCS} ${SPARSE_SRCS} ${SUBS_SRCS})
target_link_libraries(SAlgorithm_EffM core core_sp substructure hdf5 hdf5_hl)

add_executable( Custom_Server src/Custom_Server.c )
target_link_libraries(Custom_Server core core_sp substructure)

add_executable( Exact_Solution src/Exact_Solution.c src/Algorithm_Aux.c src/HDF5_Operations.c)# ${CORE_SRCS} ${SPARSE_SRCS} ${SUBS_SRCS})
target_link_libraries(Exact_Solution core core_sp substructure hdf5 hdf5_hl)

if( USE_MPI )
  add_executable( SAlgorithm_MPI src/SAlgorithm_MPI.c src/Algorithm_Aux.c src/Algorithm_Aux_MPI.c src/HDF5_Operations.c src/HDF5_Operations_MPI.c )# ${CORE_SRCS} ${SPARSE_SRCS} ${SUBS_SRCS})
  target_link_libraries(SAlgorithm_MPI core core_mpi substructure hdf5 hdf5_hl)
endif()
# Rules to build the programs related to the substructure. Each program will
# be compliant with one of the different protocols: Custom, NSEP and OpenFresco

Find_Package(Doxygen)
if (DOXYGEN_FOUND)

  set(SOURCE_DOXYGEN_PATH ${CMAKE_BINARY_DIR}/..)
  set(SOURCE_DOXYGEN "${SOURCE_DOXYGEN_PATH}/include ${SOURCE_DOXYGEN_PATH}/src ${SOURCE_DOXYGEN_PATH}/libcore/src ${SOURCE_DOXYGEN_PATH}/libcore_mpi/src ${SOURCE_DOXYGEN_PATH}/libcore_sp/src ${SOURCE_DOXYGEN_PATH}/libsubstructure/src")

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