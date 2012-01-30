# find the libraries required by MKL

set(INTEL_ROOT /opt/intel)
set(INTEL_LIB_SUBDIR /lib/intel64 /lib/em64t)
set(MKL_ROOT /opt/intel/mkl)
set(MKL_LIB_SUBDIR /lib/intel64 /lib/em64t)

find_library(MKL_INTEL_LIB
  NAMES mkl_intel_lp64
  PATHS ${MKL_ROOT}
  PATH_SUFFIXES ${MKL_LIB_SUBDIR})
find_library(MKL_CORE_LIB
  NAMES mkl_core
  PATHS ${MKL_ROOT}
  PATH_SUFFIXES ${MKL_LIB_SUBDIR})

set(MKL_REQUIRED MKL_INTEL MKL_CORE)

if(NOT USE_THREADS)
  find_library(MKL_SEQUENTIAL_LIB
    NAMES mkl_sequential
    PATHS ${MKL_ROOT}
    PATH_SUFFIXES ${MKL_LIB_SUBDIR})
  list(APPEND MKL_REQUIRED MKL_SEQUENTIAL)
endif()

if(USE_THREADS AND "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")

  find_package(OpenMP)
  find_package(Threads)
  find_library(MKL_GNU_THREAD_LIB
    NAMES mkl_gnu_thread
    PATHS ${MKL_ROOT}
    PATH_SUFFIXES ${MKL_LIB_SUBDIR})
  list(APPEND MKL_REQUIRED MKL_GNU_THREAD)
   
elseif(USE_THREADS AND "${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  find_package(OpenMP)
  find_package(Threads)
  find_library(MKL_INTEL_THREAD_LIB
    NAMES mkl_intel_thread
    PATHS ${MKL_ROOT}
    PATH_SUFFIXES ${MKL_LIB_SUBDIR})
  list(APPEND MKL_REQUIRED MKL_INTEL_THREAD)

elseif(USE_THREADS AND "${CMAKE_C_COMPILER_ID}" STREQUAL "pgi")
  find_package(OpenMP)
  find_package(Threads)
  find_library(PGI_OPENMP_LIB NAMES pgmp )
  find_library(MKL_PGI_THREAD_LIB
    NAMES mkl_pgi_thread
    PATHS ${MKL_ROOT}
    PATH_SUFFIXES ${MKL_LIB_SUBDIR})
  find_library(PGF90_LIB NAMES pgf90libs)
  list(APPEND MKL_REQUIRED PGI_OPENMP MKL_PGI_THREAD PGF90)

elseif(USE_THREADS AND "${CMAKE_C_COMPILER_ID}" STREQUAL "PathScale")
  find_package(Threads)
  find_library(MKL_INTEL_THREAD_LIB
    NAMES mkl_intel_thread
    PATHS ${MKL_ROOT}
    PATH_SUFFIXES ${MKL_LIB_SUBDIR})
  find_library(MKL_INTEL_OPENMP_LIB
    NAMES iomp5
    PATHS ${INTEL_ROOT}
    PATH_SUFFIXES ${INTEL_LIB_SUBDIR})
  list(APPEND MKL_REQUIRED MKL_INTEL_THREAD MKL_INTEL_OPENMP)
endif()

if(MPI)
  find_library(MKL_SCALAPACK_LIB
    NAMES mkl_scalapack_lp64
    PATHS ${MKL_ROOT}
    PATH_SUFFIXES ${MKL_LIB_SUBDIR})
  find_library(MKL_BLACS_INTELMPI_LIB
    NAMES mkl_blacs_intelmpi_lp64
    PATHS ${MKL_ROOT}
    PATH_SUFFIXES ${MKL_LIB_SUBDIR})
  list(APPEND MKL_REQUIRED MKL_SCALAPACK MKL_BLACS_INTELMPI)
endif(MPI)

set(MKL_FOUND TRUE)
foreach(NAME ${MKL_REQUIRED})
  if(${NAME}_LIB)
    message(STATUS "Found ${NAME}_LIB: ${${NAME}_LIB}")
    list(APPEND MATH_LIBS ${${NAME}_LIB})
  else(${NAME}_LIB)
    message(STATUS "Could not find ${NAME}_LIB") 
    set(MATH_LIBS "")
    set(MKL_FOUND FALSE)
  endif(${NAME}_LIB)
endforeach(NAME)

if(NOT MKL_FOUND)
  message(FATAL_ERROR 
    "Could not find MKL libraries. Please manually specify MATH_LIBS")
endif( )

if(USE_THREADS)
  if("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    list(APPEND MATH_LIBS m ${OpenMP_C_FLAGS})
  elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    list(APPEND MATH_LIBS ${OpenMP_C_FLAGS})
  endif()
  list(APPEND MATH_LIBS ${CMAKE_THREAD_LIBS_INIT})
endif()
