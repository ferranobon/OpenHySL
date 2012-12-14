# find the libraries required by ATLAS
set(ATLAS_ROOT /usr/lib/altas-base)

SET(ATLAS_REQUIRED LAPACK_ATLAS)

if(MPI)
  find_library(BLACSCINIT_LIB NAMES blacsCinit-openmpi)
  find_library(BLACS_LIB NAMES blacs-openmpi)
  find_library(SCALAPACK_LIB NAMES scalapack-openmpi)
  list(APPEND ATLAS_REQUIRED BLACSCINIT SCALAPACK BLACS)
endif()

# Try an atlas built lapack
find_library(LAPACK_ATLAS_LIB
  NAMES lapack
  PATHS ${ATLAS_ROOT} ${ATLAS_ROOT}/atlas)

if(USE_THREADS)
  find_package(Threads)
  find_library(ATLAS_BLAS_THREADS_LIB
    NAMES ptf77blas
    PATHS ${ATLAS_ROOT})
  list(APPEND ATLAS_REQUIRED ATLAS_BLAS_THREADS)
else()
  find_library(ATLAS_BLAS_LIB
    NAMES f77blas
    PATHS ${ATLAS_ROOT})
  list(APPEND ATLAS_REQUIRED ATLAS_BLAS)
endif()

find_library(ATLAS_LIB
  NAMES atlas
  PATHS ${ATLAS_ROOT})

list(APPEND ATLAS_REQUIRED ATLAS)

set(ATLAS_FOUND TRUE)
foreach(NAME ${ATLAS_REQUIRED})
  if(${NAME}_LIB)
    message(STATUS "Found ${NAME}_LIB: ${${NAME}_LIB}")
    list(APPEND MATH_LIBS ${${NAME}_LIB})
  else(${NAME}_LIB)
    message(STATUS "Could not find ${NAME}_LIB") 
    set(MATH_LIBS "")
    set(ATLAS_FOUND FALSE)
  endif(${NAME}_LIB)
endforeach(NAME)

if(NOT ATLAS_FOUND)
  message(FATAL_ERROR 
    "Could not find ATLAS libraries. Please manually specify MATH_LIBS")
endif( )

if(USE_THREADS)
  list(APPEND MATH_LIBS ${CMAKE_THREAD_LIBS_INIT})
endif()