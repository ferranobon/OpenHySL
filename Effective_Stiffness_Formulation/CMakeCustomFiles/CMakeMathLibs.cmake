if (USE_MKL)
  include("CMakeCustomFiles/CMakeMKL.cmake")
elseif(USE_ATLAS)
  include("CMakeCustomFiles/CMakeATLAS.cmake")
endif()

if( NOT (USE_MKL) AND
    NOT (USE_ATLAS) AND
    NOT (GOTO_BLAS) )
  if(NOT USE_MPI)
    set(REFERENCE_REQUIRED LAPACK BLAS)
  else(USE_MPI)
    set(REFERENCE_REQUIRED BLACSCINIT BLACS SCALAPACK LAPACK BLAS)
  endif( )
  # Look for BLAS/LAPACK libraries

  if(USE_MPI)
    find_library(BLACSCINIT_LIB NAMES blacsCinit-openmpi)
    find_library(BLACS_LIB NAMES blacs-openmpi)
    find_library(SCALAPACK_LIB NAMES scalapack-openmpi)
  endif()

  find_library(BLAS_LIB NAMES blas)
  find_library(LAPACK_LIB NAMES lapack reflapack) 
      
  set(REFERENCE_FOUND TRUE)
  foreach(NAME ${REFERENCE_REQUIRED})
    if(${NAME}_LIB)
      message(STATUS "Found ${NAME}_LIB: ${${NAME}_LIB}")
      list(APPEND MATH_LIBS ${${NAME}_LIB})
    else(${NAME}_LIB)
      message(STATUS "Could not find ${NAME}_LIB")
      set(MATH_LIBS "")
      set(REFERENCE_FOUND FALSE)
    endif(${NAME}_LIB)
  endforeach(NAME)
  if(REFERENCE_FOUND)
    message(STATUS "WARNING: Using reference BLAS/LAPACK.")
  else(REFERENCE_FOUND) 
    message(FATAL_ERROR 
      "Could not find BLAS/LAPACK libraries. Please manually specify MATH_LIBS")
  endif( )
endif()