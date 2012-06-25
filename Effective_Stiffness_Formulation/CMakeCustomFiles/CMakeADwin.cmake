# find ADwin libraries

set(ADWIN_ROOT /opt/adwin)
set(ADWIN_LIB_SUBDIR /lib/)
set(ADWIN_INC_SUBDIR /include/)

find_library(ADWIN_CORE_LIB
  NAMES adwin
  PATHS ${ADWIN_ROOT}
  PATH_SUFFIXES ${ADWIN_LIB_SUBDIR}
)

list(APPEND ADWIN_REQUIRED ADWIN_CORE)

set(ADWIN_FOUND TRUE)
foreach(NAME ${ADWIN_REQUIRED})
  if(${NAME}_LIB)
    message(STATUS "Found ADwin libraries: ${ADWIN_CORE_LIB}")
    list(APPEND ADWIN_LIBS ${${NAME}_LIB})
  else(ADWIN_LIB)
    message(STATUS "Could not find ADWIN_LIB")
    set(ADWIN_LIBS "")
  endif()
endforeach()

include_directories(/opt/adwin/include)

if(NOT ADWIN_FOUND)
  message(FATAL_ERROR "Could not find ADwin libraries. Please manually specify ADWIN_LIBS")
endif()