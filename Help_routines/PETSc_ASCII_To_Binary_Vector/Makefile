ALL: Ascii_to_binari

CFLAGS     =
FFLAGS     =
CPPFLAGS   =
FPPFLAGS   =
CLEANFILES = Ascii_to_binary


include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

Ascii_to_binari: Ascii_to_binary.o chkopts
	${CLINKER} -o Ascii_to_binary Ascii_to_binary.o ${PETSC_LIB}
	${RM} ex2.o