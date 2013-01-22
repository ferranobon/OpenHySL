CC         = mpicc.openmpi
CFLAGS     = -Wall
FFLAGS     =
CPPFLAGS   =
FPPFLAGS   =
CLEANFILES = Ascii_to_binary

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

OBJDIR := objdir
OBJS := $(addprefix $(OBJDIR)/, Conf_Parser.o ErrorHandling.o EndingStep.o Initiation.o MatrixVector.o ComputeU0.o Precalculations.o NSEP_Communication_Sync.o OpenFresco_Communication.o OPSocket.o Send_Receive_Data.o main.o)

all: Petsc

#Invocation
Petsc: $(OBJS) chkopts
	@echo 'Building target: $@'
	@echo 'Invoking: $(CC) Linker'
	mpicc.openmpi -L/opt/adwin/lib -o Petsc $(OBJS) $(LDFLAGS) ${PETSC_LIB}
	@echo 'Finished building target: $@'
	@echo ' '


$(OBJS): | $(OBJDIR)

$(OBJDIR):
	@echo 'Creating build directory' 
	mkdir $(OBJDIR)
	@echo ' '

$(OBJDIR)/%.o: %.c
	@echo 'Building file: $<'
	@echo 'Invoking: $(CC) Compiler'
	mpicc.openmpi $(CFLAGS) -I/usr/include/petsc -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Other Targets
clean2:
	@echo 'Removing objects and executables'
	rm -r $(OBJDIR) Petsc
	@echo ' '