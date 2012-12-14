CC     = gcc
CFLAGS = -m64 -std=c99 -pedantic -Wall -W -Wextra -Wmissing-prototypes -Wstrict-prototypes -Wbad-function-cast -Wundef -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -Wstrict-overflow=5 -fshort-enums -fno-common -g3
LDFLAGS =

OBJDIR := objdir
SRC := Full_To_Sparse.c

default: Full_To_Sparse

Full_To_Sparse: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: $(CC) Compiler'
	$(CC) $(CFLAGS) $(SRC) -o Full_To_Sparse $(OBJS) $(LDFLAGS)
	@echo 'Finished building target: $@'
	@echo ' '

.PHONY: clean
clean:
	@echo 'Removing objects and executables'
	rm -r $(OBJDIR)
	@echo ' '