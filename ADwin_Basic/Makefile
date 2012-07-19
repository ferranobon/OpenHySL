BINFILES = Substepping_Process Table_Control
BUILDDIR = build

default: $(BINFILES)

Substepping_Process: Substepping_Process.bas
	adbasic /M /SPII /P11 /PN2 /PL /O2 $< /A$(BUILDDIR)/Substepping_Process
Table_Control: Table_Control.bas
	adbasic /M /SPII /P11 /O2 $< /A$(BUILDDIR)/Table_Control

.PHONY: clean
clean:
	rm -r build