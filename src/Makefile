#!/bin/bash
FC	= mpif90
FFLAGS  = -O3
LFLAGS  = $(FFLAGS)
EXE	= DC_run
SRC     = modules.f90 \
	 main.f90 \
	 get_input_data.f90 \
	 DC.f90 \
	 PG.f90 \
	 SG.f90 \
	 acoustic_solver.f90 \
	 auxiliar.f90


# No need to edit below this line
.SUFFIXES: .f90 .o

OBJ = $(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ)

$(OBJ): $(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:
	rm -f *.o *.mod $(EXE) core
