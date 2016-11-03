FLAGS           =
FFLAGS           =
CPPFLAGS         =
FPPFLAGS         =
CLEANFILES       =

### This can be defined from outside:
# make COMP=gfortran
UNAME = ${shell uname}

COMP=gfortran
STATIC=no

F77C= $(COMP)
F90C= $(COMP)
F90OPTS=-O2
F77OPTS=-O2
LINKOPTS=-I.

MKDIR_P=mkdir -p
BIN=../bin
.PHONY: ${BIN}

all: ${BIN} IBr harm Na2_RWA3

${BIN}:
	${MKDIR_P} ${BIN}

harm: model_harm.o run_model.o
	$(F90C) $(F90OPTS) -o ${BIN}/run_harm.exe model_harm.o run_model.o \
	$(LINKOPTS)

IBr: model_IBr.o run_model.o
	$(F90C) $(F90OPTS) -o ${BIN}/run_IBr.exe model_IBr.o run_model.o \
	$(LINKOPTS)

Na2_RWA3: model_Na2_RWA_3pot.o run_model.o
	$(F90C) $(F90OPTS) -o ${BIN}/run_Na2_RWA3.exe model_Na2_RWA_3pot.o run_model.o \
	$(LINKOPTS)


clean:
	rm -f *.o *.mod

.SUFFIXES : .f90 .f

.f90.o:
	$(F90C) $(F90OPTS) -c $<

.f.o:
	$(F77C) $(F77OPTS) -c $<

