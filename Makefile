#
#Makefile for RODA program, find a better name (;
#
ODIR=obj
FC=ifort
#FC=gfortran

#FFLAGS2= -Wall -g -msse4.2 -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant
FFLAGS= -o
#FFLAGS = -g
LIB= -llapack -lblas
LIBPATH=-L/dipc/markelgtx 
SRCF90=$(wildcard *.f90)
SRCF=$(wildcard *.f)

#OBJ=$(SRCF90:.f90=.o) $(SRCF:.f=.o)
OBJ=$(patsubst %.f90,$(ODIR)/%.o,$(SRCF90)) $(patsubst %.f,$(ODIR)/%.o,$(SRCF))

$(ODIR)/%.o: %.f90 
	$(FC) $(FFLAGS) $(ODIR)/mod5.o -c mod5.f90 -O3
	$(FC) $(FFLAGS) $@ -c $< -O3

$(ODIR)/%.o: %.f
	$(FC) $(FFLAGS) $@ -c $< -O3


roda.x: $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) $@ -qmkl -O3
	rm -f *.mod

clean:
	@rm -f *.mod $(ODIR)/*.o roda	

#Testing
test: 
	./test.sh >& test.log


