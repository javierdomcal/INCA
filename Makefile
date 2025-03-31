# Makefile for RODA program
ODIR = obj
FC = ifort
FFLAGS = -O3 -g

# Libraries
LIB = -llapack -lblas
LIBPATH=-L/dipc/markelgtx

# Define all object files
OBJS = $(ODIR)/mod5.o \
       $(ODIR)/cube.o \
       $(ODIR)/properties.o \
       $(ODIR)/density_properties.o \
       $(ODIR)/pair_density.o \
       $(ODIR)/correlation_indicators.o \
       $(ODIR)/functions.o \
       $(ODIR)/derive.o \
       $(ODIR)/lebedev.o \
       $(ODIR)/intrastuff.o \
       $(ODIR)/gridpoints.o \
       $(ODIR)/integrals.o \
       $(ODIR)/read_files.o \
       $(ODIR)/c1hole.o \
       $(ODIR)/gen_cube.o \
       $(ODIR)/intracule.o \
       $(ODIR)/dm2phf.o \
       $(ODIR)/scan_utils.o \
       $(ODIR)/on_top.o \
       $(ODIR)/inca_main.o

# Create object directory

# Default target - make this depend on the obj directory
roda.x: $(ODIR) $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $@ -qmkl -O3
	rm -f *.mod

# Compilation rules with explicit dependencies
$(ODIR)/mod5.o: mod5.f90 | $(ODIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/cube.o: cube.f90 $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/properties.o: properties.f90 $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@


$(ODIR)/density_properties.o: density_properties.f90 $(ODIR)/properties.o $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/pair_density.o: pair_density.f90 $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/correlation_indicators.o: correlation_indicators.f90 $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/functions.o: functions.f90 $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/derive.o: derive.f90 $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/lebedev.o: lebedev.f | $(ODIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/intrastuff.o: intrastuff.f90 $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/gridpoints.o: gridpoints.f90 $(ODIR)/mod5.o $(ODIR)/lebedev.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/integrals.o: integrals.f90 $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/read_files.o: read_files.f90 $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/c1hole.o: c1hole.f90 $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/gen_cube.o: gen_cube.f90 $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/intracule.o: intracule.f90 $(ODIR)/mod5.o $(ODIR)/intrastuff.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/dm2phf.o: dm2phf.f90 $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/scan_utils.o: scan_utils.f90 $(ODIR)/mod5.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)/on_top.o: on_top.f90 $(ODIR)/mod5.o $(ODIR)/scan_utils.o
	$(FC) $(FFLAGS) -c $< -o $@


$(ODIR)/inca_main.o: inca_main.f90 $(ODIR)/mod5.o $(ODIR)/cube.o
	$(FC) $(FFLAGS) -c $< -o $@

# Clean up
clean:
	rm -rf $(ODIR) *.mod roda.x

.PHONY: clean
