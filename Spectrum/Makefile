# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

spectrum_objs = spectrum.o ineft.o gsmearing.o
spectrum_srcs = spectrum.f ineft.f gsmearing.f

########################################

none: spectrum

all:  spectrum

spectrum.o: spectrum.f
		 $(F77) $(FFLAGS) -c $< -o $@

ineft.o: ineft.f
	$(F77) $(FFLAGS) -c $< -o $@


gsmearing.o: gsmearing.f
	$(F77) $(FFLAGS) -c $< -o $@


spectrum: $(spectrum_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(spectrum_objs) $(OTHERLIBS)


clean:
	rm -f *.o spectrum
