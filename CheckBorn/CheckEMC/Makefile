# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib -lpacklib -lc -lm
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

EMC_SLAC_objs = EMC_SLAC.o ineft.o f2glob.o CJgrids.o
EMC_SLAC_srcs = EMC_SLAC.f ineft.f f2glob.f \
		/home/hanjie/work/MARATHON/CJ-master/downloads/CJgrids.f

########################################

none: EMC_SLAC

all:  EMC_SLAC

EMC_SLAC.o: EMC_SLAC.f
		 $(F77) $(FFLAGS) -c $< -o $@

ineft.o: ineft.f
	$(F77) $(FFLAGS) -c $< -o $@

f2glob.o: f2glob.f
	$(F77) $(FFLAGS) -c $< -o $@

CJgrids.o: /home/hanjie/work/MARATHON/CJ-master/downloads/CJgrids.f
	$(F77) $(FFLAGS) -c $< -o $@

EMC_SLAC: $(EMC_SLAC_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(EMC_SLAC_objs) $(OTHERLIBS)


clean:
	rm -f *.o EMC_SLAC
