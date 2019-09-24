# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib -lpacklib -lc -lm
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

F2_Whitlow_objs = F2_Whitlow.o f2glob.o 
F2_Whitlow_srcs = F2_Whitlow.f f2glob.f

########################################

none: F2_Whitlow

all:  F2_Whitlow

F2_Whitlow.o: F2_Whitlow.f
		 $(F77) $(FFLAGS) -c $< -o $@

f2glob.o: f2glob.f
	$(F77) $(FFLAGS) -c $< -o $@

F2_Whitlow: $(F2_Whitlow_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(F2_Whitlow_objs) $(OTHERLIBS)


clean:
	rm -f *.o F2_Whitlow
