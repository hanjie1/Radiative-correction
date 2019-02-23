# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib -lpacklib -lc -lm
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

EMCHE3_objs = EMCHE3.o ineft.o 
EMCHE3_srcs = EMCHE3.f ineft.f

########################################

none: EMCHE3

all:  EMCHE3

EMCHE3.o: EMCHE3.f
		 $(F77) $(FFLAGS) -c $< -o $@

ineft.o: ineft.f
	$(F77) $(FFLAGS) -c $< -o $@

EMCHE3: $(EMCHE3_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(EMCHE3_objs) $(OTHERLIBS)


clean:
	rm -f *.o EMCHE3
