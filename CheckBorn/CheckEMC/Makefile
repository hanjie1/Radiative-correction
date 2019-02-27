# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib -lpacklib -lc -lm
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

CheckA3_objs = CheckA3.o ineft.o 
CheckA3_srcs = CheckA3.f ineft.f

########################################

none: CheckA3

all:  CheckA3

CheckA3.o: CheckA3.f
		 $(F77) $(FFLAGS) -c $< -o $@

ineft.o: ineft.f
	$(F77) $(FFLAGS) -c $< -o $@

CheckA3: $(CheckA3_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(CheckA3_objs) $(OTHERLIBS)


clean:
	rm -f *.o CheckA3
