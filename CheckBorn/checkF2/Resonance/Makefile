# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib -lpacklib -lc -lm
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

CheckRes_objs = CheckRes.o F2NMC.o ineft.o R1998.o
CheckRes_srcs = CheckRes.f F2NMC.f ineft.f R1998.f

########################################

none: CheckRes

all:  CheckRes

CheckRes.o: CheckRes.f
		 $(F77) $(FFLAGS) -c $< -o $@

F2NMC.o: F2NMC.f
	$(F77) $(FFLAGS) -c $< -o $@

ineft.o: ineft.f
	$(F77) $(FFLAGS) -c $< -o $@

R1998.o: R1998.f
	$(F77) $(FFLAGS) -c $< -o $@

CheckRes: $(CheckRes_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(CheckRes_objs) $(OTHERLIBS)


clean:
	rm -f *.o CheckRes
