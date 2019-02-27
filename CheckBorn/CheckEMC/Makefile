# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib -lpacklib -lc -lm
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

CheckHE3_objs = CheckHE3.o ineft.o f2glob.o
CheckHE3_srcs = CheckHE3.f ineft.f f2glob.f

########################################

none: CheckHE3

all:  CheckHE3

CheckHE3.o: CheckHE3.f
		 $(F77) $(FFLAGS) -c $< -o $@

ineft.o: ineft.f
	$(F77) $(FFLAGS) -c $< -o $@

f2glob.o: f2glob.f
	$(F77) $(FFLAGS) -c $< -o $@

CheckHE3: $(CheckHE3_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(CheckHE3_objs) $(OTHERLIBS)


clean:
	rm -f *.o CheckHE3
