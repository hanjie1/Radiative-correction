# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib -lpacklib -lc -lm
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

CheckD2_objs = CheckD2.o f1f217b.o f2allm.o f2glob.o r1998.o vcoul.o ineft.o rescsp-120818.o
CheckD2_srcs = CheckD2.f f1f217b.f f2allm.f f2glob.f r1998.f vcoul.f ineft.f rescsp-120818.f

########################################

none: CheckD2

all:  CheckD2

CheckD2.o: CheckD2.f
		 $(F77) $(FFLAGS) -c $< -o $@

f2allm.o: f2allm.f
	$(F77) $(FFLAGS) -c $< -o $@

rescsp-120818.o: rescsp-120818.f
	$(F77) $(FFLAGS) -c $< -o $@

f2glob.o: f2glob.f
	$(F77) $(FFLAGS) -c $< -o $@

r1998.o: r1998.f
	$(F77) $(FFLAGS) -c $< -o $@

vcoul.o: vcoul.f
	$(F77) $(FFLAGS) -c $< -o $@

f1f217b.o: f1f217b.f
	$(F77) $(FFLAGS) -c $< -o $@

ineft.o: ineft.f
	$(F77) $(FFLAGS) -c $< -o $@

CheckD2: $(CheckD2_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(CheckD2_objs) $(OTHERLIBS)


clean:
	rm -f *.o CheckD2
