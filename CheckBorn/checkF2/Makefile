# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib -lpacklib -lc -lm
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

CompareF2d_objs = CompareF2d.o f2glob.o ineft.o F2NMC.o CJgrids.o
CompareF2d_srcs = CompareF2d.f f2glob.f ineft.f F2NMC.f /home/hanjie/work/MARATHON/CJ-master/downloads/CJgrids.f

########################################

none: CompareF2d

all:  CompareF2d

CompareF2d.o: CompareF2d.f
		 $(F77) $(FFLAGS) -c $< -o $@

f2glob.o: f2glob.f
	$(F77) $(FFLAGS) -c $< -o $@

F2NMC.o: F2NMC.f
	$(F77) $(FFLAGS) -c $< -o $@

ineft.o: ineft.f
	$(F77) $(FFLAGS) -c $< -o $@

CJgrids.o: /home/hanjie/work/MARATHON/CJ-master/downloads/CJgrids.f
	$(F77) $(FFLAGS) -c $< -o $@

CompareF2d: $(CompareF2d_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(CompareF2d_objs) $(OTHERLIBS)


clean:
	rm -f *.o CompareF2d
