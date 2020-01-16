# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib -lpacklib -lc -lm
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

F2dp_CJ_objs = F2dp_CJ.o f2glob.o ineft.o F2NMC.o CJgrids.o
F2dp_CJ_srcs = F2dp_CJ.f f2glob.f ineft.f F2NMC.f /home/hanjie/work/MARATHON/CJ-master/downloads/CJgrids.f

########################################

none: F2dp_CJ

all:  F2dp_CJ

F2dp_CJ.o: F2dp_CJ.f
		 $(F77) $(FFLAGS) -c $< -o $@

f2glob.o: f2glob.f
	$(F77) $(FFLAGS) -c $< -o $@

ineft.o: ineft.f
	$(F77) $(FFLAGS) -c $< -o $@

F2NMC.o: F2NMC.f
	$(F77) $(FFLAGS) -c $< -o $@

CJgrids.o: /home/hanjie/work/MARATHON/CJ-master/downloads/CJgrids.f
	$(F77) $(FFLAGS) -c $< -o $@

F2dp_CJ: $(F2dp_CJ_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(F2dp_CJ_objs) $(OTHERLIBS)


clean:
	rm -f *.o F2dp_CJ
