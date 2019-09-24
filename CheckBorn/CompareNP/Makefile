# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib -lpacklib -lc -lm
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

EMCRD_CJ_objs = EMCRD_CJ.o CJgrids.o
EMCRD_CJ_srcs = EMCRD_CJ.f /home/hanjie/work/MARATHON/CJ-master/downloads/CJgrids.f

########################################

none: EMCRD_CJ

all:  EMCRD_CJ

EMCRD_CJ.o: EMCRD_CJ.f
		 $(F77) $(FFLAGS) -c $< -o $@

CJgrids.o: /home/hanjie/work/MARATHON/CJ-master/downloads/CJgrids.f
	$(F77) $(FFLAGS) -c $< -o $@

EMCRD_CJ: $(EMCRD_CJ_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(EMCRD_CJ_objs) $(OTHERLIBS)


clean:
	rm -f *.o EMCRD_CJ
