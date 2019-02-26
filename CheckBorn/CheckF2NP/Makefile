# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib -lpacklib -lc -lm
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

CheckF2NP_objs = CheckF2NP.o CJgrids.o 
CheckF2NP_srcs = CheckF2NP.f /home/hanjie/work/MARATHON/CJ-master/downloads/CJgrids.f

########################################

none: CheckF2NP

all:  CheckF2NP

CheckF2NP.o: CheckF2NP.f
		 $(F77) $(FFLAGS) -c $< -o $@

CJgrids.o: /home/hanjie/work/MARATHON/CJ-master/downloads/CJgrids.f
	$(F77) $(FFLAGS) -c $< -o $@

CheckF2NP: $(CheckF2NP_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(CheckF2NP_objs) $(OTHERLIBS)


clean:
	rm -f *.o CheckF2NP
