 gfortran -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132 -o cscompq2eps cscompq2eps.f f1f217b.f f2allm.f f2glob.f r1998.f vcoul.f -L/apps/cernlib/x86_64_rhel7/2005/lib -lmathlib -lpacklib -lc -lm
# gfortran -o cscompq2eps cscompq2eps.f f1f217b.f f2allm.f f2glob.f r1998.f vcoul.f -L/user/local/cern/2005/lib -lpacklib -lmathlib -lc -lm


# gfortran -o cscompq2eps cscompq2eps.f sfcross.f csfit.f  f1f216all.f f2allm.f f2glob.f r1998.f vcoul.f -L/user/local/cern/2005/lib -lpacklib -lmathlib -lc -lm

