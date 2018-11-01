      PROGRAM CSCOMP
      IMPLICIT none

      real*8 e(20000),ep(20000),dep,th(20000),eps(20000),nu,x,q2(20000)
      real*8 w2(20000),mp,mp2,sin2,tan2,kappa,sigm,cs(20000)
      real*8 cserr(20000),sigtp,siglp,sigtn,sigln,sigt,sigl,f1d,f2d,fLd
      real*8 f1dqe,f2dqe,fLdqe,alpha,pi,pi2,gamma,flux,t1,t2,t3
      real*8 scale(12),rat,erat,q2c,dq2
      integer i,j,k,ntot,run(20000),set(20000),purge,wfn,dset
      logical doqe/.true./, first/.true./, thend/.false./
      character*40 filename
      data scale / 1.005,0.9961,0.9993,0.9936,1.0192,1.000,
     &             0.9967,1.000,1.000,1.000,1.0094,1.0064 /

      dset = 3
      wfn = 2

      
      mp = (0.938272+0.939565)/2.0d0
      mp2 = mp*mp
  
      alpha = 1./137.036
      pi = 3.141593
      pi2 = pi*pi

c      write(6,*) "here1"
      
c      first = .true.
      call SQESUB(1.0,1.0,wfn,t1,t2,t3,first)
      first = .false.

c      write(6,*) "here2"   

      ntot = 0
      i=0

c      READ(5,*) q2c,dq2,dset

c     OPEN(UNIT=15,FILE='cs_d2_tot_it2_nogrid.dat')
      
      OPEN(UNIT=15,FILE='deuterium.table')
     
      DO WHILE(.NOT.THEND)
        I = I+1
C       READ(15,*,END=3000) E(I),EP(I),TH(I),Q2(I),W2(I),NU(I),EPS(I),
C     &         GAMMA(I),X(I),CS(I),CSERR(I) 
        READ(15,*,END=2000) e(i),ep(i),th(i),w2(i),cs(i),cserr(i),t1

        set(i) = 2
        cs(i) = scale(set(i))*cs(i)
        cserr(i) = scale(set(i))*cserr(i)

        ntot = i

        doqe = .false.
       if(set(i).eq.4.OR.set(i).eq.7.OR.set(i).eq.8) then
              doqe = .true.
       endif

        if(set(i).eq.4.OR.set(i).eq.7.OR.set(i).eq.8.OR.set(i).eq.3) 
     &         doqe = .true.

        nu = e(i)-ep(i)
        sin2 = dsin(pi*th(i)/180.0/2.0)
        sin2 = sin2*sin2
        tan2 = sin2/(1.0-sin2)
        q2(i) = 4.*e(i)*ep(i)*sin2

c        write(6,*) e(i),ep(i),cs(i)

        kappa = abs((w2(i)-mp2))/2./mp
        eps = 1./(1. + 2.*(nu*nu+q2)/q2(i)*tan2)
        flux = alpha*kappa/(2.*pi*pi*q2(i))*ep(i)/e(i)/(1.-eps(i)) 
        x = q2(i)/(w2(i)-mp2+q2(i))
        gamma = sqrt(1.+4.*mp2*x*x/q2(i))   

c        write(6,*) "here1"

c        if(set(i).EQ.4) then

        call RESCSD(w2(i),q2(i),eps(i),doqe,f1d,f2d,fLd,wfn,sigm)
        sigm = flux*sigm/2.0
        rat = cs(i)/sigm
        erat = cserr(i)/cs(i)*rat
        write(6,3000) set(i),e(i),ep(i),th(i),w2(i),q2(i),eps(i),cs(i),
     &                 sigm,rat,erat 
      enddo

 2000 thend = .true.


 3000 format(1i5,6f8.3,4f13.4)

c 2000 format(1i5,6f10.4,5e12.3)


      end






