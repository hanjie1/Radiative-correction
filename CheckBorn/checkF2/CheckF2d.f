        program CheckF2d
        implicit none

        real*8 e0,ep,theta,Mp,Z,A,thr
        real*8 F1,F2,FL,sigt,sigl,sig_dis,sig_qe,sigma,eps
        real*8 Q2,Wsq,f1qe,f2qe,flqe
        real*8 f1d,f2d,fld,f1dqe,f2dqe,fldqe
        real*8 x,nu,kappa,flux
        real*8 PI,pi2,alpha,cs,sn,tn
        integer         wfn,opt,ii,i
        character*80    filename,infile,outfile
        CHARACTER*72    COMMENT
        logical   doqe,dfirst
        real*8  w1d,w2d,f1d_ineft,f2d_ineft
        real*8 emc_func_slac
        external emc_func_slac
        real*8  emctmp
        real*8  sigtn,sigln,sigtp,siglp,flp,fln
        real*8  f1p,f2p,f1n,f2n
        real*8  f2data,f2_staterr,f2_syserr


        Mp=0.93827231
        PI=3.1415927
        pi2=pi*pi
        alpha = 1./137.036
        wfn=2

        infile='../data/Marathon.dat'
        write(6,*) infile

        OPEN(UNIT=7,FILE=infile)
        READ(7,'(A72)') COMMENT

        outfile='OUT/F2d_marathon.out'
        open(unit=66,file=outfile)
        write(66,*) 'x    Q2    WSQ   F2d_data   F2d_stat_err',
     >      '  F2d_sys_err    F2d_f1f217    F2d_ineft'
     
        dfirst = .true.
        call SQESUB(1.0,1.0,wfn,f2dqe,f1dqe,fLdqe,dfirst)

        do 99 ii=1,100000
           x=0.0
           Q2=0.0
           READ(7,'(F4.3,1x,F6.3,1x,F7.5,1x,F6.3,F5.3)',END=100) x,Q2,F2data,
     >     f2_staterr,f2_syserr
c           READ(7,'(F9.4,F13.3,F14.4,F14.4,F11.4)',END=100) x,Q2,F2data,
c     >     f2_staterr,f2_syserr
           IF (x.LE.0.) goto 99
           IF (Q2.GE.14.) goto 99
 
           nu=Q2/(2.0*MP*x)
           WSQ = -Q2 + Mp**2 + 2.0*Mp*nu
           kappa = abs((wsq-Mp**2))/2./Mp

           sig_dis=0.0
           sig_qe=0.0

           eps=1.0
           doqe = .true.
           call rescsd(wsq,q2,eps,doqe,f1d,f2d,fLd,wfn,sigma)

           call ineft(Q2,sqrt(wsq),W1D,W2D,dble(2.0))
           f2d_ineft=w2d*nu


           write(66,'(3f10.5,3F13.5,2F17.5)')x,q2,wsq,f2data,f2_staterr,f2_syserr,f2d,f2d_ineft
     
99      continue
100     continue

        end

c-------------------------------------------------------------------------------------------
        real*8 function emc_func_slac(x,A)
        real*8 x,A,atemp
        real*8 alpha,C

        atemp = A
!       if(A.eq.4.or.A.eq.3) then  ! emc effect is more like C for these
!       2...
!          atemp = 12
!       endif

        alpha = -0.070+2.189*x - 24.667*x**2 + 145.291*x**3
     >        -497.237*x**4 + 1013.129*x**5 - 1208.393*x**6
     >        +775.767*x**7 - 205.872*x**8

        C = exp( 0.017 + 0.018*log(x) + 0.005*log(x)**2)
        
        emc_func_slac = C*atemp**alpha
        return
        end

