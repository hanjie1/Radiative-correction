        program CheckD2
        implicit none

        real*8 e0,ep,theta,Mp,Z,A,thr,rc_factor
        real*8 F1,F2,FL,sigt,sigl,sig_dis,sig_qe,sigma,eps,sigmott
        real*8 Q2,Wsq,f1qe,f2qe,flqe
        real*8 f1d,f2d,fld,f1dqe,f2dqe,fldqe
        real*8 x,nu,kappa,flux
        real*8 PI,pi2,alpha,cs,sn,tn
        integer         wfn,opt,ii,i
        character*80    filename,infile,outfile
        CHARACTER*72    COMMENT
        logical   doqe,dfirst
        real*8  w1d,w2d,XS_ineft
        real*8 emc_func_slac
        external emc_func_slac
        real*8  emctmp
        real*8  sigtn,sigln,sigtp,siglp,flp,fln
        real*8  f1p,f2p,f1n,f2n
        real*8  XSdata,XS_staterr,XS_syserr
        integer tmp(2)

        Z=1.0
        A=2.0
        Mp=0.93827231
        PI=3.1415927
        pi2=pi*pi
        alpha = 1./137.036
        wfn=2

        infile='../data/XS_Whitlow.dat'
        write(6,*) infile

        OPEN(UNIT=7,FILE=infile)
        READ(7,'(A72)') COMMENT

        outfile='OUT/D2_Whitlow.out'
        open(unit=66,file=outfile)
        write(66,*) 'x    Q2    WSQ   XS_data   XS_stat_err',
     >      '  XS_sys_err    XS_f1f217    XS_ineft'
     
        dfirst = .true.
        call SQESUB(1.0,1.0,wfn,f2dqe,f1dqe,fLdqe,dfirst)

        do 99 ii=1,100000
           x=0.0
           Q2=0.0
           READ(7,'(I3,I2,F8.3,2F7.3,F5.3,F7.3,F5.3,2F7.3,E11.4,2F5.3)',END=100) tmp(1),tmp(2),
     >      E0,Ep,theta,x,Q2,eps,wsq,rc_factor,XSdata,XS_staterr,XS_syserr
c           READ(7,'(F9.4,F13.3,F14.4,F14.4,F11.4)',END=100) x,Q2,F2data,
c     >     f2_staterr,f2_syserr
           IF (x.LE.0.) goto 99
           IF (Q2.GE.14.) goto 99
 
           thr = theta*PI/180.0
           sn = sin(thr/2.)
           cs = cos(thr/2.)
           tn = tan(thr/2.)
 
           nu=E0-Ep
           kappa = abs((wsq-Mp**2))/2./Mp

           sig_dis=0.0
           sig_qe=0.0

           doqe = .true.
           call rescsd(wsq,q2,eps,doqe,f1d,f2d,fLd,wfn,sigma)
           flux = alpha*kappa/(2.*pi2*Q2)*ep/e0/(1.-eps)
           sigma=1000*sigma*flux/2.0

           call ineft(Q2,sqrt(wsq),W1D,W2D,dble(2.0))
           sigmott=(19732.0/(2.0*137.036*e0*sn**2))**2*cs**2/1.d6
           XS_ineft = 1d3*sigmott*(W2D+2.0*W1D*tn**2)
           XS_ineft = XS_ineft*emc_func_slac(x, A)

           write(66,'(3f10.5,3F16.5,2F17.5)')x,q2,wsq,XSdata,XS_staterr,XS_syserr,sigma,XS_ineft
     
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

