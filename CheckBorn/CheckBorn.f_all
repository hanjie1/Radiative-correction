        program CheckBorn
        implicit none

        real*8 e0,ep,theta,Mp,Z,A,thr
        real*8 F1,F2,FL,sigt,sigl,sig_dis,sig_qe,sigma,eps
        real*8 Q2,Wsq,f1qe,f2qe,flqe
        real*8 f1d,f2d,fld,f1dqe,f2dqe,fldqe
        real*8 x,nu,kappa,flux
        real*8 PI,pi2,alpha,cs,sn,tn
        integer         wfn,opt,ii,i
        character*80    filename,infile,outfile
        CHARACTER*72    COMMENT(3)
        logical   doqe,dfirst
        real*8  w1p,w2p,w1d,w2d,w1,w2,sigmott,w1n,w2n
        real*8 emc_func_slac
        external emc_func_slac
        real*8  emctmp
        real*8  sigtn,sigln,sigtp,siglp,flp,fln
        real*8  f1p,f2p,f1n,f2n


        Z=1.0
        A=1.0
        Mp=0.93827231
        PI=3.1415927
        pi2=pi*pi
        alpha = 1./137.036
        wfn=2

        write(6,*) 'Enter the input file name (in INP directory)'
        read(5,*) filename
        infile='INP/'//trim(filename)//'.inp'
        write(6,*) infile

        OPEN(UNIT=7,FILE=infile)
        READ(7,'(A72,/,A72,/,A72,/,/,/)') (COMMENT(i),i=1,3)

        outfile='OUT/f1f217/'//trim(filename)//'.out'
        open(unit=66,file=outfile)
        write(66,*) '***x   Q2   Theta   Eprime   Sig_Born   Sig_Born_In    Sig_Born_QE'

        outfile='OUT/INEFT/'//trim(filename)//'.out'
        open(unit=77,file=outfile)
        write(77,*) '***x   Q2   Theta   Eprime   Sig_Born   Sig_Born_In    Sig_Born_QE'

        dfirst = .true.
        call SQESUB(1.0,1.0,wfn,f2dqe,f1dqe,fLdqe,dfirst)

        do 99 ii=1,100000
           e0=0
           ep=0
           theta=0
           READ(7,'(F6.3,1x,F6.4,1x,F7.4)',END=100) E0,EP,theta
           IF (E0.LE.0.) goto 100
 
           thr = theta*PI/180.0
           cs = cos(thr/2.)
           sn = sin(thr/2.)
           tn = tan(thr/2.)

           Q2 = 4.*e0*ep*sn**2
           nu=e0-ep
           WSQ = -Q2 + Mp**2 + 2.0*Mp*nu
           x = Q2/2/Mp/nu
           eps = 1.0/(1+2.*(1+Q2/(4*Mp**2*x**2))*tn**2)
           kappa = abs((wsq-Mp**2))/2./Mp
           flux = alpha*kappa/(2.*pi2*Q2)*ep/e0/(1.-eps)

           sig_dis=0.0
           sig_qe=0.0
           if((Z.eq.1.).and.(A.eq.2.)) then
               doqe = .false.
               wfn=2
               dfirst = .false.
               call SQESUB(WSQ,Q2,wfn,f2dqe,f1dqe,fLdqe,dfirst)
               sigt = 0.3894e3*f1dqe*pi2*alpha*8.0/abs(wsq-Mp**2)
               sigl = 0.3894e3*fLdqe*pi2*alpha*8.0/abs(wsq-Mp**2)/2.*abs(wsq-Mp**2+q2)/q2
               sig_qe = sigt+eps*sigl
               sig_qe=flux*sig_qe
              
               call RESCSD(WSQ,Q2,eps,doqe,f1d,f2d,fLd,wfn,sig_dis)
               sig_dis=flux*sig_dis

           elseif((Z.eq.1.).and.(A.eq.1.)) then
               call RESCSP(wsq,q2,sigtp,siglp)
               sig_dis = flux*(sigtp+eps*siglp)

               opt=1
               call SFCROSS(WSQ,Q2,A,Z,opt,sigt,sigl,f1qe,f2qe,fLqe)
               sig_qe = flux*(sigt+eps*sigl)
           else
               opt=3
               call SFCROSS(WSQ,Q2,A,Z,opt,sigt,sigl,f1,f2,fL)
               sig_dis = flux*(sigt+eps*sigl)

               opt=1
               call SFCROSS(WSQ,Q2,A,Z,opt,sigt,sigl,f1qe,f2qe,fLqe)
               sig_qe = flux*(sigt+eps*sigl)
           endif

           sigma=(sig_qe+sig_dis)*1000
           write(66,*) x,Q2,theta,Ep,sigma,sig_dis*1000,sig_qe*1000


           call ineft(Q2,sqrt(wsq),W1p,W2p,dble(1.0))
           call ineft(Q2,sqrt(wsq),W1D,W2D,dble(2.0))

           W1n=2.0*W1D-W1p
           W2n=2.0*W2D-W2p

           W1=Z*W1p+(A-Z)*W1n
           W2=Z*W2p+(A-Z)*W2n
           sigmott=(19732.0/(2.0*137.0388*e0*sn**2))**2*cs**2/1.d6
           sig_dis = 1d3*sigmott*(W2+2.0*W1*tn**2)
           emctmp=emc_func_slac(x, A)
           sig_dis = sig_dis*emctmp
           sigma=sig_qe*1000+sig_dis
           write(77,*) x,Q2,theta,Ep,sigma,sig_dis,sig_qe*1000
c           write(6,*) 'INEFT:',w1p*mp,w2p*nu,w1n*mp,w2n*nu,emctmp

     
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

