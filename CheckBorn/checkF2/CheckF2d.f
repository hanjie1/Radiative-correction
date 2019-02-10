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
        real*8 NMCF2d
        external NMCF2d
        real*8  emctmp
        real*8  sigtn,sigln,sigtp,siglp,flp,fln
        real*8  f1p,f2p,f1n,f2n
        real*8  f2data,f2_staterr,f2_syserr


        Mp=0.93827231
        PI=3.1415927
        pi2=pi*pi
        alpha = 1./137.036
        wfn=2

        infile='../data/F2d_Whitlow.dat'
        write(6,*) infile

        OPEN(UNIT=7,FILE=infile)
        READ(7,'(A72)') COMMENT

        outfile='OUT/F2d_NMC.out'
        open(unit=66,file=outfile)
        write(66,*) 'x    Q2    WSQ   F2d_data   F2d_stat_err',
     >      '  F2d_sys_err    F2d_NMC    F2d_ineft'
     
        dfirst = .true.
        call SQESUB(1.0,1.0,wfn,f2dqe,f1dqe,fLdqe,dfirst)

        do 99 ii=1,100000
           x=0.0
           Q2=0.0
           READ(7,'(F4.3,1x,F6.3,1x,F7.5,1x,F6.3,F5.3)',END=100) x,Q2,F2data,
     >     f2_staterr,f2_syserr
c           READ(7,'(F6.4,F8.2,F8.4,F9.4,F9.4)',END=100) x,Q2,F2data,
c     >     f2_syserr,f2_staterr
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
c           call rescsd(wsq,q2,eps,doqe,f1d,f2d,fLd,wfn,sigma)

           call ineft(Q2,sqrt(wsq),W1D,W2D,dble(2.0))
           f2d_ineft=w2d*nu
 
           F2D=NMCF2d(x,Q2)

           write(66,'(3f10.5,3F13.5,2F17.5)')x,q2,wsq,f2data,f2_staterr,f2_syserr,f2d,f2d_ineft
     
99      continue
100     continue

        end

c-------------------------------------------------------------------------------------------
        real*8 function NMCF2d(x,Q2)
        real*8 x,Q2,W,beta
        external beta
        real*8 F2DIS,F2RES,F2BG
        real*8 GQ2,S,sbar,xw,eps
        real*8 a /4.177/,Wthr /1.03/,Mp /0.938272/,Q02 /2.0/,lamda /0.2/
        real*8 ma2 /0.351/,mb2 /1.512/,Mdelta /1.232/,gammac /0.0728/
        real*8 b/0.5/,c/0.05/,Mpi/0.138/
        real*8 eta(4)
        real*8 aa(6),bb(4)

        DATA aa(1) /0.75966/, aa(2) /3.52/, aa(3) /0.83691/,
     *       aa(4) /12.876/, aa(5) /0.89456/, aa(6) /0.16452/
        data bb(1) /-0.18202/, bb(2) /0.46256/, bb(3) /0.97906/, bb(4) /-2.9558/

        sbar=log(log((Q2+ma2)/lamda**2)/log((Q02+ma2)/lamda**2))
        eta(1)=aa(1)+bb(1)*sbar
        eta(2)=aa(2)+bb(2)*sbar
        eta(3)=aa(3)+bb(3)*sbar
        eta(4)=aa(4)+bb(4)*sbar

        W=sqrt(Mp*Mp+Q2*(1/x-1))
        GQ2=1.0/(1.0+Q2/0.71)**2
        S=1-exp(-a*(W-Wthr))
        xw=(Q2+ma2)/(Q2/x+mb2)
        F2DIS=(5.0/18.0*3.0/beta(eta(1),eta(2)+1.0)*xw**eta(1)*(1-xw)**eta(2)
     >       +1.0/3.0*eta(3)*(1-xw)**eta(4))*S

        F2RES=aa(5)**2*sqrt(GQ2**3)*exp(-(W-Mdelta)**2/gammac**2)

        eps=sqrt(((W+c)**2+Mp**2-Mpi**2)**2/(4*(W+c)**2)-Mp**2)

        F2BG=aa(6)**2*sqrt(GQ2)*eps*exp(-b*(W-Wthr)**2)

        NMCF2d=(1-GQ2**2)*(F2DIS+F2RES+F2BG)
        return
        end

c----------------------------------------------------------------------------------------------
        real*8 function beta(z,w)
        real*8 z,w,gammln
        external gammln
        beta=exp(gammln(z)+gammln(w)-gammln(z+w))
        return
        end

        real*8 function gammln(xx)
        real*8 xx,x,tmp
        dimension cof(6)
        data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     *    -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
        data half,one,fpf/0.5d0,1.0d0,5.5d0/
        x=xx-one
        tmp=x+fpf
        tmp=(x+half)*log(tmp)-tmp
        ser=one
        do 11 j=1,6
           x=x+one
           ser=ser+cof(j)/x
11      continue
        gammln=tmp+log(stp*ser)
        return
        end

