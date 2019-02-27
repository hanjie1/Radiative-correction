        program CheckEMC
        implicit none

c       Close EMC check
        real*8 e0,ep,theta,Mp,thr
        real*8 Q2,Wsq,f1qe,f2qe,flqe
        real*8 f1d,f2d,fld,f1dqe,f2dqe,fldqe
        real*8 x,nu,kappa,flux
        real*8 PI,pi2,alpha,cs,sn,tn
        integer         wfn,opt,ii,i
        character*80    filename,infile,outfile
        CHARACTER*72    COMMENT
        logical   doqe,dfirst
        real*8  w1d,w2d
        real*8  emctmp
        real*8  muA(10),epsA,F2A,lamdaA(10)
        real*8  lamdaD /1.015/
        real*8  tmp1,Q2A
        integer AA(10)

        DATA    muA /0.6156,0.5719,0.6097,0.5893,0.5773,0.5463,0.5424,0.5463,0.5124,0.5/
        DATA    lamdaA /1.04,1.079,1.045,1.063,1.074,1.104,1.108,1.104,1.14,1.154/
        DATA    AA /3,4,6,7,9,12,16,20,27,56/

        Mp=0.93827231
        PI=3.1415927
        pi2=pi*pi
        alpha = 1./137.036
        

        outfile='OUT/EMC_Close.out'
        open(unit=66,file=outfile)
        write(66,*) 'A   x   eps   F2A   F2A/F2D '
     
        x=0.6
        Q2=10
        nu=Q2/(2.0*MP*x)
        WSQ = -Q2 + Mp**2 + 2.0*Mp*nu
        call ineft(Q2,sqrt(wsq),W1D,W2D,dble(2.0))
        f2d=w2d*nu
        write(66,'(I4,4F13.5)') 2,x,1.0,f2d,1.0

        do 99 ii=1,10
           tmp1=2*log(Q2/0.25**2)/log(muA(ii)/0.25**2) 
           epsA=(lamdaA(ii)/lamdaD)**tmp1
           Q2A=epsA*Q2 
           nu=Q2A/(2.0*MP*x)
           WSQ = -Q2A + Mp**2 + 2.0*Mp*nu

           call ineft(Q2A,sqrt(wsq),W1D,W2D,dble(2.0))
           write(6,*) AA(ii),Q2A,W1D
           F2A=w2d*nu
 
           write(66,'(I4,4F13.5)') AA(ii),x,epsA,F2A,F2A/F2D 
     
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

