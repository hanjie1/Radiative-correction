        program CheckRes
        implicit none

        real*8 Q2,xx,WSQ
        real*8 W1,W2,F2,F2_1
        integer ii
        real   x4,Qsq4,F2D4,F2Derr_lo4,F2Derr_hi4
        real*8 NMCF2d
        external NMCF2d


        open(unit=77,file='F2D_NMC.dat')

        do 10 ii=0,40
           xx=0.15+ii*0.02
           Q2=14.0*xx
           wsq=0.938**2+Q2*(1.0/xx-1)
           
           x4=xx
           Qsq4=Q2

c           call ineft(Q2,sqrt(wsq),w1,w2,dble(2.0))
           call F2NMC_new(2,x4,Qsq4,F2D4,F2Derr_lo4,F2Derr_hi4)
           F2=F2D4
           F2_1=NMCF2d(xx,Q2)

           write(77,'(3f8.3,2E15.3)')xx,Q2,Wsq,F2,F2_1

10      continue
  
        end

c----------------------------------------------------------------------------------------------

        real*8 function NMCF2d(x,Q2)
c       F2d used in NMC radiative correction 
c       NMC, P. Amaudruz et al., Nucl. Phys. B 371 (1992) 3.
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

        W=sqrt(Mp*Mp+Q2*(1./x-1.))
        GQ2=1.0/(1.0+Q2/0.71)**2
        S=1.0-exp(-a*(W-Wthr))
        xw=(Q2+ma2)/(Q2/x+mb2)
        F2DIS=(5.0/18.0*3.0/beta(eta(1),eta(2)+1.0)*xw**eta(1)*(1.-xw)**eta(2)
     >       +1.0/3.0*eta(3)*(1.-xw)**eta(4))*S

        F2RES=aa(5)**2*sqrt(GQ2**3)*exp(-(W-Mdelta)**2/gammac**2)

        eps=sqrt(((W+c)**2+Mp**2-Mpi**2)**2/(4.*(W+c)**2)-Mp**2)

        F2BG=aa(6)**2*sqrt(GQ2)*eps*exp(-b*(W-Wthr)**2)

        NMCF2d=(1.0-GQ2**2)*(F2DIS+F2RES+F2BG)
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
        real*8 xx,x,tmp,ser
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
c------------------------------------------------------------------------------------------------

