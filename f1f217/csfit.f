      SUBROUTINE CSFIT(w2,q2,A,Z,XVALC,sigt,sigL)
      IMPLICIT none

      real*8 e,ep,th,q2,w2,x,cs,flux,kappa,sin2,tan2,csmod
      real*8 f1,fl,f1qe,flqe,r,rqe,sigt,sigl,f1mec,sigm,f2,f2qe,f2mec
      real*8 alpha,pi,pi2,mp,mp2,res,veff,foc,Z,A,xvalc(40)
      real*8 psip,fy1,fy2,int1,int2,rat,f1t,f2t,fLt
      integer i,j,k,ntot 
      LOGICAL GOODFIT/.true./  
      character*40 filename
      mp = .938272
      mp2 = mp*mp

      alpha = 1./137.036
      pi = 3.141593
      pi2 = pi*pi

      x = q2/(w2-mp2+q2)
 

      int1 = 0.0D0
      int2 = 0.0D0
      rat = 1.0D0

CCC  NEXT bit only needed if fitting scaling function CCC
      do i=1,120        
        psip = -2.0+0.06*i
       FY1 = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip))
       FY2 = xvalc(7)/ (1. + xvalc(8)**2 * (psip + xvalc(9))**2) / 
     >   (1. + exp(xvalc(10) * psip))
       int1 = int1+fy1  
       int2 = int2+fy2
      enddo
      rat= int1/int2
CCC

      call gsmearing(Z,A,w2,q2,xvalc,f1,f2,fL)

      r = fL/2.0D0/x/f1

      call f1f2qe17(Z,A,q2,w2,xvalc,f1qe,f2qe)
      call MEC2017(Z,A,w2,q2,xvalc,f1mec)

      f1qe = f1qe*rat
      f2qe = f2qe*rat

      f1 = f1 + f1qe + f1mec
      f2mec = 2.*x*f1mec/(1.+4.*x*x*mp2/q2) 
      f2 = f2 + f2qe + f2mec
      fL = (1.+4.*x*x*mp2/q2)*f2-2.0*x*f1  
      sigt = f1
      sigl = fL/2./x


 2000 format(6f10.4)
          
      return

      end
      




