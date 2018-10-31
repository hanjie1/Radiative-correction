      PROGRAM CSCOMPQ2EPS
      IMPLICIT none
 
      real*8 ein(20000),e(20000),ep(20000),th(20000),nu(20000)
      real*8 q2(20000),w2(20000),cs(20000),css(20000),cserr(20000)
      real*8 flux,x(20000),eps,rat,erat,sys,kappa,sin2,tan2
      real*8 csmod(20000),A(20000),Z(20000),foc(20000)
      real*8 f1,f2,fl,f1qe,flqe,r,rqe,sigt,sigl,sigm,f2qe,f2mec  
      real*8 alpha,pi,pi2,mp,mp2,res,veff,chi2,q2min,q2max,w2max
      real*8 epsmin,epsmax

      integer i,j,k,ntot,opt,set(20000) 
      real*4 x4,a4,emcfac,fitemc
      logical thend/.false./
      logical coulomb/.true./
      logical new,wr
      LOGICAL GOODFIT/.true./  
c      real*8 norm(15)/
c     & 0.97866E+00,0.96776E+00,0.10363E+01,0.99254E+00,0.98301E+00,
c     & 0.10000E+01,0.97634E+00,0.10055E+01,0.98454E+00,0.99203E+00,
c     & 0.99209E+00,0.99553E+00,0.10000E+01,0.10125E+01,0.10060E+01 /
      real*8 norm(15)/
     &     1.0,1.0,1.0,1.0,1.0,
     &     1.0,1.0,1.0,1.0,1.0,
     &     1.0,1.0,1.0,1.0,1.0 /
      
      character*40 filename
      character*60 myfile
      integer nn

      opt = 0  !!! Include all components of cross section model  !!!
c      filename = '12C-it2n.dat'
c      filename = '27Al-it2n.dat'
c      filename = '56Fe-it2n.dat'
c      filename = '64Cu-it2n.dat'
       
c      filename = '40Ca.dat'      
c      filename = '48Ca.dat'

      filename = '4He.dat'
      
      read(5,*) q2min,q2max,epsmin,epsmax
      w2max = 100.0

      mp = .938272
      mp2 = mp*mp

      alpha = 1./137.036
      pi = 3.141593
      pi2 = pi*pi

      ntot = 0
      i=0



*****    Read in data   *****

      open(unit=15,file=filename,status='old')
      do i=1,16
        read(15,*)    !!!  Read header  !!!
      enddo

      i = 1
      dowhile(.not.thend)

       read(15,*,end=2000) Z(i),A(i),e(i),th(i),nu(i),cs(i),cserr(i),
     &   set(i)      

       ep(i) = e(i)-nu(i)

       sys = 0.03
       if(set(i).EQ.8.OR.set(i).EQ.11) then
         sys = 0.055
       elseif(set(i).EQ.12) then
         sys = 0.0225
       elseif(set(i).EQ.14) then
         sys = 0.001
       endif
       cserr(i) = sqrt(cserr(i)*cserr(i)+sys*sys*cs(i)*cs(i))

       if(coulomb) then
        call vcoul(A(i),Z(i),veff)
        foc(i) = 1.0D0 + veff/e(i)
        e(i) = e(i) + veff     
        ep(i) = ep(i) + veff
        cs(i) = cs(i)/foc(i)/foc(i)
        cserr(i) = cserr(i)/foc(i)/foc(i)
       endif
      
       q2(i) = 4.*e(i)*ep(i)*sin(th(i)*pi/180./2.)*sin(th(i)*pi/180./2.)
       w2(i) = mp2+2.*mp*nu(i)-q2(i)
       x(i) = q2(i)/2./mp/nu(i)

       ntot = i
       i=i+1 

      enddo

 2000 thend = .true.


CCC///    Now do sorting in whatever variable    ///CCC
       k = 0
       chi2 = 0.0
       do j=1,ntot
         wr = .false.            !!! if true then include data  !!!
         new = .true.
         wr = .true.
         k = k+1

         sin2 = sin(th(j)*pi/180./2.)*sin(th(j)*pi/180./2.)
         tan2 = sin2/(1.-sin2)
         eps = 1./(1. + 2.*(nu(j)*nu(j)+q2(j))/q2(j)*tan2)
         kappa = abs(w2(j)-mp2)/2./mp

         if(w2(j).LT.0.15.or.nu(j).LE.0.0.OR.q2(j).LT.q2min.OR.
     &         q2(j).GT.q2max.or.w2(j).GT.w2max.or.eps.LT.epsmin.
     &         or.eps.GT.epsmax.or.set(j).EQ.2) 
     &        wr = .false.

         if(set(j).EQ.13.AND.th(j).GE.74.0.AND.w2(j).GT.3.25) then
           wr = .false.
         endif

         flux = alpha*kappa/(2.*pi*pi*q2(j))*ep(j)/e(j)/(1.-eps) 

         cs(j) = norm(set(j))*cs(j)
         cserr(j) = norm(set(j))*cserr(j)

         cs(j) = cs(j)/flux/1000.0
         cserr(j) = cserr(j)/flux/1000.0

         if(wr) then

           call sfcross(w2(j),q2(j),A(j),Z(j),opt,sigt,sigl,f1,f2,fL)
         
           sigm = sigt+eps*sigl

           res = cs(j)-sigm
           chi2 = res*res/cserr(j)/cserr(j)

           rat = cs(j)/sigm
           erat = cserr(j)/cs(j)*rat

 
           write(6,3000) e(j),ep(j),th(j),w2(j),q2(j),
     &            eps,cs(j),cserr(j),sigm,rat,erat,set(j)
         endif  
       enddo

       nn=index(filename,'.')
       myfile='chi2_'//filename(1:nn-1)//'.out'
       open(unit=18,file=myfile,status="old", position="append")
       write(18,*) chi2


 3000  format(6f9.4,3f15.4,2f8.4,1i4)
 3200  format(6f9.4,2f13.4)
       
      end








