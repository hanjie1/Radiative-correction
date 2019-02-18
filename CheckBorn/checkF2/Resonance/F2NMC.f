!---------------------------------------------------------------------
! -*-Mode: Fortran; compile-command: "f77 -o f2nmc.o -c f2nmc.f"; -*-
       Subroutine F2NMC_new(T,x4,qsq4,F2,err_lo,err_hi)
! This subroutine returns a value for F2 from the NMC parametrisation   
! in CERN_PPE/95-138  Sept 4, 1995  Proton and Deuteron F2 Structure
! Functions 
! in deep inelastic muon scattering   P. Arneodo et. al.               
!  Published in Phys.Lett.B364:107-115,1995 
!   e-Print Archive: hep-ph/9509406 
!                   
      IMPLICIT NONE                                                          
      real x4,qsq4,f2,err_lo,err_hi ,FNP_NMC
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron
                                ! 3 = neutron

      IF(T.LE.2) THEN
       CALL F2NMC_NEW_DO(T,x4,qsq4,F2,err_lo,err_hi)
      ELSE  !neutron
       CALL F2NMC_NEW_DO(1,x4,qsq4,F2,err_lo,err_hi)  !proton
       F2 = F2 * FNP_NMC(X4,QSQ4)
       ERR_LO = ERR_LO * FNP_NMC(X4,QSQ4)
       ERR_HI = ERR_HI * FNP_NMC(X4,QSQ4)
      ENDIF
      RETURN
      END
      SUBROUTINE F2NMC_NEW_DO(T,x4,qsq4,F2,err_lo,err_hi)                                         
      IMPLICIT NONE
      real x4,qsq4,f2,f2l,f2u,err_lo,err_hi
      real*8 a(7,2)
     >  /-0.02778,  2.926,   1.0362, -1.840, 8.123, -13.074, 6.215,
     >   -0.04858,  2.863,    .8367, -2.532, 9.145, -12.504, 5.473/
      real*8 b(4,2)/ 0.285,   -2.694,   0.0188,  0.0274,
     >              -0.008,   -2.227,   0.0551,  0.0570/
      real*8 c(4,2)/-1.413,    9.366, -37.79,   47.10,
     >              -1.509,    8.553, -31.20,   39.98/

!lower limits
      real*8 al(7,2)
     >  /-0.01705,  2.851,   0.8213, -1.156, 6.836, -11.681, 5.645,
     >   -0.02732,  2.676,    .3966, -0.608, 4.946,  -7.994, 3.686/
      real*8 bl(4,2)/ 0.325,   -2.767,   0.0148,  0.0226,
     >                0.141,   -2.464,   0.0299,  0.0396/
      real*8 cl(4,2)/-1.542,   10.549, -40.81,   49.12,
     >               -2.128,   14.378, -47.76,   53.63/

!upper limits
      real*8 au(7,2)
     >  /-0.05711,  2.887,   0.9980, -1.758, 7.890, -12.696, 5.992,
     >    -0.04715,  2.814,    .7286, -2.151, 8.662, -12.258, 5.452/
      real*8 bu(4,2)/ 0.247,   -2.611,   0.0243,  0.0307,
     >               -0.048,  -2.114,   0.0672,  0.0677/
      real*8 cu(4,2)/-1.348,    8.548, -35.01,   44.43,
     >               -1.517,    9.515, -34.94,   44.42/

      real Lam/0.25/
      real*8 Qsqo/20.0/                                                   
      real*8 log_term                                                     
      real*8 A_x,B_x,C_x,dl,lam2
      real*8 x,qsq,z,z2,z3,z4,f28
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron                          
      logical first/.true./



      if(first) then
        dl =dlog(qsqo/lam**2)
        lam2 =lam**2
        first=.false.
      endif

      x=x4
      z=1.-x
      z2=z*z
      z3=z*Z2
      z4 =z3*z

      qsq=qsq4
      A_x=x**a(1,t)*z**a(2,t)*(a(3,t)+a(4,t)*z+a(5,t)*z2
     >    +a(6,t)*z3+a(7,t)*z4)
      B_x=b(1,t)+b(2,t)*x+b(3,t)/(x+b(4,t))
      C_x=c(1,t)*x+c(2,t)*x**2+c(3,t)*x**3+c(4,t)*x**4
      Log_term = dlog(qsq/lam2)/dl
      if(Log_term.lt.0) Log_term=0.     !new on 3/15/00  SER
      F28=A_x*(log_term)**B_x*(1+C_x/qsq)
      f2 = f28

      A_x=x**au(1,t)*z**au(2,t)*(au(3,t)+au(4,t)*z+au(5,t)*z2
     >    +au(6,t)*z3+au(7,t)*z4)
      B_x=bu(1,t)+bu(2,t)*x+bu(3,t)/(x+bu(4,t))                             
      C_x=cu(1,t)*x+cu(2,t)*x**2+cu(3,t)*x**3+cu(4,t)*x**4
      f2u   =A_x*(log_term)**B_x*(1+C_x/qsq)

      A_x=x**al(1,t)*z**al(2,t)*(al(3,t)+al(4,t)*z+al(5,t)*z2
     >    +al(6,t)*z3+al(7,t)*z4)
      B_x=bl(1,t)+bl(2,t)*x+bl(3,t)/(x+bl(4,t))                             
      C_x=cl(1,t)*x+cl(2,t)*x**2+cl(3,t)*x**3+cl(4,t)*x**4
      f2l   =A_x*(log_term)**B_x*(1+C_x/qsq)

      err_hi = f2u-f2
      err_lo = f2-f2l

      return
      end

!---------------------------------------------------------------------
      REAL FUNCTION FNP_NMC(X,QSQ)
!---------------------------------------------------------
! NMC FIT TO NMC,SLAC,NMC data in CERN-PPE/91-167
! No Fermi Motion Corrections
!  Steve Rock 1 Feb 1996
!-------------------------------------------------------------
      IMPLICIT NONE
      REAL X,QSQ,A,B,X2,X3

      X2 = X*X
      X3 = X2*X
      A = 0.979 -1.692*X +2.797*X2 -4.313*X3 +3.075*X3*X
C$$      B = -.171*X2 + .277*X3
      B = -.171*X2 + .244*X3  ! replaced 10/22/97 by correct value on x3
      FNP_NMC = A *(1+X2/QSQ)*(QSQ/20.)**B
      RETURN
      END


