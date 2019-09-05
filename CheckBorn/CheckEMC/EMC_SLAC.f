        program EMC_SLAC
        implicit none

        real*8 m_p /0.93827231/
        real*8 amuM,x,Q2,eps,nu,wsq,A,Z
        REAL*8 F2_NP_NMC,F2_NP_SLAC,FIS_NMC,FIS_SLAC
        real*8 emciso,emccor,F_IS,F2_NP,CJ_f2n,CJ_f2p
        REAL*8 EMC_NMC,EMCSLAC
        character*80  outfile
        integer ii
        real*8 emc_func_slac
        external emc_func_slac
        real*8 CJsfn
        external CJsfn
        real*8 F2NP_NMC
        external F2NP_NMC

        outfile='OUT/Compare_EMC_H3.out'
        open(unit=66,file=outfile)
        write(66,*) 'x   F2A/F2D_CJ   NMC     SLAC'

        A=3.0
        Z=1.0
        call setCJ(600)

        do 99 ii=0,70
           x=0.15+0.01*ii
           Q2=14.0*x

           emciso=emc_func_slac(x,A)

c           F2_NP=1-0.8*x
           CJ_f2p=CJsfn(1,x,sqrt(Q2))
           CJ_f2n=CJsfn(2,x,sqrt(Q2))
           F2_NP=CJ_f2n/CJ_f2p

           F2_NP_SLAC=1-0.8*x
           F2_NP_NMC=F2NP_NMC(x,Q2)

           F_IS=(1+F2_NP)/(Z+(A-Z)*F2_NP)
           emccor=emciso/F_IS

           FIS_NMC=(1+F2_NP_NMC)/(Z+(A-Z)*F2_NP_NMC)
           emc_NMC=emciso/FIS_NMC

           FIS_SLAC=(1+F2_NP_SLAC)/(Z+(A-Z)*F2_NP_SLAC)
           emcSLAC=emciso/FIS_SLAC
           

           write(66,'(5F13.5)') x,emciso,EMCcor,emc_nmc,emcslac

99      continue


        return
        end

c-------------------------------------------------------------------------------------------
        real*8 function emc_func_slac(x,A)
        real*8 x,A,atemp
        real*8 alpha,C
!       Javier EMC fit for isoscalar nuclei, Phys. Rev. D 49 (4348) 1994

        atemp = A
!       if(A.eq.4.or.A.eq.3) then  ! emc effect is more like C for these
!       2...
!          atemp = 12
!       endif

        alpha = -0.070+2.189*x - 24.667*x**2 + 145.291*x**3
     >         -497.237*x**4 + 1013.129*x**5 - 1208.393*x**6
     >         +775.767*x**7 - 205.872*x**8

        C = exp( 0.017 + 0.018*log(x) + 0.005*log(x)**2)
        
        emc_func_slac = C*atemp**alpha
        return
        end


!---------------------------------------------------------------------
!---------------------------------------------------------------------
        real*8 function F2NP_NMC(x,Q2)
        real*8 x,Q2,AX,BX

        AX=0.979-1.692*x+2.797*x**2-4.313*x**3+3.075*x**4
        BX=-0.171*x+0.244*x**2
        F2NP_NMC=AX*((Q2/20.0)**BX)*(1+x**2/Q2)

        return
        end


