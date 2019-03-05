        program EMC_SLAC
        implicit none

        real*8 m_p /0.93827231/
        real*8 amuM,x,Q2,eps,nu,wsq,A,Z
        real*8 emciso,emccor,F_IS,F2_NP,CJ_f2n,CJ_f2p
        character*80  outfile
        integer ii
        real*8 emc_func_slac
        external emc_func_slac
        real*8 CJsfn
        external CJsfn

        outfile='OUT/SLAC_EMC_He3_CJ15.out'
        open(unit=66,file=outfile)
        write(66,*) 'x   F2A/F2D '

        A=3.0
        Z=2.0
        call setCJ(600)

        do 99 ii=1,80
           x=0.1+0.01*ii
           Q2=14.0*x

           emciso=emc_func_slac(x,A)

c           F2_NP=1-0.8*x
           CJ_f2p=CJsfn(1,x,sqrt(Q2))
           CJ_f2n=CJsfn(2,x,sqrt(Q2))
           F2_NP=CJ_f2n/CJ_f2p


           F_IS=(1+F2_NP)/(Z+(A-Z)*F2_NP)
           emccor=emciso/F_IS

           write(66,'(3F13.5)') x,emciso,EMCcor

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

