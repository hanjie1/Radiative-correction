        program CheckA3
        implicit none

c       This is to check whether we should include resonance in F2D when
c       construct F2(He3) by comparing the EMC from model with Hall C
c       EMC data

        real*8 A,Z,x,Q2,Wsq
        real*8 W1D,W2D,F2D,F2HE3
        character*80 infile,COMMENT,outfile
        integer ii
        real*8 xi,EMC_HALLC,Estat,Esys,Fis,EMC_model
        real*8 EMC_KP
        EXTERNAL EMC_KP

        infile='22deg_he3_newiso_emc_rebin_x.dat'
        write(6,*) infile

        OPEN(UNIT=7,FILE=infile)
        READ(7,'(A80)') COMMENT

        outfile='OUT/CheckA3_He3_22deg.out'
        open(unit=66,file=outfile)
        write(66,*) 'x    Q2    WSQ   EMC_hallC   EMC_model'

        Z=2.0
        A=3.0

        do 99, ii=1,28
           x=0.0
           Q2=0.0
           READ(7,'(F7.4,3F9.4,4E12.4)',END=100) x,xi,Q2,WSQ,EMC_HALLC,Estat,Esys,Fis

           call ineft(Q2,sqrt(wsq),W1D,W2D,dble(2.0)) 
           F2D=W2D
           call ineft(Q2,sqrt(wsq),W1D,W2D,dble(3.0)) 
           F2HE3=W2D*EMC_KP(x,Z,A)*2.0/A

           EMC_model=F2HE3/F2D

           write(66,'(F7.4,2F9.4,5E12.4)') x,Q2,WSQ,EMC_HALLC,Estat,Esys,Fis,EMC_model         

99      continue
100     continue



        end

c---------------------------------------------------------------------------------------------

        real*8 function emc_KP(x,Z,A)
        real*8 x,Z,A

        emc_KP = 1.0
        if((Z .EQ. 1.0) .and. (A .eq. 3.0)) then
          emc_KP = 1.07251-1.61648*x+12.3626*x**2-65.5932*x**3+213.311*x**4
     >           -423.943*x**5+499.994*x**6-321.304*x**7+87.0596*x**8
          emc_KP = emc_KP*A/2.0
        endif

        if((Z .EQ. 2.0) .and. (A .eq. 3.0)) then
          emc_KP = 1.02967-0.135929*x+3.92009*x**2-21.2861*x**3
     >           +64.7762*x**4-129.928*x**5+169.609*x**6
     >           -127.386*x**7+41.0723*x**8
          emc_KP = emc_KP*A/2.0
        endif

        return
        end

