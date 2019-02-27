        program CheckHE3
        implicit none

c       This is to check whether we should include resonance in F2D when
c       construct F2(He3) by comparing the He3 cross section from Hall C 
c       data

        real*8 A,Z,x,Q2,Wsq,Ep,R,sig,esig,F2A,eF2A
        real*8 W1D,W2D,F2D,F2HE3
        character*80 infile1,COMMENT,outfile
        integer ii
        real*8 xi,F2_model,F2_KP,nu
        real*8 m_p /0.93827231/
        real*8 EMC_KP
        EXTERNAL EMC_KP

        infile1='18deg_all_he3.xbin.xsec'
        write(6,*) infile1

        OPEN(UNIT=7,FILE=infile1)
        READ(7,'(A80)') COMMENT

        outfile='OUT/XS_18deg.out'
        open(unit=66,file=outfile)
        write(66,*) 'x    Q2    WSQ   F2(He3)   err   F2A_model    F2A_KP'

        Z=2.0
        A=3.0

        do 99, ii=1,28
           x=0.0
           Q2=0.0
           READ(7,'(F7.4,4F9.4,5E12.4)',END=100) Ep,x,xi,WSQ,Q2,sig,esig,R,F2A,eF2A

           call ineft(Q2,sqrt(wsq),W1D,W2D,dble(2.0)) 
           nu=Q2/(2.0*m_p*x)
           F2_KP=W2D*EMC_KP(x,Z,A)*2.0*nu
c           call ineft(Q2,sqrt(wsq),W1D,W2D,dble(3.0)) 
           call F2GLOB(X,Q2,'D',12,F2D)
c           F2_model=W2D*EMC_KP(x,Z,A)*2.0*nu
           F2_model=F2D*EMC_KP(x,Z,A)*2.0
          


           write(66,'(F7.4,2F9.4,5E12.4)') x,Q2,WSQ,F2A,eF2A,F2_model,F2_KP      

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
          emc_KP = EMC_KP*3.0/2.0
        endif

        if((Z .EQ. 2.0) .and. (A .eq. 3.0)) then
          emc_KP = 1.02967-0.135929*x+3.92009*x**2-21.2861*x**3
     >           +64.7762*x**4-129.928*x**5+169.609*x**6
     >           -127.386*x**7+41.0723*x**8
          emc_KP = EMC_KP*3.0/2.0
        endif

        return
        end

