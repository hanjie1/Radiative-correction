        program EMCHE3
        implicit none

        real*8 m_p /0.93827231/
        real*8 amuM,x,Q2,eps,nu,wsq
        real*8 lamdaD /1.015/,lamdaHE /1.04/,nuHE2 /0.6156/
        real*8 A1,newQ2,newW2,newNU,EMCcor
        real*8 w1d,w2d,F2D,F2A
        character*80  outfile
        integer ii

        outfile='OUT/EMC_He3.out'
        open(unit=66,file=outfile)
        write(66,*) 'x   F2A/F2D '

        do 99 ii=1,99
           x=0.+0.01*ii
           Q2=14.0*x
           wsq=m_p**2+Q2*(1.0/x-1.0)

           nu=Q2/(2.0*m_p*x)
           call ineft(Q2,sqrt(wsq),W1D,W2D,dble(2.0))
           F2D=nu*W2D

           A1=2*log(Q2/0.25**2)/log(nuHE2/0.25**2)
           eps=(lamdaHE/lamdaD)**A1
           newQ2=eps*Q2

           newNU=newQ2/(2.0*m_p*x)
           newW2=m_p**2+newQ2*(1.0/x-1.0)
           call ineft(newQ2,sqrt(newW2),W1D,W2D,dble(2.0))
           F2A=W2D*newNU

           EMCcor=F2A/F2D
           write(66,'(2F13.5)') x,EMCcor

99      continue


        return
        end


