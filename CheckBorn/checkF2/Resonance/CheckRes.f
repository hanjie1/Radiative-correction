        program CheckRes
        implicit none

        real*8 Q2,xx,WSQ
        real*8 W1,W2
        integer ii

        do 10 ii=0,40
           xx=0.1+ii*0.02
           Q2=3+ii*0.1
           wsq=0.938**2+Q2*(1.0/xx-1)

           call ineft(Q2,sqrt(wsq),w1,w2,dble(2.0))

10      continue
  
        end
