        program spectrum
        implicit none

        real*8 e0,ep,theta,Mp,Z,A
        real*8 F1,F2
        real*8 Qsq,Wsq,W1P,W2P,W1D,W2D,W1n,W2n,w1,w2
        real*8 xbj
        real*8 PI
        integer flag

        flag=1

        Z=1.0
        A=1.0
        e0=10.59
        ep=3.1
        Mp=0.93827231
        w1p=0.0
        w2p=0.0
        w1d=0.0
        w2d=0.0
        PI=3.1415927

c        open(unit=11,file='Bodek_F2d.dat')
        open(unit=11,file='Bodek_test.dat')
c        write(11,*) 'Ep    Theta    x    F2    F1'
        write(11,*) 'Qsq    Wsq    x    F2    F1'

c        do 100 theta=15, 36, 1
        do 100 xbj=0.18, 0.24, 0.01
c           Qsq=4.0*e0*ep*sin(theta*PI/180./2.)*sin(theta*PI/180./2.)
           Qsq=3.06
           Wsq=mp*mp+Qsq*(1./xbj-1)
           ep=e0-Qsq/(2.0*Mp*xbj)
c           xbj=Qsq/(2.*Mp*(e0-ep))
c           write(6,*) theta,Qsq,W2,xbj
           if(flag .eq. 1) then
            call ineft(Qsq,sqrt(Wsq),W1p,W2p,dble(1.0))
            call ineft(Qsq,sqrt(Wsq),W1D,W2D,dble(2.0))
            W1n=2.0*W1D-W1p
            W2n=2.0*W2D-W2p
            W1=Z*W1p+(A-Z)*W1n
            W2=Z*W2p+(A-Z)*W2n
            F1=W1*Mp
            F2=W2*(e0-ep)
           endif
           if(flag .eq. 2) then
            call gsmearing(Z,A,WSQ,QSQ,F1,F2)
            W1=F1/Mp
            W2=F2/(e0-ep)
           endif
c           write(11,*) ep,theta,xbj,F2,F1 
           write(11,*) Qsq,Wsq,xbj,F2,F1 
100     continue

        end
