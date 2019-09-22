        program F2dp_Bodek
        implicit none

        real*8 xx(8),Q2(8)
        real*8 tmpX,tmpQ2,Wsq,Mp
        real*8 F2d,F2p,F1D,F1P
        real*8 F2dp,F2dp_err
        integer ii
        character*20 outfile
        
        DATA XX/0.17,0.19,0.22,0.25,0.29,0.32,0.34,0.38/
        DATA Q2/2.42,2.72,2.72,3.08,3.09,3.15,3.35,3.46/

        outfile='OUT/F2dp_Bodek.out'
        open(unit=66,file=outfile)
        write(66,*) 'x    Q2    F2d/F2p    error' 
            
        Mp=0.93827231
        do 99 ii=1,8
           tmpX=xx(ii)
           tmpQ2=Q2(ii)
           Wsq=Mp**2+tmpQ2*(1./tmpX-1.)

           call ineft(tmpQ2,sqrt(WSQ),F1p,F2p,dble(1.0))
           call ineft(tmpQ2,sqrt(WSQ),F1d,F2d,dble(2.0))

           F2dp=2.0*F2d/F2p
           F2dp_err=0.0

           write(66,'(2F7.2,2F10.5)') tmpX,tmpQ2,F2dp,F2dp_err
99      continue

        end
