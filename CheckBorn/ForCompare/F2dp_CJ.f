        program F2dp_CJ
        implicit none

        real*8 xx(8),Q2(8)
        real*8 tmpX,tmpQ2
        real*8 F2d,F2p,F2d_err,F2p_err
        real*8 F2dp,F2dp_err
        integer ii
        character*20 outfile
        real*8 CJsfn
        external CJsfn
        
        DATA XX/0.17,0.19,0.22,0.25,0.29,0.32,0.34,0.38/
        DATA Q2/2.42,2.72,3.10,3.48,4.02,4.50,4.91,5.32/


        outfile='OUT/F2dp_CJ.out'
        open(unit=66,file=outfile)
        write(66,*) 'x    Q2    F2d/F2p    error' 
        
        call setCJ(600)    
        do 99 ii=1,8
           tmpX=xx(ii)
           tmpQ2=Q2(ii)

           f2p=CJsfn(1,tmpx,sqrt(tmpQ2))
           f2d=CJsfn(3,tmpx,sqrt(tmpQ2))


           F2dp=2.0*F2d/F2p
           F2dp_err=0.0

           write(66,'(2F7.2,2F10.5)') tmpX,tmpQ2,F2dp,F2dp_err
99      continue

        end
