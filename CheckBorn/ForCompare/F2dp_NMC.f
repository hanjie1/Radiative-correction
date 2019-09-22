        program F2dp_NMC
        implicit none

        real xx(8),Q2(8)
        real tmpX,tmpQ2
        real F2d,F2p,F2d_loerr,F2d_hierr,F2p_loerr,F2p_hierr
        real F2dp,F2dp_loerr,F2dp_hierr
        integer ii
        character*20 outfile
        
        DATA XX/0.17,0.19,0.22,0.25,0.29,0.32,0.34,0.38/
        DATA Q2/2.42,2.72,2.72,3.08,3.09,3.15,3.35,3.46/

        outfile='OUT/F2dp_NMC.out'
        open(unit=66,file=outfile)
        write(66,*) 'x    Q2    F2d/F2p    low_error    high_error' 
            
        do 99 ii=1,8
           tmpX=xx(ii)
           tmpQ2=Q2(ii)

           call F2NMC_new(1,tmpX,tmpQ2,F2p,F2p_loerr,F2p_hierr)
           call F2NMC_new(2,tmpX,tmpQ2,F2d,F2d_loerr,F2d_hierr)
           
           F2dp=2.0*F2d/F2p
           F2dp_loerr=F2dp*sqrt((F2p_loerr/F2p)**2+(F2d_loerr/F2d)**2)
           F2dp_hierr=F2dp*sqrt((F2p_hierr/F2p)**2+(F2d_hierr/F2d)**2)

           write(66,'(2F7.2,3F10.5)') tmpX,tmpQ2,F2dp,F2dp_loerr,F2dp_hierr
99      continue

        end
