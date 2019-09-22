        program F2dp_Whitlow
        implicit none

        real xx(8),Q2(8)
        real tmpX,tmpQ2
        real F2d,F2p,F2d_err,F2p_err,ST,SY
        real F2dp,F2dp_err
        REAL slope,dslope
        integer MODEL,ii
        character*20 outfile
        logical GOODFIT
        
        DATA XX/0.17,0.19,0.22,0.25,0.29,0.32,0.34,0.38/
        DATA Q2/2.42,2.72,2.72,3.08,3.09,3.15,3.35,3.46/

        MODEL=12
        outfile='OUT/F2dp_Whitlow.out'
        open(unit=66,file=outfile)
        write(66,*) 'x    Q2    F2d/F2p    error' 
            
        do 99 ii=1,8
           tmpX=xx(ii)
           tmpQ2=Q2(ii)

           call F2GLOB(tmpX,tmpQ2,'P',MODEL,F2p,ST,SY,slope,dslope,GOODFIT) 
           F2p_err=sqrt(ST**2+SY**2)           
           call F2GLOB(tmpX,tmpQ2,'D',MODEL,F2d,ST,SY,slope,dslope,GOODFIT) 
           F2d=F2d*2.0
           F2d_err=2.0*sqrt(ST**2+SY**2)           

           F2dp=F2d/F2p
           F2dp_err=F2dp*sqrt((F2p_err/F2p)**2+(F2d_err/F2d)**2)

           write(66,'(2F7.2,2F10.5)') tmpX,tmpQ2,F2dp,F2dp_err
99      continue

        end
