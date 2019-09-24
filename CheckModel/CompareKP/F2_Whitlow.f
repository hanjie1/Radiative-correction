        program F2dp_Whitlow
        implicit none

        real*8 xx,Q2
        real tmpX,tmpQ2
        real F2d,F2p,F2d_err,F2p_err,ST,SY
        real F2d_data,F2p_data,F2n_data
        REAL slope,dslope
        integer MODEL,ii
        character*100 outfile,infile
        logical GOODFIT
        
        MODEL=12
        infile='/home/hanjie/work/MARATHON/RadCor/KP_R/F2dis_os1tm1ht1mec1_Dav18_He3Salme'
        open(unit=77,file=infile)
        outfile='OUT/F2dp_Whitlow.out'
        open(unit=66,file=outfile)
        write(66,*) 'x    Q2    F2d/F2p    error' 
            
        do 99 ii=1,1200
           read(77,'(E11.4,4E12.4)',END=100)tmpX,tmpQ2,F2p_data,F2n_data,F2d_data
           if((tmpX .gt.0.8).or.(tmpX .lt. 0.15))goto 99  
           
           call F2GLOB(tmpX,tmpQ2,'P',MODEL,F2p,ST,SY,slope,dslope,GOODFIT) 
           F2p_err=sqrt(ST**2+SY**2)           
           call F2GLOB(tmpX,tmpQ2,'D',MODEL,F2d,ST,SY,slope,dslope,GOODFIT) 
           F2d=F2d
           F2d_err=sqrt(ST**2+SY**2)           


           write(66,'(2F9.4,6F10.5)')tmpX,tmpQ2,F2p,F2p_err,F2d,F2d_err,F2p_data,F2d_data
99      continue

100     end
