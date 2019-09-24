        program NP_CJ
        implicit none

        real*8 xx(8),Q2(8)
        real*8 tmpX,tmpQ2
        real*8 F2n,F2p,F2n_err,F2p_err
        real*8 F2np,F2np_err
        integer ii
        character*20 outfile
        real*8 CJsfn
        external CJsfn
        
        DATA XX/0.17,0.19,0.22,0.25,0.29,0.32,0.34,0.38/
        DATA Q2/2.42,2.72,2.72,3.08,3.09,3.15,3.35,3.46/

        outfile='OUT/NP_CJ_14x.out'
        open(unit=66,file=outfile)
        write(66,*) 'x    Q2    F2n/F2p    error' 
        
        call setCJ(600)    
        do 99 ii=1,8
           tmpX=xx(ii)
c          tmpQ2=Q2(ii)
           tmpQ2=xx(ii)*14.0

           f2p=CJsfn(1,tmpx,sqrt(tmpQ2))
           f2n=CJsfn(2,tmpx,sqrt(tmpQ2))


           F2np=F2n/F2p
           F2np_err=0.0

           write(66,'(2F7.2,2F10.5)') tmpX,tmpQ2,F2np,F2np_err
99      continue

        end
