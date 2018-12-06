      PROGRAM PARMCONVERT
      
      real*8 temp(10),xval(100)
      real*8 parms(5)
      integer i,j,k


      open(unit=15,file='init_parms.dat',status='old')
      do i=1,100
        read(15,*) temp(1)  ! par #
        read(15,*) xval(i) ! starting value
        read(15,*) temp(2)   ! initial step (0 means fixed parm)
        read(15,*) temp(3) ! low limit
        read(15,*) temp(4) ! high limit
      enddo
      close(15)

      j=0 
      k=0
      do i=1,100
        if(mod(i,5).EQ.0.) then
           k = j*5+1
          write(6,2001) "     & ",xval(k),",",xval(k+1),",",xval(k+2)
     &       ,",",xval(k+3),",",xval(k+4),","

          j = j+1
        endif
      enddo 

      write(6,*)
      write(6,*)


 2001  format(a7,1e11.5,a1,1e11.5,a1,1e11.5,a1,1e11.5,a1,1e11.5,a1)

      end
