      SUBROUTINE RESCSP(w2,q2,sigtp,siglp)
CCCCC  test version from 12/08/18 to address Q^2 > 9 behavior at W^2 < 9 CCCCC
      
      IMPLICIT none

      real*8 w2,q2,sigtp,siglp
      real*8 xvalp(100),xval1(50),xvalL(50)
      Integer i
      real*8 sigtdis,sigLdis,w2max,w2min

      data xvalp / 
     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     & 0.10936E+00,0.37944E+00,0.77805E+01,0.42291E+01,0.12598E+01,
     & 0.21242E+01,0.63351E+01,0.68232E+04,0.33521E+05,0.25686E+01,
     & 0.60347E+00,0.21240E+02,0.55746E-01,0.24886E+01,0.23305E+01,
     & -.28789E+00,0.18607E+00,0.63534E-01,0.19790E+01,-.56175E+00,
     & 0.38964E+00,0.54883E+00,0.22506E-01,0.46213E+03,0.19221E+00,
     & 0.19141E+01,0.24606E+03,0.67469E-01,0.13501E+01,0.12054E+00,
     & -.89360E+02,0.20977E+00,0.15715E+01,0.90736E-01,-.38495E-02,
     & 0.10362E-01,0.19341E+01,0.38000E+00,0.34187E+01,0.14462E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.29414E+02,0.19910E+02,0.22587E+00,
     & 0.00000E+00,0.00000E+00,0.38565E+04,0.65717E+00,0.00000E+00,
     & 0.15792E+03,0.97046E+02,0.31042E+00,0.00000E+00,0.42160E+01,
     & 0.38200E-01,0.12182E+01,0.00000E+00,0.13764E+02,0.31393E+00,
     & 0.29997E+01,0.00000E+00,0.55124E+01,0.53743E-01,0.13091E+01,
     & 0.00000E+00,0.86746E+02,0.40864E-04,0.40294E+01,0.31285E+01,
     & 0.33403E+00,0.49623E+01,0.00000E+00,0.00000E+00,0.11000E+02,
     & 0.18951E+01,0.51376E+00,0.00000E+00,0.42802E+01,0.00000E+00 /

      do i=1,50
        xval1(i) = xvalp(i)
        xvalL(i) = xvalp(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)  !!! use same masses for L as T !!!
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)

      if(q2.LE.8) then
         w2max = 8.0
         w2min = 6.0
      else
         w2max = 8.0*(1.0-0.1*(q2-8.0))
         w2min = 6.0*(1.0-0.1*(q2-8.0))
      endif
      
      call resmodp(1,w2,q2,xval1,sigTp)
      call resmodp(2,w2,q2,xvalL,sigLp)
c      call disp(w2,q2,sigTdis,sigLdis)
      if(w2.GT.w2max) then 
c        write(6,*) "resmod: ",w2,q2,sigTp,sigLp
        call disp(w2,q2,sigTdis,sigLdis)
c        write(6,*) "dismod: ",w2,q2,sigtdis,sigLdis
        if(w2.LE.w2max) then
         sigTp = (w2max-w2)*sigTp+(w2-w2min)*sigTdis
         sigLp = (w2max-w2)*sigLp+(w2-w2min)*sigLdis
         sigTp = sigTp/2.0
         sigLp = sigLp/2.0
         else
          sigTp = sigTdis
          sigLp = sigLdis
        endif
      endif

      return
      end






