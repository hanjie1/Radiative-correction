###     GNUPLOT script
###     "Super ratio" R=R32/R31 where R32=F2(3He)/(2F2p+F2n) and R31=F2(3H)/(F2p+2F2n)
reset
set macros
set term pdfcairo enhanced

#data="F2dis_os1tm1ht1_Salme"
#data0="F2dis_os1tm1ht0_Salme"
#data00="F2dis_os1tm0ht0_Salme"
#data000="F2dis_os0tm0ht0_Salme"

data="F2dis_os1tm1ht1mec1_Dav18_He3Salme"
data0="F2dis_os1tm1ht0mec1_Dav18_He3Salme"
data00="F2dis_os1tm0ht0mec1_Dav18_He3Salme"
data000="F2dis_os0tm0ht0mec1_Dav18_He3Salme"
data0000="F2dis_os0tm0ht0mec0_Dav18_He3Salme"

#set xlabel "{/Times-Italic*1.2 x}"
set xlabel "x_{Bj}" offset 0,1
set ylabel "Super-ratio {/Times-Italic*1.15 R}" offset 1,0
#set key left Left reverse spacing 1.5 
set key left Left spacing 1.2 
set grid
set xrange[0.15:0.95]
set yrange[0.975:1.025]
#set title "The ratio R_{32}/R_{31} computed with Q^2=14*x_{Bj} GeV^2 and different assumptions on F_2^{p,n}"
set title "Ratio R=R_{32}({}^3He)/R_{31}({}^3H) computed with Q^2=14*x_{Bj} and different assumptions on F@_2^{p,n}"
set output 'super_r.pdf'
plot data using 1:( ($7/(2.0*$3+$4))/($6/($3+2.0*$4)) ) title 'Full model' with lines lt 1 lc rgb "red" lw 1.5,\
     data0 u 1:( ($7/(2.0*$3+$4))/($6/($3+2.0*$4)) ) title 'No HT correction' with lines lt 1 lc rgb "black" lw 1.5 dt "--",\
    data00 u 1:( ($7/(2.0*$3+$4))/($6/($3+2.0*$4)) ) title 'No TMC and HT' with lines lt 1 lc rgb "black" lw 1.5 dt "-.-",\
   data000 u 1:( ($7/(2.0*$3+$4))/($6/($3+2.0*$4)) ) title 'No OS, TMC and HT' with lines lt 1 lc rgb "black" lw 1.8 dt "...",\
  data0000 u 1:( ($7/(2.0*$3+$4))/($6/($3+2.0*$4)) ) title 'No OS, TMC, HT and MEC' with lines lt 1 lc rgb "violet" lw 1.8 dt "..."

### Plot #2
#set title "The ratio R=R_{32}({}^3He)/R_{31}({}^3H) computed with different models of 3He/3H spect. func."
#dat1="F2dis_os1tm1ht1_SS"
#dat2="F2dis_os1tm1ht1_KPSV"
#dat3="F2dis_os1tm1ht1_Salme"
#set output "super_r_spfn.pdf"
#plot dat1 u 1:( ($7/(2.0*$3+$4))/($6/($3+2.0*$4)) ) title 'SS spec. func.' with lines lt 1 lc rgb "black" lw 1.5 dt "--",\
#     dat3 u 1:( ($7/(2.0*$3+$4))/($6/($3+2.0*$4)) ) title 'Salme spec. func.' with lines lt 1 lc rgb "red" lw 1.5

