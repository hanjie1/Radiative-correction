###     GNUPLOT script
###     EMC ratios R2=F2D/F20 and R32=F2(3He)/(2F2p+F2n) and R31=F2(3H)/(F2p+2F2n)
reset
set macros
set title "Ratios R_2, R_{31} and R_{32} computed for Q^2=14*x"
set term pdfcairo enhanced

set output 'r2r3.pdf'
data="F2dis_os1tm1ht1mec1_Dav18_He3Salme"
data0="F2dis_os0tm0ht0mec0_Dav18_He3Salme"

set xlabel "x_{Bj}" offset 0,0.5
set ylabel "Ratio" offset 1,0
set key left Left reverse spacing 1.5
set grid
set xrange [0.15:0.9]
set yrange [0.9:1.3]
plot \
 data u 1:(2*$5/($3+$4))  title 'R_2 = F@_{/*0.8 2}^{2H}/(F@_{/*0.8 2}^p+F@_{/*0.8 2}^n)' w lines lt 1 lw 1 lc rgb "black",\
 data u 1:(3*$6/($3+2*$4)) title 'R_{31} = F_{/*0.8 2}^{3H}/(F@_{/*0.8 2}^p+2F@_{/*0.8 2}^n)' w lines lt 1 lw 1 lc rgb "blue",\
 data u 1:(3*$7/(2*$3+$4)) title 'R_{32} = F_{/*0.8 2}^{3He}/(2F@_{/*0.8 2}^p+F@_{/*0.8 2}^n)' w lines lt 1 lw 1 lc rgb "red"
#\
# data0 u 1:(2*$5/($3+$4))  title 'R_2 without TMC and HT' w lines lt 1 lw 1 lc rgb "black" dt '--',\
# data0 u 1:(3*$6/($3+2*$4)) title 'R_{31} without TMC and HT' w lines lt 1 lw 1 lc rgb "blue" dt '--',\
# data0 u 1:(3*$7/(2*$3+$4)) title 'R_{32} without TMC and HT' w lines lt 1 lw 1 lc rgb "red" dt '--'

### Plot #2
#set title "EMC effect in {}^3H and {}^3He with different spectral functions"
#dat1="F2dis_os1tm1ht1_SS"
#dat2="F2dis_os1tm1ht1_Salme"
#set output "r3_spfn.pdf"
#plot \
# dat1 u 1:(3*$6/($3+2*$4)) title 'R_{31} with Hannover spec. func.' w lines lt 1 lw 1 lc rgb "blue",\
# dat2 u 1:(3*$6/($3+2*$4)) title 'R_{31} with Salme spec. func.' w lines lt 1 lw 1 lc rgb "blue" dt '--',\
# dat1 u 1:(3*$7/(2*$3+$4)) title 'R_{32} - Hannover' w lines lt 1 lw 1 lc rgb "red",\
# dat2 u 1:(3*$7/(2*$3+$4)) title 'R_{32} - Salme' w lines lt 1 lw 1 lc rgb "red" dt '--'
