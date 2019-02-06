###     GNUPLOT script
###     EMC ratio F2(3He)/F2(D)
###
reset
set macros
# Isoscalarity correction factor
# b=(N-Z)/A
# r=F2n/F2p
fis(b,r)=1.0/(1.0 - b*(1.0-r)/(1.0+r))
#
set term pdfcairo enhanced
set title "Ratios {}^3He/D and {}^3H/D computed with Q^2=14*x \
\nNote that F@_2^A_{is} = F@_2^A(A/2)(F@_2^p + F@_2^n)/(Z F@_2^p + N F@_2^n) with A=Z+N the nucleon number"

dat="F2dis_os1tm1ht1mec1_Dav18_He3Salme"
set output "H3_He3_D.pdf"

set xlabel "x_{Bj}" offset 0,1
set ylabel "F@_2^A / F@_2^D" offset 1,0
set key left Left reverse spacing 1.5
set grid
set xrange[0.15:0.95]
set yrange[0.8:1.3]
b32=(1.0 - 2.0)/3.0
b31=(2.0 - 1.0)/3.0
plot \
  dat u 1:($6/$5)  title 'F@_2^{3H} / F@_2^D' with lines lt 1 lw 1 lc rgb "blue" dt "--",\
  dat u 1:(fis(b31,$4/$3)*$6/$5) title 'F@_{2is}^{3H} / F@_2^D' w l lt 1 lw 1 lc rgb "blue",\
  dat u 1:($7/$5)  title 'F@_2^{3He} / F@_2^D' with lines lt 1 lw 1 lc rgb "red" dt "--",\
  dat u 1:(fis(b32,$4/$3)*$7/$5) title 'F@_{2is}^{3He} / F@_2^D' w l lt 1 lw 1 lc rgb "red"
