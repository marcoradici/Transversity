 set term postscript enhanced dashed color 'Helvetica,18' size 10,7
 set out "comparison.eps"

 set style fill transparent solid 0.5

 set xlabel "x"
 set ylabel "u@^v_{NLO}(x,Q) / u@^v_{LO}(x,Q)"

 f(x) = 1

 set logscale x

 set title "Transversity collinear evolution Q = 10 GeV"
 plot [0.01:0.9][0.7:1.1] \
 "comparison.dat" index 0 u 1:($3/$2) smooth bezier dt 1 lw 4 lc rgb "red" t "Time-like (FFs)", \
 "comparison.dat" index 1 u 1:($3/$2) smooth bezier dt 1 lw 4 lc rgb "blue" t "Space-like (PDFs)", f(x) lw 4 lc rgb "black" notitle
