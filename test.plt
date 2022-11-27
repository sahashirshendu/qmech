set terminal cairolatex standalone pdf size 15cm,10cm
set out 'test.tex'
set style data lines
set xlabel '$x$'
set ylabel '$\psi(x)$'
set xrange [-4:4]
set yrange [0:]
plot 'hon.txt' t ''
