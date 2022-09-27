set term post
set out 'ho.ps'
set yrange [0:2]
plot 'ho.txt' u 1:2 w l, 'ho.txt' u 1:3 w l
