set term post color
set out 'ho.ps'
set yrange [0:2]
plot 'ho.txt' u 1:2 w l title 'Ground State', 'ho.txt' u 1:3 w l title 'Potential'
