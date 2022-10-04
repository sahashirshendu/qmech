set term post color
set out 'hon.ps'
set yrange [0:2]
plot 'hon.txt' u 1:2 w l title 'Ground State', 'hon.txt' u 1:3 w l title 'Potential'
