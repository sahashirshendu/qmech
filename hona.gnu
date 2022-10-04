set term post color
set out 'hona.ps'
set yrange [0:5]
plot 'hona.txt' u 1:2 w l title 'Ground State', 'hona.txt' u 1:3 w l title '1st Excited State', 'hona.txt' u 1:4 w l title 'Potential'
