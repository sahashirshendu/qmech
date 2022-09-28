set term post color
set out 'hoa.ps'
set yrange [0:5]
plot 'hoa.txt' u 1:2 w l title 'Ground State', 'hoa.txt' u 1:3 w l title '1st Excited State', 'hoa.txt' u 1:4 w l title 'Potential'
