set term post color
set out 'hoa.ps'
plot 'hoa.txt' u 1:2 w l title 'Ground State', 'hoa.txt' u 1:3 w l title 'Potential'
