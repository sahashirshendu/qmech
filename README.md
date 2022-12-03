# QMech
Quantum Mechanics Lab Code

- [Stern-Gerlach Experiment](./deflstern.f)
- [Sequential Stern-Gerlach Experiment](./seqstern.f)
- [1D Schrodinger Equation for Sine Potential](./schrost.f)
  - Plot
  ```gnuplot
  set term post
  set out 'schrost.ps'
  plot 'dschrost' w l, 'dschrost0' w l
  ```
- [Harmonic Oscillator](./ho.f)
  - Plot
  ```gnuplot
  set term post color
  set out 'ho.ps'
  set yrange [0:2]
  plot 'ho.txt' u 1:2 w l title 'Ground State', 'ho.txt' u 1:3 w l title 'Potential'
  ```
- [Anharmonic Oscillator](./hoa.f)
  - Plot
  ```gnuplot
  set term post color
  set out 'hoa.ps'
  plot 'hoa.txt' w l title 'Ground State'
  ```
- [ε(r) vs r](./hoe.f)
  - Plot
  ```gnuplot
  set term post color
  set out 'hoe.ps'
  plot 'hoe.txt' w l
  ```
- [Harmonic Oscillator (Numerov Method)](./hon.f)
  - Plot
  ```gnuplot
  set term post color
  set out 'hon.ps'
  plot 'hon.txt' u 1:2 w l
  ```
