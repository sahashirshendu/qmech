      program deflection
c Deflection in Stern-Gerlach Experiment
      double precision k, bm, bg, l, r, T, g, d
      k = 1.381e-23 ! Boltzmann Constant
      bm = 9.274e-24 ! Bohr Magneton
      print *, "Magnetic field gradient (T/m) ="
      read *, bg
      print *, "Length of magnet (m) :"
      read *, r
      print *, "Distance (m) :"
      read *, l
      print *, "Oven temperature (K) :"
      read *, T
      print *, "g-factor :"
      read *, g
      d = 0.5 * bg * bm * g * r / (4 * k * T) * (r / 2 + l)
      print *, "Deflection (m) =", d
      end
