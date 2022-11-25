      program numerov
      integer, parameter :: n = 100
      real, dimension(n) :: x, wfn
      real en, emin, emax, energy, enval, k2, dx, xmin, L
      common x, wfn, en, dx
      external enval

      open(1,file='hon.txt')
      L = 8
      xmin = - 0.5 * L
      emin = 0.0
      emax = 1.0

      dx = L / (n - 1)
      do i = 1, n
         x(i) = xmin + (i - 1) * dx
      end do

      call bisection(enval, emin, emax, energy)
      write(*,*) 'Energy =', energy

      do i = 1, n
      write (1, *) x(i), wfn(i)
      end do
      stop
      end

      function enval(energy)
      integer, parameter :: n = 100
      real, dimension(n) :: x, wfn, k2
      common x, wfn, en, dx

      en = energy
      do i = 1, n
        k2(i) = 2.0 * (en - 0.5 * x(i) ** 2)
      end do
      wfn(1) = 0.0
      wfn(2) = 1.0e-5
      do i = 2, n - 1, 1
        cm = 1.0 + dx ** 2 / 12.0 * k2(i - 1)
        cz = 2.0 * (1.0 - 5.0 / 12.0 * dx ** 2 * k2(i))
        cp = 1.0 + dx ** 2 / 12.0 * k2(i + 1)
        wfn(i + 1) = (cz * wfn(i) - cm * wfn(i - 1)) / cp
      end do
      enval = wfn(n)
      return
      end

      subroutine bisection(f, a, b, x)
      do 11 i = 1, 100
      x = (a + b) / 2
      if ((f(x) * f(b)) .le. 0.0) then
      a = x
      else
      b = x
      end if
      err = abs((b - a) / b)
      if (err .lt. 0.000001) exit
11    continue
      return
      end
