      program numerov
      integer, parameter :: n = 100
      real, dimension(n) :: x, wfn
      real en, e1, e2, energy, enval, k2, e, dx, L
      common x, wfn, en, dx
      external enval

      open(1,file='hon.txt')
      L = 8
      e1 = 0
      e2 = 1

      dx = L / (n - 1)
      do i = 1, n
         x(i) = - .5 * L + (i - 1) * dx
      end do

      do 11 i = 1, 100
      e = (e1 + e2) / 2
      if ((enval(e) * enval(e2)) .le. 0.0) then
      e1 = e
      else
      e2 = e
      end if
      if (abs((e2 - e1) / e2) .lt. 0.000001) exit
11    continue
      write(*,*) 'Energy =', e

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
