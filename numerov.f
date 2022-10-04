      program lho
        implicit none
        integer np, ns, ierr, i, j, r
        real L, xm, dx, p, pot, de, hbar, m, k
        parameter(np = 100, ns = 10, L = 10.)
        parameter(hbar = 1d0, m = 1d0, k = 1d0)
        real a(np, np), H(np, np), sum, csint, normc
        real idz(np, np), idm(np, np), idp(np, np)
        real am(np, np), bm(np, np), bi(np, np), vm(np, np)
        real d(np), e(np), v(np), xr(np)
        external pot

        xm = - L / 2d0
        dx = L / (np - 1)
        p = hbar ** 2 / (2 * m)

        call idn(np, -1, idm)
        call idn(np, 0, idz)
        call idn(np, +1, idp)

        am = (idm - 2 * idz + idp) / dx ** 2
        bm = (idm + 10 * idz + idp) / 12
        vm = 0

        call inm(bm, bi, np)

        do i = 1, np
          xr(i) = xm + (i - 1) * dx
          v(i) = pot(xr(i), k)
          vm(i, i) = v(i)
        end do

        ! H = - p * am + vm
        H = - p * matmul(bi, am) + vm
        a = H

c       Diagonalization & Computing Eigenvalues, Eigenvectors
        call tred2(a, np, np, d, e)
        call tql2(np, np, d, e, a, ierr)
        do i = 1, ns
71        format (A,1X,I3,1X,A,1X,G23.16)
          print 71, "Eigenvalue", i, "=", d(i)
        end do

c       Normalization of the Ground State using Composite Simpson Integral
        sum = a(1, 1) ** 2 + a(np, 1) ** 2
        do i = 2, np - 1
          if (i .eq. (i / 2) * 2) then
            sum = sum + 2 * a(i, 1) ** 2
          else
            sum = sum + 4 * a(i, 1) ** 2
          end if
        end do
        csint = dx / 3 * sum
        normc = sqrt(dble(1d0 / csint)) ! Normalization Constant
        a = normc * a
        print '(/,A,1X,A,1X,G23.16)', "Normalization Constant",
     . "of the ground state wavefunction =", normc

c ----- GNUPlot Plots ------------------------------------------------
! 75      format(A,E20.10,A,E20.10,A)

!         open(10, file='LHO.dat')
! c       Data
!         write(10, '(11X,A,24X,A,17X,999(A,F11.8,9X))') "X", "V",
!      .    ("# E =", d(i), i = 1, ns)
!         do i = 1, np
!           write(10, '(999(E20.10,5X))') xr(i), v(i),
!      .      (a(i, j), j = 1, ns)
!         end do
!         write(10, '(999(E20.10,5X))') xr(1), v(1), (a(1, j), j = 1, ns)

! c       Plots
!         open(11, file='LHO.plt')
!         write(11, '(A,/)') "set term post enhanced color 'Monospace' 12"

! cc      First 10 Eigenstates
!         write(11, '(A)') "set output 'LHO.ps'"
!         write(11, '(A, 1X, A)') "set title",
!      .    "'Linear Harmonic Oscillator'"
!         write(11, '(A)') "set style data lines"
!         write(11, '(A)') "set xlabel 'x'"
!         write(11, '(A)') "set ylabel '{/Symbol y}(x)'"
!         ! Spacing of energy levels
!         de = (d(ns) - d(1)) / dble(ns - 1) * 2d0 / normc
!         write(11, 75) "set xrange [", xr(1), ":", xr(np), "]"
!         write(11, 75) "set yrange [", minval(v), ":", d(ns)+de, "]"
!         write(11, '(A,A)') "plot 'LHO.dat' u 1: 2 t 'Potential', ",
!      .    char(92)
!         do i = 1, ns - 1
!           write(11, '(A,E20.10,A,E20.10,A,I3,A,I3,A,A)') "'' u 1:(",
!      .    d(i), "+", de, "*$", i + 2, ") t 'Eigenstate ", i,
!      .    "',", char(92)
!         end do
!         write(11, '(A,E20.10,A,E20.10,A,I3,A,I3,A)') "'' u 1:(",
!      .    d(ns), "+", de, "*$", ns + 2, ") t 'Eigenstate ", ns, "'"

! cc      Ground State
!         write(11, '(/,A)') "set output 'LHOG.ps'"
!         write(11, '(A, 1X, A)') "set title 'Ground State of",
!      .    "Linear Harmonic Oscillator'"
!         write(11, '(A)') "set style data lines"
!         write(11, '(A)') "set xlabel 'x'"
!         write(11, '(A)') "set ylabel '{/Symbol y}(x)'"
!         de = 1.0 ! Average Spacing
!         write(11, 75) "set xrange [", xr(1), ":", xr(np), "]"
!         write(11, 75) "set yrange [", minval(v), ":", d(1) + de, "]"
!         write(11, '(A,A)') "plot 'LHO.dat' u 1:2 t 'Potential', ",
!      .    char(92)
!         write(11, 75) "'' u 1:(", d(1), "+",
!      .    de, "*$3) t 'Ground State (Normalized)"

!         call execute_command_line("gnuplot LHO.plt")
!         call execute_command_line("ps2pdf LHO.ps && rm LHO.ps")
!         call execute_command_line("ps2pdf LHOG.ps && rm LHOG.ps")
      end

      subroutine idn(n, k, idnk)
c     nxn matrix that has 1s along kth diagonal
        integer n, k
        real idnk(n, n)
        idnk = 0.0
        do i = 1, n
          idnk(i, i - k) = 1.0
        end do
        return
      end

      function pot(x, k)
c     The Potential
        real x, k, pot
        pot = 0.5 * k * x ** 2
        return
      end

c inverse of a matrix 'a' is 'c'
      subroutine inm(a,c,n)
      real a(n,n), c(n,n), L(n,n), U(n,n), b(n), d(n), x(n)
      L=0.0
      U=0.0
      b=0.0
      do k=1, n-1
      do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
      a(i,j) = a(i,j)-coeff*a(k,j)
      end do
      end do
      end do
      do i=1,n
      L(i,i) = 1.0
      end do
      do j=1,n
      do i=1,j
      U(i,j) = a(i,j)
      end do
      end do
      do k=1,n
      b(k)=1.0
      d(1) = b(1)
      do i=2,n
      d(i)=b(i)
      do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
      end do
      end do
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
      x(i) = d(i)
      do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
      end do
      x(i) = x(i)/u(i,i)
      end do
      do i=1,n
      c(i,k) = x(i)
      end do
      b(k)=0.0
      end do
      end

c converts real symmetric matrix to tridiagonal matrix
      subroutine tred2(a,n,np,d,e)
      dimension a(np,np),d(np),e(np)
      if(n.gt.1)then
        do 18 i=n,2,-1  
          l=i-1
          h=0.
          scale=0.
          if(l.gt.1)then
            do 11 k=1,l
              scale=scale+abs(a(i,k))
11          continue
            if(scale.eq.0.)then
              e(i)=a(i,l)
            else
              do 12 k=1,l
                a(i,k)=a(i,k)/scale
                h=h+a(i,k)**2
12            continue
              f=a(i,l)
              g=-sign(sqrt(h),f)
              e(i)=scale*g
              h=h-f*g
              a(i,l)=f-g
              f=0.
              do 15 j=1,l
                a(j,i)=a(i,j)/h
                g=0.
                do 13 k=1,j
                  g=g+a(j,k)*a(i,k)
13              continue
                if(l.gt.j)then
                  do 14 k=j+1,l
                    g=g+a(k,j)*a(i,k)
14                continue
                endif
                e(j)=g/h
                f=f+e(j)*a(i,j)
15            continue
              hh=f/(h+h)
              do 17 j=1,l
                f=a(i,j)
                g=e(j)-hh*f
                e(j)=g
                do 16 k=1,j
                  a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
16              continue
17            continue
            endif
          else
            e(i)=a(i,l)
          endif
          d(i)=h
18      continue
      endif
      d(1)=0.
      e(1)=0.
      do 23 i=1,n
        l=i-1
        if(d(i).ne.0.)then
          do 21 j=1,l
            g=0.
            do 19 k=1,l
              g=g+a(i,k)*a(k,j)
19          continue
            do 20 k=1,l
              a(k,j)=a(k,j)-g*a(k,i)
20          continue
21        continue
        endif
        d(i)=a(i,i)
        a(i,i)=1.
        if(l.ge.1)then
          do 22 j=1,l
            a(i,j)=0.
            a(j,i)=0.
22        continue
        endif
23    continue
      return
      end

c finds the eigenvalues and eigenvectors of a tridiagonal matrix
      subroutine tql2(nm,n,d,e,z,ierr)
      real d(n),e(n),z(nm,n)
      ierr = 0
      if (n .eq. 1) go to 1001
      do 100 i = 2, n
  100 e(i-1) = e(i)
      f = 0.0e0
      tst1 = 0.0e0
      e(n) = 0.0e0
      do 240 l = 1, n
        j = 0
        h = abs(d(l)) + abs(e(l))
        if (tst1 .lt. h) tst1 = h
        do 110 m = l, n
          tst2 = tst1 + abs(e(m))
          if (tst2 .eq. tst1) go to 120
  110   continue
  120   if (m .eq. l) go to 220
  130   if (j .eq. 30) go to 1000
        j = j + 1
        l1 = l + 1
        l2 = l1 + 1
        g = d(l)
        p = (d(l1) - g) / (2.0e0 * e(l))
        r = pythag(p,1.0e0)
        d(l) = e(l) / (p + sign(r,p))
        d(l1) = e(l) * (p + sign(r,p))
        dl1 = d(l1)
        h = g - d(l)
        if (l2 .gt. n) go to 145
        do 140 i = l2, n
  140   d(i) = d(i) - h
  145   f = f + h
        p = d(m)
        c = 1.0e0
        c2 = c
        el1 = e(l1)
        s = 0.0e0
        mml = m - l
        do 200 ii = 1, mml
          c3 = c2
          c2 = c
          s2 = s
          i = m - ii
          g = c * e(i)
          h = c * p
          r = pythag(p,e(i))
          e(i+1) = s * r
          s = e(i) / r
          c = p / r
          p = c * d(i) - s * g
          d(i+1) = h + s * (c * g + s * d(i))
          do 180 k = 1, n
            h = z(k,i+1)
            z(k,i+1) = s * z(k,i) + c * h
            z(k,i) = c * z(k,i) - s * h
  180     continue
  200   continue
        p = -s * s2 * c3 * el1 * e(l) / dl1
        e(l) = s * p
        d(l) = c * p
        tst2 = tst1 + abs(e(l))
        if (tst2 .gt. tst1) go to 130
  220   d(l) = d(l) + f
  240 continue
      do 300 ii = 2, n
        i = ii - 1
        k = i
        p = d(i)
        do 260 j = ii, n
          if (d(j) .ge. p) go to 260
          k = j
          p = d(j)
  260   continue
        if (k .eq. i) go to 300
        d(k) = d(i)
        d(i) = p
        do 280 j = 1, n
          p = z(j,i)
          z(j,i) = z(j,k)
          z(j,k) = p
  280   continue
  300 continue
      go to 1001
 1000 ierr = l
 1001 return
      end

c finds sqrt(a**2+b**2)
      real function pythag(a,b)
      p = amax1(abs(a),abs(b))
      if (p .eq. 0.0e0) go to 20
      r = (amin1(abs(a),abs(b))/p)**2
   10 continue
         t = 4.0e0 + r
         if (t .eq. 4.0e0) go to 20
         s = r/t
         u = 1.0e0 + 2.0e0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
