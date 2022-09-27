      program ho
      parameter(n=99)
      real a(n,n),h(n,n),x(n),d(n),e(n),f(n),sum,int
      call hmg(h,x,d,e,n,dx)
      a = 0.
      do 12 i=1,n
        a(i,i)=1.
12    continue
      call tql2(n,n,d,e,a,ierr)
      do 16 i=1,n
        do 14 j=1,n
          f(j)=0.0
          do 13 k=1,n
            f(j)=f(j)+h(j,k)*a(k,i)
13        continue
14      continue
        write(*,*) "Eigenvalue", i, " =", d(i)
16    continue

      sum = 0
      do i=2,n-1
        sum = sum + 2*a(i,1)**2
      enddo
      int = dx/2. * (a(1,1)**2 + sum + a(n,1)**2)
      cons = sqrt(1.0/int)
      a = cons * a
      write(*,*) "Normalization Constant =", cons

      open(11,file='ho.txt')
      do 18 i = 1, n
        write(11,*) x(i), a(i,1)+d(1), 0.5*x(i)**2
18    continue
      end

      subroutine hmg(h,x,d,e,n,a)
      real L,h(n,n),d(n),e(n),v(n),x(n),a
      x0 = -5.
      L = 10.
      a = L/(n - 1)
      pi = 4.*atan(1.)
      do 11 i = 1,n
        x(i) = x0+(i-1)*a
        v(i) = 0.5*1.*(x0+(i-1)*a)**2
        d(i) = v(i) + 1./a**2
        e(i) = - 1./(2*a**2)
11    continue
      h = 0.
      do 12 i = 1,n
        do 13 j = 1,n
          if (i.eq.j) then
            h(i,j) = d(i)
            if (j.gt.1) then
              h(i,j-1) = e(i)
            end if
            if (j.lt.n) then
              h(i,j+1) = e(i)
            end if
          end if
  13    continue
12    continue
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
