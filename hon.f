      program hon
      parameter(n=99)
      real a(n,n),h(n,n),x(n),d(n),e(n),sum,int
      call hmg(h,x,n,dx)
      do 12 i=1,n
        do 13 j=1,n
          a(i,j)=h(i,j)
13      continue
12    continue
      call tred2(a,n,n,d,e)
      call tql2(n,n,d,e,a,ierr)
      do 16 i=1,n
        do 14 j=1,n
14      continue
        write(*,*) "Eigenvalue", i, " =", d(i)
16    continue

      sum = 0
      do i=2,n-1
        sum = sum+2*a(i,1)**2
      enddo
      int = dx/2. * (a(1,1)**2 + sum + a(n,1)**2)
      cons = sqrt(1.0/int)
      a = cons * a
      write(*,*) "Normalization Constant =", cons

      open(11,file='hon.txt')
      do 18 i = 1, n
        write(11,*) x(i), a(i,1)+d(1), 0.5*x(i)**2
18    continue

      open(12,file='hona.txt')
      do 19 i = 1, n
        write(12,*) x(i), 0.75*a(i,1)+d(1), 0.75*a(i,2)+d(2),
     &  0.5*x(i)**2
19    continue
      end

      subroutine hmg(h,x,n,a)
      real L,f(n,n),g(n,n),gi(n,n),h(n,n),v(n,n),x(n),a
      x0 = -5.
      L = 10.
      a = L/(n - 1)
      v = 0.
      f = 0.
      g = 0.
      do 11 i = 1,n
        x(i) = x0+(i-1)*a
        v(i,i) = 0.5*1.*(x0+(i-1)*a)**2
11    continue
      do 12 i = 1,n
        do 13 j = 1,n
          if (i.eq.j) then
            f(i,j) = -2.
            g(i,j) = 10.
            if (j.gt.1) then
              f(i,j-1) = 1.
              g(i,j-1) = 1.
            end if
            if (j.lt.n) then
              f(i,j+1) = 1.
              g(i,j+1) = 1.
            end if
          end if
  13    continue
12    continue
      f=f/(a**2)
      g=g/12.
      call inm(g,gi,n)
      h = -1./2.*matmul(gi,f) + v
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

