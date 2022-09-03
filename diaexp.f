      program diaexp
      real h,g,s,p,r,hh,pi,Le,q,x0
      integer n,i,j,k,l,m,iter
      parameter(np=100,tiny=1.0e-38)
      dimension  X(np),u(np),a(np,np)
      dimension c(np,np),d(np),e(np),f(np),z(np,np)

c	specified potential v(x)=sin(pi*x/L),L=2*a,-a<=x<=a,a=1

      open(1,file="gsi")
      pi=4*atan(1.0)
      Le=10.0
      delta = (Le/(np-1))
      q=1.0/delta**2    !2.41e-37
      x0=0 !-a
c	getting potential array
      do 100 i=1,np
      ! X(i)=i*(Le/np)
      ! u(i)=sin(pi*X(i)/Le)
      X(i)=x0+(i-1)*delta
      u(i) = sin(pi*((i-1)*delta)/Le)
      ! u(i) = 0.0
      d(i) = 2*q+u(i) ! <---->
      e(i) = - q
!	write(*,*)u(i)
100   continue

c	making all the elements of H to be 0	
      do 98 i=1,np
      do 99 j=1,np
      a(i,j)=0.0
99    continue
98    continue
        
c giving input for the matrix H

      do 101 i=1,np
      do 102 j=1,np
      if(j==i)then
      a(i,j)=2*q + u(i)
      endif
      if(j/=i)then
      if(i>1)then
      a(i,i-1)=-q
      a(i,i+1)=-q
      else
      a(i,i+1)=-q
      end if
      end if
102   continue
101   continue

c	printing elements of H
      do 103 i=1,np
      ! write(*,*)(a(i,j),j=1,np)
103   continue
	
c	starting diagonalization	
	     do 12 i=1,np
        do 11 j=1,np
         c(i,j)=a(i,j)
11      continue
12    continue
      ! call tred2(a,np,np,d,e)
      ! call tqli(d,e,np,np,a)
      call tql2(np,np,d,e,a,ierr)
      ! write (*,*) 'eigenvectors for a real symmetric matrix'
      do 16 i=1,np
        do 14 j=1,np
          f(j)=0.0
          do 13 k=1,np
            f(j)=f(j)+c(j,k)*a(k,i)
13        continue
14      continue
        ! write (*,*) 'eigenvalue',i,' =',d(i)
        ! write (*,*) ' vector','            mtrx*vect.','       ratio'
        do 15 j=1,np
          if (abs(a(j,i)).lt.tiny) then
            ! write (*,*) a(j,i),f(j),'div. by 0'
          else
            ! write (*,*) a(j,i),f(j),f(j)/a(j,i)
          endif
15      continue
        ! write (*,*) ''
16    continue
c writing the ground state wave function and x values

      do 505 i=1,np
      write(1,*)X(i),a(i,1)
505   continue
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
