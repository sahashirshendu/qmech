      program seod
      parameter(n=10)
      real a(n,n),h(n,n),d(n),e(n),f(n)
      call hmg(h,d,e,n)
      write(*,*) "Hamiltonian -"
      do 11 i = 1,n
        write(*,*) (h(i,j), j=1,n)
11    continue
      a = 0.
      do 12 i=1,n
        a(i,i)=1.
12    continue
      call tqli(d,e,n,n,a)
      do 16 i=1,n
        do 14 j=1,n
          f(j)=0.0
          do 13 k=1,n
            f(j)=f(j)+h(j,k)*a(k,i)
13        continue
14      continue
        write(*,*) "Eigenvalue", i, " =", d(i)
        write(*,*) "    Vector", "         Mtx * Vec", "          Ratio"
        do 15 j=1,n
          write(*,*) a(j,i),f(j),f(j)/a(j,i)
15      continue
        write(*,*) ""
16    continue
      end

      subroutine hmg(h,d,e,n)
      real L,h(n,n),d(n),e(n),v(n)
      x0 = 0.
      L = 10.
      a = L/(n - 1)
      pi = 4.*atan(1.)
      do 11 i = 1,n
        v(i) = sin(pi * (x0 + (i-1)*a) / L)
        d(i) = v(i) + 2./a**2
        e(i) = - 1./a**2
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

c TQLI finds the eigenvalues and eigenvectors of a tridiagonal matrix
      subroutine tqli(d,e,n,np,z)
      dimension d(np),e(np),z(np,np)
      if (n.gt.1) then
        do 11 i=2,n
          e(i-1)=e(i)
11      continue
        e(n)=0.
        do 15 l=1,n
          iter=0
1         do 12 m=l,n-1
            dd=abs(d(m))+abs(d(m+1))
            if (abs(e(m))+dd.eq.dd) go to 2
12        continue
          m=n
2         if(m.ne.l)then
            if(iter.eq.30)stop 'too many iterations'
            iter=iter+1
            g=(d(l+1)-d(l))/(2.*e(l))
            r=sqrt(g**2+1.)
            g=d(m)-d(l)+e(l)/(g+sign(r,g))
            s=1.
            c=1.
            p=0.
            do 14 i=m-1,l,-1
              f=s*e(i)
              b=c*e(i)
              if(abs(f).ge.abs(g))then
                c=g/f
                r=sqrt(c**2+1.)
                e(i+1)=f*r
                s=1./r
                c=c*s
              else
                s=f/g
                r=sqrt(s**2+1.)
                e(i+1)=g*r
                c=1./r  
                s=s*c
              endif
              g=d(i+1)-p
              r=(d(i)-g)*s+2.*c*b
              p=s*r
              d(i+1)=g+p
              g=c*r-b
              do 13 k=1,n
                f=z(k,i+1)
                z(k,i+1)=s*z(k,i)+c*f
                z(k,i)=c*z(k,i)-s*f
13            continue
14          continue
            d(l)=d(l)-p
            e(l)=g
            e(m)=0.
            go to 1
          endif
15      continue
      endif
      return
      end
