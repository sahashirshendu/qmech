      program diag
      parameter(np=10,tiny=1.0e-6)
      dimension a(np,np),c(np,np),d(np),e(np),f(np)
      data c/5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0,-4.0,
     *            4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0,
     *            3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,
     *            2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,
     *            1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,
     *            0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,
     *            -1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,
     *            -2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,
     *            -3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,
     *            -4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0/
      do 12 i=1,np
        do 11 j=1,np
          a(i,j)=c(i,j)
11      continue
12    continue
      call tred2(a,np,np,d,e)
      call tqli(d,e,np,np,a)
      write (*,*) 'eigenvectors for a real symmetric matrix'
      do 16 i=1,np
        do 14 j=1,np
          f(j)=0.0
          do 13 k=1,np
            f(j)=f(j)+c(j,k)*a(k,i)
13        continue
14      continue
        write (*,*) 'eigenvalue',i,' =',d(i)
        write (*,*) ' vector','            mtrx*vect.','       ratio'
        do 15 j=1,np
          if (abs(a(j,i)).lt.tiny) then
            write (*,*) a(j,i),f(j),'div. by 0'
          else
            write (*,*) a(j,i),f(j),f(j)/a(j,i)
          endif
15      continue
        write (*,*) ''
16    continue
      end

c TRED2 reduces the input matrix to tridiagonal form
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
            if(iter.eq.30)pause 'too many iterations'
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
