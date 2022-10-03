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
