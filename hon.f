      program hon
      integer,parameter::n=100
      real,dimension(n)::x,psi
      real f,L
      common x,psi,h
      external f

      open(1,file='hon.txt')
      L = 8
      e1 = 0
      e2 = 1
      h = L/(n-1)
      do 11 i = 1,n
      x(i) = -.5*L+(i-1)*h
11    continue
      do 12 i = 1,100
      e = (e1+e2)/2
      if ((f(e)*f(e2)).le.0) then
      e1 = e
      else
      e2 = e
      end if
      if (abs((e2-e1)/e2).lt.0.000001) exit
12    continue
      write(*,*) 'Energy =',e
      do i = 1,n
      write(1,*) x(i),psi(i)
      end do
      stop
      end

      function f(e)
      integer,parameter::n=100
      real,dimension(n)::x,psi,k2
      common x,psi,h
      do 13 i = 1,n
        k2(i) = 2.*(e-.5*x(i)**2)
13    continue
      psi(1) = 0
      psi(2) = 1e-5
      do 14 i = 2,n-1,1
        psi(i+1) = ((2.-5.*h**2/6.*k2(i))*psi(i)
     & -(1.+h**2/12.*k2(i-1))*psi(i-1))/(1.+h**2/12.*k2(i+1))
14    continue
      f = psi(n)
      return
      end
