C234567
C matrix A
      dimension a(6,6),c(6,6)
      n=6
      a(1,1)=3.333E-8
      a(1,2)=-1.E-8
      a(1,3)=-1.E-8
      a(1,4)=0.0
      a(1,5)=0.0
      a(1,6)=0.0
      a(2,1)=-1.E-8
      a(2,2)=3.333E-8
      a(2,3)=-1.E-8
      a(2,4)=0.0
      a(2,5)=0.0
      a(2,6)=0.0
      a(3,1)=-1.E-8
      a(3,2)=-1.e-8
      a(3,3)=3.333E-8
      a(3,4)=0.0
      a(3,5)=0.0
      a(3,6)=0.0
      a(4,1)=0.0
      a(4,2)=0.0
      a(4,3)=0.0
      a(4,4)=9.09E-6
      a(4,5)=0.0
      a(4,6)=0.0
      a(5,1)=0.0
      a(5,2)=0.0
      a(5,3)=0.0
      a(5,4)=0.0
      a(5,5)=9.09E-6
      a(5,6)=0.0
      a(6,1)=0.0
      a(6,2)=0.0
      a(6,3)=0.0
      a(6,4)=0.0
      a(6,5)=0.0
      a(6,6)=9.09E-6
      write (*,200)
      do 10 j=1,n
   10 write (*,201) a(j,1),a(j,2),a(j,3),a(j,4),a(j,5),a(j,6)
      call inverse(a,c,n)
      write (*,202)
      do 20 j=1,n
   20 write (*,201) c(j,1),c(j,2),c(j,3),c(j,4),c(j,5),c(j,6)
      sigmax=c(1,1)*.003+c(1,2)*.001+c(1,3)*.001
      sigmay=c(2,1)*.001+c(2,2)*.001+c(2,3)*.001
      sigmaz=c(3,1)*.001+c(3,2)*.001+c(3,3)*.001
      sigmayz=c(4,4)*.000006
      sigmazx=c(5,5)*.000006
      sigmaxy=c(6,6)*.0006
      write(*,*) ' SIGMx:',sigmax
      write(*,*) ' SIGMy:',sigmay
      write(*,*) ' SIGMz:',sigmaz
      write(*,*) ' SIGMyz:',sigmayz
      write(*,*) ' SIGMzx:',sigmazx
      write(*,*) ' SIGMxy:',sigmaxy
  200 format (' Computing Inverse matrix ',/,/,' Matrix A')
  201 format (6(1x,e12.6,2x))
  202 format (/,' Inverse matrix A^{-1}')
      end
      subroutine inverse(a,c,n)
      dimension xL(6,60),b(6),a(6,6),c(6,6),U(6,6),d(6),x(6)
C23456
      do 11 i=1,n
      do 10 j=1,n
      xL(i,j)=0.0
      U(i,j)=0.0
   10 b(j)=0.0
   11 continue
      do 20 k=1, n-1
      do 30 i=k+1,n
      coeff=a(i,k)/a(k,k)
      xL(i,k) = coeff
      do 40 j=k+1,n
      a(i,j) = a(i,j)-coeff*a(k,j)
   40 continue
   30 continue
   20 continue
      do 50 i=1,n
      xL(i,i) = 1.0
   50 continue
      do 60 j=1,n
      do 70 i=1,j
      U(i,j) = a(i,j)
   70 continue
   60 continue
      do 80 k=1,n
      b(k)=1.0
      d(1) = b(1)
      do 90 i=2,n
      d(i)=b(i)
      do 91 j=1,i-1
      d(i) = d(i) - xL(i,j)*d(j)
   91 continue
   90 continue
      x(n)=d(n)/U(n,n)
      do 92 i = n-1,1,-1
      x(i) = d(i)
      do 93 j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
   93 continue
      x(i) = x(i)/u(i,i)
   92 continue
      do 94 i=1,n
      c(i,k) = x(i)
   94 continue
      b(k)=0.0
   80 continue
      end subroutine inverse
