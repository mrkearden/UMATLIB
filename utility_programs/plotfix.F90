    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER :: i
    REAL(dp) :: time,strain,stress,pretime
    real(dp) :: wtime,wstrain,wstress
    open(10,file='fort.89')  
read(10,*) pretime,strain,stress
Do
6 read(10,*,end=5) time,strain,stress
  if (time.eq.pretime) then
   wtime=time
   wstrain=strain
   wstress=stress
   goto 6
  endif
  if (time.gt.pretime) then
   write(20,fmt='(f8.4,es12.4,es12.5)') wtime,wstrain,wstress
   pretime=time
   goto 6
  endif
end do
5 close(10)
  close(20)
end

 
 
 

