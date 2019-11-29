    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER :: noel,npt,done
    REAL(dp) :: time,strain,stress(4),pretime,estrain,pstrain
    Real(dp) :: prinelas,prinplas,prin,miseselas,misesplas
    real(dp) :: wtime,wstrain,wstress,E
    open(10,file='fort.91')  
read(10,*) pretime,done,stress(1),stress(2),stress(3),stress(4)
Do
6 read(10,*,end=5) time,done,E,stress(1),stress(2),stress(3),stress(4)
  if (time.eq.pretime) then
   goto 6
  endif
  if (time.gt.pretime) then
   write(22,fmt='(f8.4,1x,i3,5(es12.5,1x))') time,done,E,stress(1),&
   stress(2),stress(3),stress(4)
   pretime=time
   goto 6
  endif
end do
5 close(10)
  close(22)
end

 
 
 

