    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER :: noel,npt,done
    REAL(dp) :: time,strain,stress,pretime,estrain,pstrain
    Real(dp) :: prinelas,prinplas,prin,miseselas,misesplas
    real(dp) :: wtime,wstrain,wstress,E
    open(10,file='fort.90')  
read(10,*) pretime,noel,npt,estrain,pstrain,strain,prinelas,prinplas,prin,&
miseselas,misesplas,stress
Do
6 read(10,*,end=5) time,noel,npt,estrain,pstrain,strain,prinelas,prinplas,prin,&
miseselas,misesplas,stress
  if (time.eq.pretime) then
   goto 6
  endif
  if (time.gt.pretime) then
   write(21,fmt='(f8.4,2(1x,i3),9(es12.5,1x))') time,noel,npt,estrain,pstrain,strain,prinelas,prinplas,prin,&
miseselas,misesplas,stress
   pretime=time
   goto 6
  endif
end do
5 close(10)
  close(21)
end

 
 
 

