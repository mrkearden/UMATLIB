    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER :: noel,npt,done,wnoel,wnpt
    REAL(dp) :: time,strain,stress,pretime,estrain,pstrain,westrain,wpstrain
    Real(dp) :: prinelas,prinplas,prin,miseselas,misesplas
    Real(dp) :: wprinelas,wprinplas,wprin,wmiseselas,wmisesplas
    real(dp) :: wtime,wstrain,wstress,E
    open(10,file='fort.90')  
read(10,*) pretime,noel,npt,estrain,pstrain,strain,prinelas,prinplas,prin,&
miseselas,misesplas,stress
Do
6 read(10,*,end=5) time,noel,npt,estrain,pstrain,strain,prinelas,prinplas,prin,&
miseselas,misesplas,stress
  if (time.eq.pretime) then
   wtime=time
   wnoel=noel
   wnpt=npt
   westrain=estrain
   wpstrain=pstrain
   wstrain=strain
   wprinelas=prinelas
   wprinplas=prinplas
   wprin=prin
   wmiseselas=miseselas
   wmisesplas=misesplas
   wstress=stress
   goto 6
  endif
  if (time.gt.pretime) then
   write(21,fmt='(f8.4,2(1x,i3),9(es12.5,1x))') wtime,wnoel,wnpt,westrain,wpstrain,wstrain,wprinelas,wprinplas,wprin,&
wmiseselas,wmisesplas,wstress
   pretime=time
   goto 6
  endif
end do
5 write(21,fmt='(f8.4,2(1x,i3),9(es12.5,1x))') time,noel,npt,estrain,pstrain,strain,prinelas,prinplas,prin,&
miseselas,misesplas,stress
  close(10)
  close(21)
end

 
 
 

