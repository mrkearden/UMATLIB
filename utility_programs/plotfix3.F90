    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER :: noel,npt,done,wdone,itens
    REAL(dp) :: time,strain,stress(4),pretime,estrain,pstrain
    Real(dp) :: prinelas,prinplas,prin,miseselas,misesplas
    real(dp) :: wtime,wstrain,wstress(4),E,WE
    character(6) :: fname1
    character(1) :: fname2
    character(7) :: fname3
    fname1="fort.9"
    write(*,*) " Which stress tensor? "
    read(*,fmt='(a1)') fname2
    fname3=fname1 // fname2
    open(10,file=fname3)  
read(10,*) pretime,done,stress(1),stress(2),stress(3),stress(4)
Do
6 read(10,*,end=5) time,done,E,stress(1),stress(2),stress(3),stress(4)
  if (time.eq.pretime) then
   wtime=time
   wdone=done
   WE=E
   wstress = stress
   goto 6
  endif
  if (time.gt.pretime) then
   write(22,fmt='(f8.4,1x,i3,5(es12.5,1x))') wtime,wdone,WE,wstress(1),&
   wstress(2),wstress(3),wstress(4)
   pretime=time
   goto 6
  endif
end do
5 write(22,fmt='(f8.4,1x,i3,5(es12.5,1x))') time,done,E,stress(1),&
   stress(2),stress(3),stress(4)
  close(10)
  close(22)
end

 
 
 

