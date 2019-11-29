!------------------------------------------------------------------------------
  SUBROUTINE bi_linear(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
!------------------------------------------------------------------------------
!
! Bi-linear plastic material model for 2D elements 
! 7 Props input nu,strain at yield,yield stress, second strain point,stress,
! write data option, element data to write (-1 = all)
!
!------------------------------------------------------------------------------
    USE Types
    IMPLICIT NONE

    REAL(KIND=dp), INTENT(INOUT) :: STRESS(NTENS)
    ! Requirement for Elmer: At the time of calling the Cauchy stress T_n before
    ! the time/load increment is given
    ! Requirement for umat:  The stress T_{n+1}^{(k)} corresponding to the 
    ! current approximation of the strain increment (DSTRAN) must be returned. 
    ! If the strain increment is defined to be zero in the beginning of the
    ! nonlinear iteration, Elmer will generate a candidate for the strain increment
    ! by assuming purely elastic increment characterized by DDSDDE.

    REAL(KIND=dp), INTENT(INOUT) :: STATEV(NSTATEV)
    ! Requirement for Elmer: The state variables Q_n as specified at the 
    ! previous time/load level for converged solution are given.
    ! Requirement for umat:  The state variables Q_{n+1}^{(k)} corresponding to 
    ! the current approximation of the strain increment must be returned. If 
    ! convergence is attained, these values will be saved and associated with the 
    ! converged solution (cf. the input values)

    REAL(KIND=dp), INTENT(OUT) :: DDSDDE(NTENS,NTENS)
    ! The derivative of (Cauchy) stress response function with respect to the 
    ! strain evaluated for the current approximation must be returned

    REAL(KIND=dp), INTENT(INOUT) :: SSE, SPD, SCD
    ! Requirement for Elmer: Provide specific strain energy (sse), plastic 
    ! dissipation (spd) and creep dissipation (scd) at the previous time/load 
    ! level (these are supposed to be declared to be state variables)
    ! Requirement for umat:  The values of the energy variables corresponding to 
    ! the current approximation may be returned

    REAL(KIND=dp), INTENT(OUT) :: rpl
    ! The mechanical heating power (volumetric)

    REAL(KIND=dp), INTENT(OUT) :: ddsddt(NTENS), drplde(NTENS), drpldt

    REAL(KIND=dp), INTENT(IN) :: STRAN(NTENS)
    ! This gives the strains before the time/load increment.
    ! The strain can be computed from the deformation gradient, so this
    ! argument can be considered to be redundant. Elmer provides
    ! this information anyway. Abaqus assumes that the logarithmic strain 
    ! is used, but Elmer may also use other strain measures.

    REAL(KIND=dp), INTENT(IN) :: DSTRAN(NTENS)
    ! The current candidate for the strain increment to obtain the current 
    ! candidate for the stress. In principle this could be computed from the 
    ! deformation gradient; cf. the variable stran.

    REAL(KIND=dp), INTENT(IN) :: TIME(2)
    ! Both entries give time before the time/load increment (the time for the last
    ! converged solution

    REAL(KIND=dp), INTENT(IN) :: DTIME
    ! The time increment

    REAL(KIND=dp), INTENT(IN) :: TEMP
    ! Temperature before the time/load increment

    REAL(KIND=dp), INTENT(IN) :: dtemp
    ! Temperature increment associated wíth the time/load increment. Currently
    ! Elmer assumes isothermal conditions during the load increment.

    REAL(KIND=dp), INTENT(IN) :: predef(1), dpred(1)
    ! These are just dummy variables for Elmer

    CHARACTER(len=80), INTENT(IN) :: CMNAME
    ! The material model name

    INTEGER, INTENT(IN) :: NDI
    ! The number of direct stress components

    INTEGER, INTENT(IN) :: NSHR
    ! The number of the engineering shear strain components

    INTEGER, INTENT(IN) :: NTENS 
    ! The size of the array containing the stress or strain components

    INTEGER, INTENT(IN) :: NSTATEV
    ! The number of state variables associated with the material model

    REAL(KIND=dp), INTENT(IN) :: PROPS(NPROPS)
    ! An array of material constants

    INTEGER, INTENT(IN) :: NPROPS
    ! The number of the material constants

    REAL(KIND=dp), INTENT(IN) :: coords(3)
    ! The coordinates of the current point could be specified

    REAL(KIND=dp), INTENT(IN) :: drot(3,3)
    ! No support for keeping track of rigid body rotations 
    ! (the variable is initialized to the identity)

    REAL(KIND=dp), INTENT(INOUT) :: pnewdt
    ! Currently, suggesting a new size of time increment does not make any impact

    REAL(KIND=dp), INTENT(IN) :: celent
    ! The element size is not yet provided by Elmer

    REAL(KIND=dp), INTENT(IN) :: DFRGRD0(3,3)
    ! The deformation gradient before the time/load increment (at the previous 
    ! time/load level for converged solution)

    REAL(KIND=dp), INTENT(IN) :: DFRGRD1(3,3)
    ! The deformation gradient corresponding to the current approximation
    ! (cf. the return value of STRESS variable) 

    INTEGER, INTENT(IN) :: NOEL
    ! The element number

    INTEGER, INTENT(IN) :: NPT
    ! The integration point number

    INTEGER, INTENT(IN) :: layer, kspt, kstep, kinc
    ! kstep and kinc could be provided to give information on the incrementation
    ! procedure
!
!  Treating as 2D Plane stress element Stress XX, YY, XY
!  
!  Bi-linear stress-strain curve
!
!  Material Constants 7
!  1) Density, 2) strain at yield, 3) yield stress, 4)strain at final, 5) stress at final
!  6) Write output file 1.0=yes, 7) element number to write (-1.0) is all and all npt)
!  8) NPT to write (-1.0) for all
!------------------------------------------------------------------------------
    ! Local variables:
    INTEGER :: i,iwrite,ielw,done(ntens),iunit,j,istatdone, inpt
    REAL(KIND=dp) :: nu, E, E2, LambdaLame, MuLame
    LOGICAL :: exists
    REAl(KIND=dp) :: mises, term1, term2, term3, term4, elaeprinmax, elaeprinmin
    REAL(KIND=dp) :: plaeprinmax,plaeprinmin,plasprinmax,plasprinmin, elamises, plamises
    REAL(KIND=dp) :: elastrain(ntens),plastrain(ntens),elastress(ntens),plastress(ntens),tresca
    real(kind=dp) :: edprinmax,edprinmin,elasprinmax,elasprinmin,tempstran(ntens)
    real(kind=dp) :: strain1,stressc1,strain2,stressc2,term,eprint(ntens)
    REAL(KIND=dp) :: totstrain(ntens),totstress(ntens)
    REAL(KIND=dp) :: dstranelas(NTENS),dstranplas(ntens)
!------------------------------------------------------------------------------

    ! Get Young's modulus and the Poisson ratio:
    E = Props(3)/Props(2)
    E2 = Props(5)/Props(4)
    nu = Props(1)
    iwrite = int(Props(6))
    ielw = int(Props(7))
    inpt = int(Props(8))
    strain1 = Props(2)
    stressc1 = Props(3)
    strain2 = Props(4)
    stressc2 = Props(5)
    ddsdde = 0.0d0
!
!  Could be on first slope done = 1
!  Could be on second slope done = 2
!  DSTRAN could be spanning the change in slope done = 3
!  could be negative done = 0
!  start ntens do loop
!
!  for each ntens calculate ddsdde
!
do i=1,ntens
!
! reset variables for next ntens
!
!
! If dstran is negative not doing unloading or compression
!
    E = Props(3)/Props(2)
    E2 = Props(5)/Props(4)
    eprint(i)=E
    done(i)=0
    istatdone=0
    tempstran(i)=stran(i)+dstran(i)
!
! in the interration loop? No lots of elements and NPTs
! if ((noel.eq.ielw).and.(npt.eq.5)) then
!    write(*,*) time(1),dtime,noel,npt,tempstran(i),Stran(i),dstran(i)
! endif
! start main if positive strain
!
if (dstran(i).ge.0.0) then
!
    if (tempstran(i).le.strain1) then
     elastrain(i)=tempstran(i)
     plastrain(i)=0.0
     totstrain(i)=tempstran(i)
     dstranelas(i)=dstran(i)
     dstranplas(i)=0.0
     elastress(i)=elastrain(i)*E
     plastress(i)=plastrain(i)*E2
     totstress(i)=elastress(i)+plastress(i)
     stress(i)=totstress(i)
     done(i)=1
     istatdone=1
     eprint(i)=E
     LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E / (2.0d0 * (1.0d0 + nu))
     ddsdde(1:ndi,1:ndi) = LambdaLame
     ddsdde(i,i) = ddsdde(i,i) + MuLame
    endif

   if (stran(i).gt.strain1) then
     elastrain(i)=strain1
     plastrain(i)=tempstran(i)-strain1
     totstrain(i)=tempstran(i)
     dstranelas(i)=0.0
     dstranplas(i)=dstran(i)
     elastress(i)=elastrain(i)*E
     plastress(i)=plastrain(i)*E2
     totstress(i)=elastress(i)+plastress(i)
     stress(i)=totstress(i)
     done(i)=2
     istatdone=2
     eprint(i)=E2
     LambdaLame = E2 * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E2 / (2.0d0 * (1.0d0 + nu))
     ddsdde(1:ndi,1:ndi) = LambdaLame
     ddsdde(i,i) = ddsdde(i,i) + MuLame
    endif

    if ((tempstran(i).gt.strain1).and.(stran(i).lt.strain1)) then
     elastrain(i)=strain1
     plastrain(i)=tempstran(i)-strain1
     totstrain(i)=tempstran(i)
     dstranelas(i)=strain1-stran(i)
     dstranplas(i)=tempstran(i)-strain1
     elastress(i)=elastrain(i)*E
     plastress(i)=plastrain(i)*E2
     totstress(i)=elastress(i)+plastress(i)
     stress(i)=totstress(i)
     done(i)=3
     istatdone=3
     eprint(i)=E2
     LambdaLame = E2 * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E2 / (2.0d0 * (1.0d0 + nu))
     ddsdde(1:ndi,1:ndi) = LambdaLame
     ddsdde(i,i) = ddsdde(i,i) + MuLame
    endif
!
!  main else if dstran is negative
!
else
     elastrain(i)=tempstran(i)
     plastrain(i)=0.0
     totstrain(i)=tempstran(i)
     dstranelas(i)=dstran(i)
     dstranplas(i)=0.0
     elastress(i)=elastrain(i)*E
     plastress(i)=plastrain(i)*E2
     totstress(i)=elastress(i)+plastress(i)
     stress(i)=totstress(i)
     done(i)=0
     istatdone=0
     eprint(i)=E
     LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E / (2.0d0 * (1.0d0 + nu))
     ddsdde(1:ndi,1:ndi) = LambdaLame
     ddsdde(i,i) = ddsdde(i,i) + MuLame
!
! end main if
!
end if 
 
END DO

if (istatdone.lt.3) then 
        
    DO i=1,ndi
      ddsdde(i,i) = ddsdde(i,i) + MuLame
    END DO
!    stress = stress + MATMUL(ddsdde,dstran)

endif

if (istatdone.ge.3) then
    LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
    MuLame = E / (2.0d0 * (1.0d0 + nu))
    ddsdde = 0.0d0
!    ddsdde(1:ndi,1:ndi) = LambdaLame
         
      ddsdde(i,i) = ddsdde(i,i) + MuLame
       
    DO i=1,ndi
      ddsdde(i,i) = ddsdde(i,i) + MuLame
    END DO
!    stress = stress + MATMUL(ddsdde,dstranelas) 

    LambdaLame = E2 * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
    MuLame = E2 / (2.0d0 * (1.0d0 + nu))
!   ddsdde = 0.0d0
    ddsdde(1:ndi,1:ndi) = LambdaLame
         
      ddsdde(i,i) = ddsdde(i,i) + MuLame
       
    DO i=1,ndi
      ddsdde(i,i) = ddsdde(i,i) + MuLame
    END DO
!    stress = stress + MATMUL(ddsdde,dstranplas)
!Do i=1,ntens
!stress(i)=totstress(i)
!end do

endif
! 
!    Calculate invariants 
     term1 = (elastrain(1)+elastrain(2))/2.
     term2 = ((elastrain(1)-elastrain(2))/2)**2
     term3 = (elastrain(4)/2)**2
     elaeprinmax = term1+sqrt(term2+term3)
     elaeprinmin = term1-sqrt(term2+term3)

     term1 = (plastrain(1)+plastrain(2))/2.
     term2 = ((plastrain(1)-plastrain(2))/2)**2
     term3 = (plastrain(4)/2)**2
     plaeprinmax = term1+sqrt(term2+term3)
     plaeprinmin = term1-sqrt(term2+term3)

     term=elastress(1)+elastress(2)
     term1 =  0.5*(term)
     term=elastress(1)-elastress(2) 
     term2 = (0.5*term)**2
     term3 = (stress(4))**2
     elasprinmax = term1+sqrt(term2+term3)
     elasprinmin = term1-sqrt(term2+term3)

     term=plastress(1)+plastress(2)
     term1 =  0.5*(term)
     term=plastress(1)-plastress(2) 
     term2 = (0.5*term)**2
     term3 = (plastress(4))**2
     plasprinmax = term1+sqrt(term2+term3)
     plasprinmin = term1-sqrt(term2+term3)

     elamises = sqrt(elastress(1)**2-elastress(1)*elastress(2)+elastress(2)**2+3.*elastress(4)**2)
     plamises = sqrt(plastress(1)**2-plastress(1)*plastress(2)+plastress(2)**2+3.*plastress(4)**2)
     mises = sqrt(totstress(1)**2-totstress(1)*totstress(2)+totstress(2)**2+3.*totstress(4)**2)
!    if 1     
     if (iwrite.eq.1) then
!      if 2
       if (ielw.lt.0) then
       
        write(90,fmt='(f8.3,1x,i7,1x,i7,9(1x,eS12.5))') time(1)+dtime,noel,npt,elaeprinmax,&
        plaeprinmax,elaeprinmax+plaeprinmax,elasprinmax,plasprinmax,elasprinmax+plasprinmax,&
        elamises,plamises,mises

       else
!       if 3       
        if ((noel.eq.ielw).and.(npt.eq.inpt)) then
       
 !      write(90,*) ' Time Element IntPnt done E  ElasStrain PlasStrain TotStrain ElasPrin ElasPrin &
 !      TotPrin ElasMisis PlasMises TotMises (max then min prins)'
        write(90,fmt='(f8.3,1x,i7,1x,i7,9(1x,eS12.5))') time(1)+dtime,noel,npt,elaeprinmax,&
        plaeprinmax,elaeprinmax+plaeprinmax,elasprinmax,plasprinmax,elasprinmax+plasprinmax,&
        elamises,plamises,mises
!
! write plot file if one element and one npt
!
       write(89,fmt='(es12.5,a1,es12.5,a1,es12.5)') time(1)+dtime,",",elaeprinmax+plaeprinmax,",",mises
!
! write each ntens to a unit of one element and one npt
!        
        do j=1,ntens
         iunit=90+j 
         write(iunit,fmt='(es12.5,1x,i1,1x,5(es12.5,1x))') time(1)+dtime,done(j),eprint(j),dstran(j),stran(j),&
         tempstran(j),stress(j)
!
        end do
! 
!       Else all NPTs
       
!       end 3
        endif

        if ((noel.eq.ielw).and.(npt.lt.0)) then
         write(90,fmt='(f8.3,1x,i7,1x,i7,9(1x,eS12.5))') time(1)+dtime,noel,npt,elaeprinmax,&
         plaeprinmax,elaeprinmax+plaeprinmax,elasprinmax,plasprinmax,elasprinmax+plasprinmax,&
         elamises,plamises,mises
        endif

!      end 2
       endif
!   end 1
    endif

!------------------------------------------------------------------------------
  END SUBROUTINE bi_linear
!------------------------------------------------------------------------------

