!------------------------------------------------------------------------------
  SUBROUTINE bi_linear(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
!------------------------------------------------------------------------------
!
! XX, YY, and XY plane b-linear plastic model, tension only
!  see bi_linear_sym for symetric compression and tension
! Bi-linear plastic material model for 2D elements 
! 8 Props input nu,strain at yield,yield stress, second strain point,stress,
! write data option, element data to write (-1 = all), npt to write
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
    REAL(KIND=dp) :: nu, E, E2, LambdaLame, MuLame, G, G2
    LOGICAL :: exists,f1exist,f2exist
    REAl(KIND=dp) :: mises, term1, term2, term3, term4, elaeprinmax, elaeprinmin
    REAL(KIND=dp) :: plaeprinmax,plaeprinmin,plasprinmax,plasprinmin, elamises, plamises
    REAL(KIND=dp) :: elastrain(ntens),plastrain(ntens),elastress(ntens),plastress(ntens),tresca
    real(kind=dp) :: edprinmax,edprinmin,elasprinmax,elasprinmin,tempstran(ntens)
    real(kind=dp) :: strain1,stressc1,strain2,stressc2,term,eprint(ntens)
    REAL(KIND=dp) :: totstrain(ntens),totstress(ntens)
    REAL(KIND=dp) :: dstranelas(NTENS),dstranplas(ntens),astatev(16)
!------------------------------------------------------------------------------

    ! Get Young's modulus and the Poisson ratio:
    E = Props(3)/Props(2)
    E2 = Props(5)/Props(4)
    
    nu = Props(1)

    G = E/(2*(1+nu))
    G2 = E2/(2*(1+nu))

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
    G = E/(2*(1+nu))
    G2 = E2/(2*(1+nu))
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
if (dstran(i).gt.0.0) then
!
    if (tempstran(i).le.strain1) then
     elastrain(i)=tempstran(i)
     plastrain(i)=0.0
     totstrain(i)=tempstran(i)
     dstranelas(i)=dstran(i)
     dstranplas(i)=0.0
     done(i)=1
     istatdone=1
     eprint(i)=E
     LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E / (2.0d0 * (1.0d0 + nu))
     ddsdde(1:ndi,1:ndi) = LambdaLame
     ddsdde(i,i) = ddsdde(i,i) + MuLame
    DO j=1,ntens
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
    endif

   if (stran(i).gt.strain1) then
     elastrain(i)=strain1
     plastrain(i)=tempstran(i)-strain1
     totstrain(i)=tempstran(i)
     dstranelas(i)=0.0
     dstranplas(i)=dstran(i)
     done(i)=2
     istatdone=2
     eprint(i)=E2
     LambdaLame = E2 * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E2 / (2.0d0 * (1.0d0 + nu))
     ddsdde(1:ndi,1:ndi) = LambdaLame
     ddsdde(i,i) = ddsdde(i,i) + MuLame
    DO j=1,ntens
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
    endif

    if ((tempstran(i).gt.strain1).and.(stran(i).lt.strain1)) then
     elastrain(i)=strain1
     plastrain(i)=tempstran(i)-strain1
     totstrain(i)=tempstran(i)
     dstranelas(i)=strain1-stran(i)
     dstranplas(i)=tempstran(i)-strain1
     done(i)=3
     istatdone=3
     eprint(i)=E2
     LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E / (2.0d0 * (1.0d0 + nu))
     ddsdde = 0.0d0
     ddsdde(1:ndi,1:ndi) = LambdaLame
         
     ddsdde(i,i) = ddsdde(i,i) + MuLame
       
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO


    LambdaLame = E2 * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
    MuLame = E2 / (2.0d0 * (1.0d0 + nu))
!   ddsdde = 0.0d0
    ddsdde(1:ndi,1:ndi) = LambdaLame
         
      ddsdde(i,i) = ddsdde(i,i) + MuLame
       
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO

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
     done(i)=4
     istatdone=4
     eprint(i)=E
     LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E / (2.0d0 * (1.0d0 + nu))

    ddsdde = 0.0d0
    ddsdde(1:ndi,1:ndi) = LambdaLame
    DO j=1,ntens
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
!
! end main if
!
end if 
 
END DO


        
     elastress(4)=2*elastrain(4)*G
     elastress(1)=(E*elastrain(1)-nu*elastress(4))/(1+nu)
     elastress(2)=(E*elastrain(2)+nu*elastress(1))
     elastress(3)=0.0
     plastress(4)=2*plastrain(4)*G2
     plastress(1)=(E2*plastrain(1)-nu*plastress(4))/(1+nu)
     plastress(2)=(E2*plastrain(2)+nu*plastress(1))
     plastress(3)=0.0

! 
!    Calculate invariants
     do j=1,ntens
     totstress(j)=elastress(j)+plastress(j)
     stress(j)=totstress(j) 
     end do 
     term1 = (elastrain(1)+elastrain(2))/2.
     term2 = ((elastrain(1)-elastrain(2))/2.)**2
     term3 = (elastrain(4)/2.)**2
     elaeprinmax = term1+sqrt(term2+term3)
     elaeprinmin = term1-sqrt(term2+term3)

     term1 = (plastrain(1)+plastrain(2))/2.
     term2 = ((plastrain(1)-plastrain(2))/2.)**2
     term3 = (plastrain(4)/2.)**2
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

     astatev(1)=elaeprinmax
     astatev(2)=plaeprinmax
     astatev(3)=elaeprinmax+plaeprinmax
     astatev(4)=elasprinmax
     astatev(5)=plasprinmax
     astatev(6)=elasprinmax+plasprinmax
     astatev(7)=elamises
     astatev(8)=plamises
     astatev(9)=mises
     astatev(10)=elaeprinmin
     astatev(11)=plaeprinmin
     astatev(12)=elaeprinmin+plaeprinmin
     astatev(13)=elasprinmin
     astatev(14)=plasprinmin
     astatev(15)=elasprinmin+plasprinmin
     astatev(16)=time(1)

!    if 1     
     if (iwrite.eq.1) then
     inquire(file="maxprindata.txt",EXIST=f1exist)
     inquire(file="minprindata.txt",EXIST=f2exist)
     if(.not.f2exist) then
      open(unit=90,file="minprindata.txt",status='new',recl=145)
      write(90,*) '  Time   Element  IntPnt  ElasStrain  PlasStrain   TotStrain    ElasPrin     PlasPrin    &
  TotPrin       ElasMisis    PlasMises   TotMises'  
     endif
     if(.not.f1exist) then
      open(unit=91,file="maxprindata.txt",status='new',recl=145)
      write(91,*) '  Time   Element  IntPnt  ElasStrain  PlasStrain   TotStrain    ElasPrin     PlasPrin    &
  TotPrin       ElasMisis    PlasMises   TotMises'  
     endif


!      if 2
       if (ielw.lt.0) then
       
        write(90,fmt='(f8.3,1x,i7,1x,i7,9(1x,eS12.5))') astatev(16),noel,npt,astatev(10),&
        astatev(11),astatev(12),astatev(13),astatev(14),astatev(15),astatev(7),astatev(8),astatev(9)
        write(91,fmt='(f8.3,1x,i7,1x,i7,6(1x,eS12.5))') astatev(16),noel,npt,astatev(1),&
        astatev(2),astatev(3),astatev(4),astatev(5),astatev(6),astatev(7),astatev(8),astatev(9)

       else
!       if 3       
        if ((noel.eq.ielw).and.(npt.eq.inpt)) then
 
 !      write(90,*) ' Time Element IntPnt done E  ElasStrain PlasStrain TotStrain ElasPrin ElasPrin &
 !      TotPrin ElasMisis PlasMises TotMises (max then min prins)'
        write(90,fmt='(f8.3,1x,i7,1x,i7,9(1x,eS12.5))') astatev(16),noel,npt,astatev(10),&
        astatev(11),astatev(12),astatev(13),astatev(14),astatev(15),astatev(7),astatev(8),astatev(9)
        write(91,fmt='(f8.3,1x,i7,1x,i7,9(1x,eS12.5))') astatev(16),noel,npt,astatev(1),&
        astatev(2),astatev(3),astatev(4),astatev(5),astatev(6),astatev(7),astatev(8),astatev(9)
       
         
! write each ntens to a unit of one element and one npt
!        
!        do j=1,ntens
!         iunit=92+j 
!         write(iunit,fmt='(es12.5,1x,i1,1x,5(es12.5,1x))') time(1)+dtime,done(j),eprint(j),dstran(j),stran(j),&
!         tempstran(j),stress(j)
!
!        end do
       
!       end 3
        endif

        if ((noel.eq.ielw).and.(npt.lt.0)) then
         write(90,fmt='(f8.3,1x,i7,1x,i7,9(1x,eS12.5))') astatev(16),noel,npt,astatev(10),&
         astatev(11),astatev(12),astatev(13),astatev(14),astatev(15),astatev(7),astatev(8),astatev(9)
         write(91,fmt='(f8.3,1x,i7,1x,i7,9(1x,eS12.5))') astatev(16),noel,npt,astatev(1),&
         astatev(2),astatev(3),astatev(4),astatev(5),astatev(6),astatev(7),astatev(8),astatev(9)
        endif

!      end 2
       endif
!   end 1
    endif
!------------------------------------------------------------------------------
  END SUBROUTINE bi_linear
!------------------------------------------------------------------------------
! The call from ElasticSolve
!CALL UMATusersubrtn(UMATSubrtn, StressVec(1:ntens), StateV, StressDer(1:ntens,1:ntens), EnergyElast, &
!          EnergyPlast, EnergyVisc, rpl, ddsddt(1:ntens), drplde(1:ntens), drpldt, &
!          stran(1:ntens), dstran(1:ntens), TimeAtStep, dtime, Temp, dTemp, &
!          predef, dpred, cmname, ndi, nshr, ntens, NStateV, InProps, NrInProps, coords, &
!          drot, pnewdt, celent, DefG0, DefG, ElementIndex, t, layer, kspt, kstep, kinc)
!------------------------------------------------------------------------------
  SUBROUTINE bi_linear_sym(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
!------------------------------------------------------------------------------
!
! XX, YY, and XY plane b-linear plastic model symmetric compression and tension
!
! Bi-linear plastic material model for 2D elements 
! 8 Props input nu,strain at yield,yield stress, second strain point,stress,
! write data option, element data to write (-1 = all), npt to write
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
    REAL(KIND=dp) :: nu, E, E2, LambdaLame, MuLame, G, G2
    LOGICAL :: f1exist,f2exist
    REAl(KIND=dp) :: mises, term1, term2, term3, term4, elaeprinmax, elaeprinmin
    REAL(KIND=dp) :: plaeprinmax,plaeprinmin,plasprinmax,plasprinmin, elamises, plamises
    REAL(KIND=dp) :: elastrain(ntens),plastrain(ntens),elastress(ntens),plastress(ntens),tresca
    real(kind=dp) :: edprinmax,edprinmin,elasprinmax,elasprinmin,tempstran(ntens)
    real(kind=dp) :: strain1,stressc1,strain2,stressc2,term,eprint(ntens)
    REAL(KIND=dp) :: totstrain(ntens),totstress(ntens)
    REAL(KIND=dp) :: dstranelas(NTENS),dstranplas(ntens),astatev(16)
!------------------------------------------------------------------------------

    ! Get Young's modulus and the Poisson ratio:
    E = Props(3)/Props(2)
    E2 = Props(5)/Props(4)
    
    nu = Props(1)

    G = E/(2*(1+nu))
    G2 = E2/(2*(1+nu))

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
! If dstran is negative assume compression symetric to tension
!
    E = Props(3)/Props(2)
    E2 = Props(5)/Props(4)
    G = E/(2*(1+nu))
    G2 = E2/(2*(1+nu))
    eprint(i)=E
    done(i)=0
    istatdone=0
    tempstran(i)=stran(i)+dstran(i)

if (dstran(i).gt.0.0) then
!! DSTRAN is loading if dstran is positive
    if (tempstran(i).le.strain1) then
     elastrain(i)=tempstran(i)
     plastrain(i)=0.0
     totstrain(i)=tempstran(i)
     dstranelas(i)=dstran(i)
     dstranplas(i)=0.0
     done(i)=1
     istatdone=1
     eprint(i)=E
     LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E / (2.0d0 * (1.0d0 + nu))
     ddsdde(1:ndi,1:ndi) = LambdaLame
     ddsdde(i,i) = ddsdde(i,i) + MuLame
    DO j=1,ntens
      ddsdde(j,j) = ddsdde(j,j) + MuLame

    END DO
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
    endif

   if (stran(i).gt.strain1) then
     elastrain(i)=strain1
     plastrain(i)=tempstran(i)-strain1
     totstrain(i)=tempstran(i)
     dstranelas(i)=0.0
     dstranplas(i)=dstran(i)
     done(i)=2
     istatdone=2
     eprint(i)=E2
     LambdaLame = E2 * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E2 / (2.0d0 * (1.0d0 + nu))
     ddsdde(1:ndi,1:ndi) = LambdaLame
     ddsdde(i,i) = ddsdde(i,i) + MuLame
    DO j=1,ntens
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
    endif

    if ((tempstran(i).gt.strain1).and.(stran(i).lt.strain1)) then
     elastrain(i)=strain1
     plastrain(i)=tempstran(i)-strain1
     totstrain(i)=tempstran(i)
     dstranelas(i)=strain1-stran(i)
     dstranplas(i)=tempstran(i)-strain1
     done(i)=3
     istatdone=3
     eprint(i)=E2
     LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E / (2.0d0 * (1.0d0 + nu))
     ddsdde = 0.0d0
     ddsdde(1:ndi,1:ndi) = LambdaLame
         
     ddsdde(i,i) = ddsdde(i,i) + MuLame
       
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO


    LambdaLame = E2 * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
    MuLame = E2 / (2.0d0 * (1.0d0 + nu))
!   ddsdde = 0.0d0
    ddsdde(1:ndi,1:ndi) = LambdaLame
         
      ddsdde(i,i) = ddsdde(i,i) + MuLame
       
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO

    endif
!
!  main else if dstran is negative
!
else
!! dstran is negative could be unloading or compression
if (tempstran(i).lt.0.0) then
! compression
    if (abs(tempstran(i)).le.strain1) then
     elastrain(i)=tempstran(i)
     plastrain(i)=0.0
     totstrain(i)=tempstran(i)
     dstranelas(i)=dstran(i)
     dstranplas(i)=0.0
     done(i)=5
     istatdone=5
     eprint(i)=E
     LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E / (2.0d0 * (1.0d0 + nu))
     ddsdde(1:ndi,1:ndi) = LambdaLame
     ddsdde(i,i) = ddsdde(i,i) + MuLame
    DO j=1,ntens
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
    endif

   if (abs(stran(i)).gt.strain1) then
     elastrain(i)=-1.*strain1
     plastrain(i)=tempstran(i)+strain1
     totstrain(i)=tempstran(i)
     dstranelas(i)=0.0
     dstranplas(i)=dstran(i)
     done(i)=6
     istatdone=6
     eprint(i)=E2
     LambdaLame = E2 * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E2 / (2.0d0 * (1.0d0 + nu))
     ddsdde(1:ndi,1:ndi) = LambdaLame
     ddsdde(i,i) = ddsdde(i,i) + MuLame
    DO j=1,ntens
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
    endif

    if ( (abs(tempstran(i)) .gt. strain1) .and. (abs(stran(i)) .lt. strain1) ) then
     elastrain(i)=-1*strain1
     plastrain(i)=tempstran(i)+strain1
     totstrain(i)=tempstran(i)
     dstranelas(i)=-1*strain1+stran(i)
     dstranplas(i)=tempstran(i)+strain1
     done(i)=7
     istatdone=7
     eprint(i)=E2
     LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E / (2.0d0 * (1.0d0 + nu))
     ddsdde = 0.0d0
     ddsdde(1:ndi,1:ndi) = LambdaLame
         
     ddsdde(i,i) = ddsdde(i,i) + MuLame
       
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO


    LambdaLame = E2 * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
    MuLame = E2 / (2.0d0 * (1.0d0 + nu))
!   ddsdde = 0.0d0
    ddsdde(1:ndi,1:ndi) = LambdaLame
         
      ddsdde(i,i) = ddsdde(i,i) + MuLame
       
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO

    endif
else
! unloading
     elastrain(i)=tempstran(i)
     plastrain(i)=0.0
     totstrain(i)=tempstran(i)
     dstranelas(i)=dstran(i)
     dstranplas(i)=0.0
     done(i)=4
     istatdone=4
     eprint(i)=E
     LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
     MuLame = E / (2.0d0 * (1.0d0 + nu))

    ddsdde = 0.0d0
    ddsdde(1:ndi,1:ndi) = LambdaLame
    DO j=1,ntens
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
    DO j=1,ndi
      ddsdde(j,j) = ddsdde(j,j) + MuLame
    END DO
endif
!
! end main if
!
end if 
 
END DO


        
     elastress(4)=2*elastrain(4)*G
     elastress(1)=(E*elastrain(1)-nu*elastress(4))/(1+nu)
     elastress(2)=(E*elastrain(2)+nu*elastress(1))
     elastress(3)=0.0
     plastress(4)=2*plastrain(4)*G2
     plastress(1)=(E2*plastrain(1)-nu*plastress(4))/(1+nu)
     plastress(2)=(E2*plastrain(2)+nu*plastress(1))
     plastress(3)=0.0

! 
!    Calculate invariants
     do j=1,ntens
     totstress(j)=elastress(j)+plastress(j)
     stress(j)=totstress(j) 
     end do 
     term1 = (elastrain(1)+elastrain(2))/2.
     term2 = ((elastrain(1)-elastrain(2))/2.)**2
     term3 = (elastrain(4)/2.)**2
     elaeprinmax = term1+sqrt(term2+term3)
     elaeprinmin = term1-sqrt(term2+term3)

     term1 = (plastrain(1)+plastrain(2))/2.
     term2 = ((plastrain(1)-plastrain(2))/2.)**2
     term3 = (plastrain(4)/2.)**2
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

     astatev(1)=elaeprinmax
     astatev(2)=plaeprinmax
     astatev(3)=elaeprinmax+plaeprinmax
     astatev(4)=elasprinmax
     astatev(5)=plasprinmax
     astatev(6)=elasprinmax+plasprinmax
     astatev(7)=elamises
     astatev(8)=plamises
     astatev(9)=mises
     astatev(10)=elaeprinmin
     astatev(11)=plaeprinmin
     astatev(12)=elaeprinmin+plaeprinmin
     astatev(13)=elasprinmin
     astatev(14)=plasprinmin
     astatev(15)=elasprinmin+plasprinmin
     astatev(16)=time(1)


!    if 1     
     if (iwrite.eq.1) then
     inquire(file="maxprindata.txt",EXIST=f1exist)
     inquire(file="minprindata.txt",EXIST=f2exist)
     if(.not.f2exist) then
      open(unit=90,file="minprindata.txt",status='new',recl=145)
      write(90,*) '  Time   Element  IntPnt  ElasStrain  PlasStrain   TotStrain    ElasPrin     PlasPrin    &
  TotPrin       ElasMisis    PlasMises   TotMises'  
     endif
     if(.not.f1exist) then
      open(unit=91,file="maxprindata.txt",status='new',recl=145)
      write(91,*) '  Time   Element  IntPnt  ElasStrain  PlasStrain   TotStrain    ElasPrin     PlasPrin    &
  TotPrin       ElasMisis    PlasMises   TotMises'  
     endif


!      if 2
       if (ielw.lt.0) then
       
        write(90,fmt='(f8.3,1x,i7,1x,i7,9(1x,eS12.5))') astatev(16),noel,npt,astatev(10),&
        astatev(11),astatev(12),astatev(13),astatev(14),astatev(15),astatev(7),astatev(8),astatev(9)
        write(91,fmt='(f8.3,1x,i7,1x,i7,6(1x,eS12.5))') astatev(16),noel,npt,astatev(1),&
        astatev(2),astatev(3),astatev(4),astatev(5),astatev(6),astatev(7),astatev(8),astatev(9)

       else
!       if 3       
        if ((noel.eq.ielw).and.(npt.eq.inpt)) then
 
 !      write(90,*) ' Time Element IntPnt done E  ElasStrain PlasStrain TotStrain ElasPrin ElasPrin &
 !      TotPrin ElasMisis PlasMises TotMises (max then min prins)'
        write(90,fmt='(f8.3,1x,i7,1x,i7,9(1x,eS12.5))') astatev(16),noel,npt,astatev(10),&
        astatev(11),astatev(12),astatev(13),astatev(14),astatev(15),astatev(7),astatev(8),astatev(9)
        write(91,fmt='(f8.3,1x,i7,1x,i7,9(1x,eS12.5))') astatev(16),noel,npt,astatev(1),&
        astatev(2),astatev(3),astatev(4),astatev(5),astatev(6),astatev(7),astatev(8),astatev(9)
       
         
! write each ntens to a unit of one element and one npt
!        
!        do j=1,ntens
!         iunit=92+j 
!         write(iunit,fmt='(es12.5,1x,i1,1x,5(es12.5,1x))') time(1)+dtime,done(j),eprint(j),dstran(j),stran(j),&
!         tempstran(j),stress(j)
!
!        end do
       
!       end 3
        endif

        if ((noel.eq.ielw).and.(npt.lt.0)) then
         write(90,fmt='(f8.3,1x,i7,1x,i7,9(1x,eS12.5))') astatev(16),noel,npt,astatev(10),&
         astatev(11),astatev(12),astatev(13),astatev(14),astatev(15),astatev(7),astatev(8),astatev(9)
         write(91,fmt='(f8.3,1x,i7,1x,i7,9(1x,eS12.5))') astatev(16),noel,npt,astatev(1),&
         astatev(2),astatev(3),astatev(4),astatev(5),astatev(6),astatev(7),astatev(8),astatev(9)
        endif

!      end 2
       endif
!   end 1
    endif

!------------------------------------------------------------------------------
  END SUBROUTINE bi_linear_sym
!------------------------------------------------------------------------------

