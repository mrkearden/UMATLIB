!------------------------------------------------------------------------------
  SUBROUTINE bi_linear(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
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
!  6) Write output file 1.0=yes, 7) element number to write (-1.0) is all)
!------------------------------------------------------------------------------
    ! Local variables:
    INTEGER :: i,iwrite,ielw
    REAL(KIND=dp) :: nu, E, E2, LambdaLame, MuLame
    LOGICAL :: exists
    REAl(KIND=dp) :: mises, term1, term2, term3, term4, eprinmax, eprinmin
    REAL(KIND=dp) :: elastrain,plastrain,elastress,plastress,tresca
    real(kind=dp) :: edprinmax,edprinmin,sprinmax,sprinmin
    real(kind=dp) :: strain1,stressc1,strain2,stressc2,term
    REAL(KIND=dp) :: Stresselas1(ntens),stressplas(ntens)
    REAL(KIND=dp) :: dstranelas(NTENS),dstranpla(ntens)
!------------------------------------------------------------------------------

    ! Get Young's modulus and the Poisson ratio:
    E = Props(3)/Props(2)
    E2 = Props(5)/Props(4)
!   write(*,*) E2
    nu = Props(1)
    iwrite = int(Props(6))
    ielw = int(Props(7))
    strain1 = Props(2)
    stressc1 = Props(3)
    strain2 = Props(4)
    stressc2 = Props(5)
    if (iwrite.eq.1) then
     if (noel.eq.ielw) then
      Write(*,fmt='(a5,Es13.6)') " E = ",E
      write(*,fmt='(a5,Es13.6)') " E2 =",E2
      Write(*,*) "stress xx     stress yy    stress  xy"
      write(*,fmt='(3(1x,Es13.6))') stress(1),stress(2),stress(4)
     endif
    endif
 
!
! check to see if mises over yield
!
    mises = sqrt(stress(1)**2-stress(1)*stress(2)+stress(2)**2+3.*stress(4)**2)
    if (mises.gt.stressc1) then
!
!   calculate elastic strain, plastic, strain, total strain, 
!   increment of elastic strain, increment of plastic strain,
!   elastic stress, plastic stress, total stress, and set E to new value
    do i=1,ntens
    if (stran(i).gt.strain1) then
     dstranelas(i)=0.0
     dstranplas(i)=dstran(i)
     E=E2
    endif
   if (iwrite.eq.1) then
    if (noel.eq.ielw) then
      Write(*,fmt='(a9,Es13.6)') " Mises = ",mises
      write(*,fmt='(a10,Es13.6)') " Stressc1:",stressc1
      write(*,fmt='(a5,Es13.6)') " E =",E
    endif
   endif
!
    LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
    MuLame = E / (2.0d0 * (1.0d0 + nu))
    ddsdde = 0.0d0
    ddsdde(1:ndi,1:ndi) = LambdaLame
    DO i=1,ntens
     
      ddsdde(i,i) = ddsdde(i,i) + MuLame
       
    END DO
    DO i=1,ndi
      ddsdde(i,i) = ddsdde(i,i) + MuLame
    END DO
!      write(93,fmt='(10a13)') "     TIME    ","   ELEMENT   ","  INT POINT  ",&
!     "  STRAN "," DSTRAN  "," DDSDDE "
    !
    ! We have a linear response function, so the following update is precise
    ! (no higher-order terms related to the notion of differentiability).
    ! Note that we could also define
    !
    !        stress = stress_response_function(stran + dstran)
    ! or
    !        stress = stress_response_function(dfrgrd1)
    !
    ! which may be the precise definition of the functionality required. 

     stress = stress + MATMUL(ddsdde,dstran)
!     stress(3)=0.0
!     stress1 = stress + MATMUL(ddsdde,dstran1)
!    
!    Calculate invariants 
     term1 = (stran(1)+stran(2))/2.
     term2 = ((stran(1)-stran(2))/2)**2
     term3 = (stran(4)/2)**2
     eprinmax = term1+sqrt(term2+term3)
     eprinmin = term1-sqrt(term2+term3)

     term1 = (dstran(1)+dstran(2))/2.
     term2 = ((dstran(1)-dstran(2))/2)**2
     term3 = (dstran(4)/2.)**2
     edprinmax = term1+sqrt(term2+term3)
     edprinmin = term1-sqrt(term2+term3)

     term=stress(1)+stress(2)
     term1 =  0.5*(term)
     term=stress(1)-stress(2) 
     term2 = (0.5*term)**2
     term3 = (stress(4))**2
     sprinmax = term1+sqrt(term2+term3)
     sprinmin = term1-sqrt(term2+term3)

     mises = sqrt(stress(1)**2-stress(1)*stress(2)+stress(2)**2+3.*stress(4)**2)
     tresca=0.5*(sprinmax-sprinmin)
     
     if (iwrite.eq.1) then

       if (ielw.lt.0) then
 
        write(92,fmt='(f8.3,3x,i7,6x,i7,3x,7(1x,eS12.5))') time(1),noel,npt,eprinmax,&
        eprinmin,edprinmax,tresca,sprinmax,sprinmin,mises
       else
       
        if (noel.eq.ielw) then
 
         write(92,fmt='(f8.3,3x,i7,6x,i7,3x,7(1x,eS12.5))') time(1)+dtime,noel,npt,eprinmax+edprinmax,&
         eprinmin+edprinmin,edprinmax,tresca,sprinmax,sprinmin,mises

         write(93,fmt='(f8.3,1x,i8,1x,i8,4(1x,ES13.5))') time(1)+dtime,noel,npt,stress(1),stress(2),stress(3),stress(4)

        endif
       endif
    endif
 
    ! So, for this model, the other way to return the stress:
    !stress = MATMUL(ddsdde,stran+dstran)

!------------------------------------------------------------------------------
  END SUBROUTINE bi_linear
!------------------------------------------------------------------------------

