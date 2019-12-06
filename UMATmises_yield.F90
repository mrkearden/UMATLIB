!------------------------------------------------------------------------------
  SUBROUTINE mises_yield(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
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
!  Material Constants 8
!  1) Density, 2) strain at yield, 3) yield stress, 4)strain at final, 5) stress at final
!  6) Write output file 1.0=yes, 7) element number to write (-1.0) is all and all npt)
!  8) NPT to write (-1.0) for all
!------------------------------------------------------------------------------
    ! Local variables:
    INTEGER :: i,iwrite,ielw,done(ntens),iunit,j,istatdone, inpt
    REAL(KIND=dp) :: nu, E1, E2, LambdaLame, MuLame, G, G1, G2, E
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
    E1 = Props(3)/Props(2)
    E2 = (Props(5)-Props(3))/(Props(4)-Props(2))
    
    nu = Props(1)

    G1 = E1/(2*(1+nu))
    G2 = E2/(2*(1+nu))

    iwrite = int(Props(6))
    ielw = int(Props(7))
    inpt = int(Props(8))
    strain1 = Props(2)
    stressc1 = Props(3)
    strain2 = Props(4)
    stressc2 = Props(5)
!
! check mises
!
 mises = sqrt(stress(1)**2-stress(1)*stress(2)+stress(2)**2+3.*stress(4)**2)
!
!
!
if (mises.gt.props(3)) then
!
     E = E2   
     
else
     
    E = E1
     
end if 
!
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
    ! So, for this model, the other way to return the stress:
    !stress = MATMUL(ddsdde,stran+dstran)
       if (ielw.lt.0) then
       
       write(20,*) "strain xx",stran(1),"mises:",mises,"stress_xx:",stress(1)," E:",E

       else
        
        if ((noel.eq.ielw).and.(npt.eq.inpt)) then
 
        write(20,*) "strain xx",stran(1),"mises:",mises,"stress_xx:",stress(1)," E:",E
        
        endif

       endif

!------------------------------------------------------------------------------
  END SUBROUTINE mises_yield
!------------------------------------------------------------------------------
