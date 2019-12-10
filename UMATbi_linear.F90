!------------------------------------------------------------------------------
  SUBROUTINE bi_linear(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
! the call from Elmer
!CALL UMATusersubrtn(UMATSubrtn, StressVec(1:ntens), StateV, StressDer(1:ntens,1:ntens), EnergyElast, &
!          EnergyPlast, EnergyVisc, rpl, ddsddt(1:ntens), drplde(1:ntens), drpldt, &
!          stran(1:ntens), dstran(1:ntens), TimeAtStep, dtime, Temp, dTemp, &
!          predef, dpred, cmname, ndi, nshr, ntens, NStateV, InProps, NrInProps, coords, &
!          drot, pnewdt, celent, DefG0, DefG, ElementIndex, t, layer, kspt, kstep, kinc)
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
!  Material Constants 5
!  1) Nu, 2) strain at yield, 3) yield stress, 4)strain at final, 5) stress at final
!  State Variables 7
! 3 prinipal strains, 3 prinipal stress, mises
!------------------------------------------------------------------------------
    ! Local variables:
    REAL(KIND=dp) :: nu, E1, E2, LambdaLame, MuLame, G, G1, G2, E
    LOGICAL :: exists,f1exist,f2exist
    REAl(KIND=dp) :: mises, term1, term2, term3, term4, eprinmax, eprinmin
    REAL(KIND=dp) :: eprinmid,sprinmax,sprinmid,sprinmin
    real(kind=dp) :: strain1,stressc1,strain2,stressc2,term
    REAL(KIND=dp) :: totstran(ntens)
    Real*8 AA(3,3), ZZ(3,3)
    integer :: i
!------------------------------------------------------------------------------
  ! Get Young's modulus and the Poisson ratio:
    E1 = Props(3)/Props(2)
    E2 = (Props(5)-Props(3))/(Props(4)-Props(2))
    
    nu = Props(1)

    G1 = E1/(2*(1+nu))
    G2 = E2/(2*(1+nu))

    strain1 = Props(2)
    stressc1 = Props(3)
    strain2 = Props(4)
    stressc2 = Props(5)
    
!
! check mises
E=E1
!
if (ntens.eq.4) then
 mises = sqrt(stress(1)**2-stress(1)*stress(2)+stress(2)**2+3.*stress(4)**2)
     term1=stress(1)-stress(2) 
     term2 = (0.5*term1)**2
     term3 = (stress(4))**2
     sprinmax = term1+sqrt(term2+term3)
     sprinmin = term1-sqrt(term2+term3)
end if
if (ntens.eq.6) then
 term1=(stress(1)-stress(2))**2
 term2=(stress(2)-stress(3))**2
 term3=(stress(3)-stress(1))**2
 term4=6.0*(stress(4)**2+stress(5)**2+stress(6)**2)
 mises=sqrt(term1+term2+term3+term4)
endif
!
!
if (mises.gt.props(3)) then
!
     E = E2   
     
else
     
    E = E1
     
end if 
if (mises.gt.props(6)) then
  E = 0.0
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
   stress = stress + MATMUL(ddsdde,dstran)
   totstran = stran + dstran
!
! Calculate invariants for statev
!
Select case(ntens)
     case(4) 
     mises = sqrt(stress(1)**2-stress(1)*stress(2)+stress(2)**2+3.*stress(4)**2)
     term1=stress(1)-stress(2) 
     term2 = (0.5*term1)**2
     term3 = (stress(4))**2
     sprinmax = term1+sqrt(term2+term3)
     sprinmin = term1-sqrt(term2+term3)
     term1=totstran(1)-totstran(2) 
     term2 = (0.5*term1)**2
     term3 = (totstran(4))**2
     eprinmax = term1+sqrt(term2+term3)
     eprinmin = term1-sqrt(term2+term3)
     statev(1) = eprinmax
     statev(2) = 0.0
     statev(3) = eprinmin
     statev(4) = sprinmax
     statev(5) = 0.0
     statev(6) = sprinmin
     statev(7) = mises
     case(6)
     term1=(stress(1)-stress(2))**2
    term2=(stress(2)-stress(3))**2
    term3=(stress(3)-stress(1))**2
    term4=6.0*(stress(4)**2+stress(5)**2+stress(6)**2)
    mises=sqrt(term1+term2+term3+term4)
    AA(1,1)=stress(1)
    AA(2,2)=stress(2)
    AA(3,3)=stress(3)
    AA(1,2)=stress(4)
    AA(1,3)=stress(6)
    AA(2,3)=stress(5) 
    AA(2,1)=AA(1,2)
    AA(3,1)=AA(1,3)
    AA(3,2)=AA(2,3)
     call jacobi2(AA,3,ZZ,1.E-6,100)
     sprinmax=AA(1,1)
     sprinmid=AA(2,2)
     sprinmin=AA(3,3)
    AA(1,1)=totstran(1)
    AA(2,2)=totstran(2)
    AA(3,3)=totstran(3)
    AA(1,2)=totstran(4)
    AA(1,3)=totstran(6)
    AA(2,3)=totstran(5) 
    AA(2,1)=AA(1,2)
    AA(3,1)=AA(1,3)
    AA(3,2)=AA(2,3)
     call jacobi2(AA,3,ZZ,1.E-6,100)
     eprinmax=AA(1,1)
     eprinmid=AA(2,2)
     eprinmin=AA(3,3)
     statev(1) = eprinmax
     statev(2) = eprinmid
     statev(3) = eprinmin
     statev(4) = sprinmax
     statev(5) = sprinmid
     statev(6) = sprinmin
     statev(7) = mises
 end select
!------------------------------------------------------------------------------
  END SUBROUTINE bi_linear
!------------------------------------------------------------------------------
!     
!
      SUBROUTINE JACOBI2 (AO,N,A,EPS,ITMAX)
      Real*8 :: AO(3,3), A(3,3),Z,Y,dif,tang,cose,sine,zz,yy
      Integer iter,i,j,itmax,N,nm1,ip1,irow,icol
      Real eps
      ITER=0
      DO I=1,N
       DO J=1,N
        A(I,J)=0.0
        A(I,I)=1.0
       end do
      end do
    
    do while (iter.lt.itmax)
      Z=0.
      NM1=N-1 

      DO I=1,NM1
       IP1=I+1
       DO J=IP1,N
        IF (ABS(AO(I,J)) .gt. Z) then
         Z=ABS(AO(I,J))
         IROW=I
         ICOL=J
        end if
       end do
      end do
   
   IF (ITER .EQ. 0) Y=Z*EPS
 
     IF (Z .gt. Y) then
       DIF=AO(IROW,IROW)-AO(ICOL,ICOL)
       TANG= (-DIF+SQRT(DIF**2+4.0*Z**2))/(2.0*AO(IROW,ICOL))
       COSE=1.0/SQRT(1.0+TANG**2)
       SINE=COSE*TANG
       DO I=1,N
        ZZ=A(I,IROW)
        A(I,IROW)=COSE*ZZ+SINE*A(I,ICOL)
        A(I,ICOL)=COSE*A(I,ICOL)-SINE*ZZ
       end do
       I=1
       do while (I .ne. IROW)  
        YY=AO(I,IROW)
        AO(I,IROW)=COSE*YY+SINE*AO(I,ICOL)
        AO(I,ICOL)=COSE*AO(I,ICOL)-SINE*YY
        I=I+1
       end do
      I=IROW+1
      do while (I .ne. ICOL) 
       YY=AO(IROW,I)
       AO(IROW,I)=COSE*YY+SINE*AO(I,ICOL)
       AO(I,ICOL)=COSE*AO(I,ICOL)-SINE*YY
       I=I+1
      end do
      I=ICOL+1
      do while (I .le. N)  
      ZZ=AO(IROW,I)
      AO(IROW,I)=COSE*ZZ+SINE*AO(ICOL,I)
      AO(ICOL,I)=COSE*AO(ICOL,I)-SINE*ZZ
      I=I+1
      end do
      YY=AO(IROW,IROW)
      AO(IROW,IROW)=YY*COSE**2+AO(IROW,ICOL)*2.0*COSE*SINE+&
AO(ICOL,ICOL)*SINE**2
      AO(ICOL,ICOL)=AO(ICOL,ICOL)*COSE**2+YY*SINE**2-AO(IROW,ICOL)&
*2.0*COSE*SINE
      AO(IROW,ICOL)=0.0
      
   else
   iter=itmax 
   end if
 
 iter=iter+1 
 end do 
 END subroutine jacobi2
