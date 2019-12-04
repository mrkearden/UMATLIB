# UMATLIB
UMAT LIB for Elmer. 
UMAT files for use in Elmer ElasticSolver.F90.

  UMATLib.F90 is the UMATLib provieded with Elmer and contains a linear istropic model
  and a hencky_stvenant_kirchhoff model.
  
  UMATLib_test.F90 contains a 2D bi-linear plastic material model. Only working for XX and YY tensile, compression is assumed to be non-yielding. 2D elements.
 
UMAT_Plastic.F90 2D bi-linear plastic material model for XX, YY, XY stress. bi_linear UMAT is Tensile yielding only. bi_linear_sym UMAT assumes the bi-linear model is symetric in compression.

The binary plastic model is input with two sets of strain-stress points yeild is assumed at the first point
after yielding the stiffness is the slope of the second point relative to the first
Example
0.002 57000. 0.162 63025.
has a yield of 57000 at a strain 0f 0.002 E=28.5E7 then E=3.89E5 after that

full input constants are poisson ratio,strain1,stress1,strain2,stress2,write flag
1=yes anything else is no, element to write -1=all elements
integration point to write -1 =all


Example in the test cases
Material 1
Name = "Material 1"
!  Poisson ratio = 0.3
!  Youngs modulus = 2.85E7
Number of Material Constants = Integer 8
Number of State Variables = Integer 0
! List material constants as {poisson strain1 stress1 strain2 stress2,
!  write output 1 is yes element to write -1=all npt -1=all}:
Material Constants(8) = Real 0.3 0.002 57000. 0.162 63025. 1.0 901.0 5.0
Density = 1.0
Reference Temperature = 293.0
UMAT Subroutine = File "UMATplastic" "bi_linear"
Name = "bi_linear"  ! This specifies the CMNAME argument of UMAT
End
