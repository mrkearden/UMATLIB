# UMATLIB

umat_tests.zip are the test files for myUMATlib.F90

myUMATlib.F90 combines all of the material models
see materialLIBS for descriptions

UMAT LIB for Elmer. 
UMAT files for use in Elmer ElasticSolver.F90.
   
  UMATmises_yield.F90 (works for 2D and 3D elements) checks the mises stress for yielding using the first pair of stress strain points
Works with Elmer ElasticSolver.
ElasticSolver has now been modified to output state variables to vtu file.
UMATbi-linear is used with ElasticSolver and calculates invariants and puts them into state variables for output to vtu.
See example in cut test case.  Added ultimate strain for perfectly plastic behavior.

The binary plastic model is input with two sets of strain-stress points yeild is assumed at the first point
after yielding the stiffness is the slope of the second point relative to the first, Fifth entry is ultimate strain.

Example

0.002 57000. 0.162 63025. 75000.

has a yield of 57000 at a strain 0f 0.002 E=28.5E7 then E=3.89E5 after that

full input constants are poisson ratio,strain1,stress1,strain2,stress2,write flag
1=yes anything else is no, element to write -1=all elements
integration point to write -1 =all


Example in the test cases
Material 1

Name = "Material 1"

!  Poisson ratio = 0.3

!  Youngs modulus = 2.85E7

Number of Material Constants = Integer 6

Number of State Variables = Integer 7

! List material constants as {poisson strain1 stress1 strain2 stress2, ultimate strain

Material Constants(6) = Real 0.3 0.002 57000. 0.162 63025. .25

Density = 1.0

Reference Temperature = 293.0

UMAT Subroutine = File "myUMATLib" "bi_linear"

Name = "bi_linear"  ! This specifies the CMNAME argument of UMAT

End

