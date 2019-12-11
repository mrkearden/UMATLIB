# UMATLIB
UMAT LIB for Elmer. 
UMAT files for use in Elmer ElasticSolver.F90.
   
  UMATmises_yield.F90 (works for 2D and 3D elements) checks the mises stress for yielding using the first pair of stress strain points
Works with Elmer ElasticSolver.
ElasticSolver has now been modified to output state variables to vtu file.
UMATbi-linear is used with ElasticSolver and calculates invariants and puts them into state variables for output to vtu.
See example in cut test case.  Added ultimate stress for perfectly plastic behavior.

PlasticSolver in this repository can be used in place of elasticsolver, Plastic Solver is the same but outputs a file called output.txt
That has prinipal stress and mises for one element and one NPT chosen by the sif input.

The binary plastic model is input with two sets of strain-stress points yeild is assumed at the first point
after yielding the stiffness is the slope of the second point relative to the first, Fifth entry is ultimate stress.

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
Number of Material Constants = Integer 8
Number of State Variables = Integer 0
! List material constants as {poisson strain1 stress1 strain2 stress2,
!  write output 1 is yes element to write -1=all npt -1=all}:
Material Constants(8) = Real 0.3 0.002 57000. 0.162 63025. 1.0 901.0 5.0
Density = 1.0
Reference Temperature = 293.0
UMAT Subroutine = File "UMATmises" "mises_yield"
Name = "mises_yield"  ! This specifies the CMNAME argument of UMAT
End
