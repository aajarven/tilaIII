function acosh ( x )
!
!! ACOSH computes the arc hyperbolic cosine.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4C
!***TYPE      SINGLE PRECISION (ACOSH-S, DACOSH-D, CACOSH-C)
!***KEYWORDS  ACOSH, ARC HYPERBOLIC COSINE, ELEMENTARY FUNCTIONS, FNLIB,
!             INVERSE HYPERBOLIC COSINE
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! ACOSH(X) computes the arc hyperbolic cosine of X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  ACOSH

  real, parameter :: aln2 = 0.69314718055994530942E+00
  real x
  real, save :: xmax = 0.0E+00

!***FIRST EXECUTABLE STATEMENT  ACOSH

  if ( xmax == 0.0E+00 ) then
    xmax = 1.0E+00 / sqrt ( r1mach(3) )
  end if

  if ( x < 1.0E+00 ) then
    call XERMSG ( 'SLATEC', 'ACOSH', 'X LESS THAN 1', 1, 2 )
  end if

  if ( x < xmax ) then
    acosh = log ( x + sqrt ( x * x - 1.0E+00 ) )
  end if

  if ( x >= xmax ) then
    acosh = aln2 + log ( x )
  end if

  return
end
