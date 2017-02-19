function GAMR (X)
!
!! GAMR computes the reciprocal of the Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      SINGLE PRECISION (GAMR-S, DGAMR-D, CGAMR-C)
!***KEYWORDS  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! GAMR is a single precision function that evaluates the reciprocal
! of the gamma function for single precision argument X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  ALGAMS, GAMMA, XERCLR, XGETF, XSETF
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  GAMR
  EXTERNAL GAMMA
!***FIRST EXECUTABLE STATEMENT  GAMR
  GAMR = 0.0
  if (X <= 0.0 .AND. AINT(X) == X) RETURN
!
  call XGETF (IROLD)
  call XSETF (1)
  if (ABS(X) > 10.0) go to 10
  GAMR = 1.0/GAMMA(X)
  call XERCLR
  call XSETF (IROLD)
  return
!
 10   call ALGAMS (X, ALNGX, SGNGX)
  call XERCLR
  call XSETF (IROLD)
  GAMR = SGNGX * EXP(-ALNGX)
  return
!
end
