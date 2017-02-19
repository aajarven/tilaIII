FUNCTION CGAMR (Z)
!
!! CGAMR computes the reciprocal of the Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      COMPLEX (GAMR-S, DGAMR-D, CGAMR-C)
!***KEYWORDS  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CGAMR(Z) calculates the reciprocal gamma function for COMPLEX
! argument Z.  This is a preliminary version that is not accurate.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CLNGAM, XERCLR, XGETF, XSETF
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CGAMR
  COMPLEX CGAMR
  COMPLEX Z, CLNGAM
!***FIRST EXECUTABLE STATEMENT  CGAMR
  CGAMR = (0.0, 0.0)
  X = REAL (Z)
  if (X <= 0.0 .AND. AINT(X) == X .AND. AIMAG(Z) == 0.0) RETURN
!
  call XGETF (IROLD)
  call XSETF (1)
  CGAMR = CLNGAM(Z)
  call XERCLR
  call XSETF (IROLD)
  CGAMR = EXP (-CGAMR)
!
  return
end
