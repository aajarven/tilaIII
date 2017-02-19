FUNCTION CLBETA (A, B)
!
!! CLBETA computes the natural logarithm of the complete Beta function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7B
!***TYPE      COMPLEX (ALBETA-S, DLBETA-D, CLBETA-C)
!***KEYWORDS  FNLIB, LOGARITHM OF THE COMPLETE BETA FUNCTION,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CLBETA computes the natural log of the complex valued complete beta
! function of complex parameters A and B.  This is a preliminary version
! which is not accurate.
!
! Input Parameters:
!       A   complex and the real part of A positive
!       B   complex and the real part of B positive
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CLNGAM, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  CLBETA
  COMPLEX CLBETA
  COMPLEX A, B, CLNGAM
!***FIRST EXECUTABLE STATEMENT  CLBETA
  if (REAL(A)  <=  0.0 .OR. REAL(B)  <=  0.0) call XERMSG ('SLATEC', &
     'CLBETA', 'REAL PART OF BOTH ARGUMENTS MUST BE GT 0', 1, 2)
!
  CLBETA = CLNGAM(A) + CLNGAM(B) - CLNGAM(A+B)
!
  return
end
