FUNCTION CLNREL (Z)
!
!! CLNREL evaluates ln(1+X) accurate in the sense of relative error.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4B
!***TYPE      COMPLEX (ALNREL-S, DLNREL-D, CLNREL-C)
!***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CLNREL(Z) = LOG(1+Z) with relative error accuracy near Z = 0.
! Let   RHO = ABS(Z)  and
!       R**2 = ABS(1+Z)**2 = (1+X)**2 + Y**2 = 1 + 2*X + RHO**2 .
! Now if RHO is small we may evaluate CLNREL(Z) accurately by
!       LOG(1+Z) = CMPLX  (LOG(R), CARG(1+Z))
!                 = CMPLX  (0.5*LOG(R**2), CARG(1+Z))
!                 = CMPLX  (0.5*ALNREL(2*X+RHO**2), CARG(1+Z))
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  ALNREL, CARG, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  CLNREL
  COMPLEX CLNREL
  COMPLEX Z
  SAVE SQEPS
  DATA SQEPS /0.0/
!***FIRST EXECUTABLE STATEMENT  CLNREL
  if (SQEPS == 0.) SQEPS = SQRT (R1MACH(4))
!
  if (ABS(1.+Z)  <  SQEPS) call XERMSG ('SLATEC', 'CLNREL', &
     'ANSWER LT HALF PRECISION BECAUSE Z TOO NEAR -1', 1, 1)
!
  RHO = ABS(Z)
  if (RHO > 0.375) CLNREL = LOG (1.0+Z)
  if (RHO > 0.375) RETURN
!
  X = REAL(Z)
  CLNREL = CMPLX (0.5*ALNREL(2.*X+RHO**2), CARG(1.0+Z))
!
  return
end
