function EXPREL (X)
!
!! EXPREL calculates the relative error exponential (EXP(X)-1)/X.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4B
!***TYPE      SINGLE PRECISION (EXPREL-S, DEXPRL-D, CEXPRL-C)
!***KEYWORDS  ELEMENTARY FUNCTIONS, EXPONENTIAL, FIRST ORDER, FNLIB
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate  EXPREL(X) = (EXP(X) - 1.0) / X.   For small ABS(X) the
! Taylor series is used.  If X is negative, the reflection formula
!         EXPREL(X) = EXP(X) * EXPREL(ABS(X))
! may be used.  This reflection formula will be of use when the
! evaluation for small ABS(X) is done by Chebyshev series rather than
! Taylor series.  EXPREL and X are single precision.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH
!***REVISION HISTORY  (YYMMDD)
!   770801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  EXPREL
  LOGICAL FIRST
  SAVE NTERMS, XBND, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  EXPREL
  if (FIRST) THEN
     ALNEPS = LOG(R1MACH(3))
     XN = 3.72 - 0.3*ALNEPS
     XLN = LOG((XN+1.0)/1.36)
     NTERMS = XN - (XN*XLN+ALNEPS)/(XLN+1.36) + 1.5
     XBND = R1MACH(3)
  end if
  FIRST = .FALSE.
!
  ABSX = ABS(X)
  if (ABSX > 0.5) EXPREL = (EXP(X) - 1.0) / X
  if (ABSX > 0.5) RETURN
!
  EXPREL = 1.0
  if (ABSX < XBND) RETURN
!
  EXPREL = 0.0
  DO 20 I=1,NTERMS
    EXPREL = 1.0 + EXPREL*X/(NTERMS+2-I)
 20   CONTINUE
!
  return
end
