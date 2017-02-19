  DOUBLE PRECISION FUNCTION DEXPRL (X)
!
!! DEXPRL calculates the relative error exponential (EXP(X)-1)/X.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4B
!***TYPE      DOUBLE PRECISION (EXPREL-S, DEXPRL-D, CEXPRL-C)
!***KEYWORDS  ELEMENTARY FUNCTIONS, EXPONENTIAL, FIRST ORDER, FNLIB
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate  EXPREL(X) = (EXP(X) - 1.0) / X.   For small ABS(X) the
! Taylor series is used.  If X is negative the reflection formula
!         EXPREL(X) = EXP(X) * EXPREL(ABS(X))
! may be used.  This reflection formula will be of use when the
! evaluation for small ABS(X) is done by Chebyshev series rather than
! Taylor series.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   770801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DEXPRL
  DOUBLE PRECISION X, ABSX, ALNEPS, XBND, XLN, XN,  D1MACH
  LOGICAL FIRST
  SAVE NTERMS, XBND, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DEXPRL
  if (FIRST) THEN
     ALNEPS = LOG(D1MACH(3))
     XN = 3.72D0 - 0.3D0*ALNEPS
     XLN = LOG((XN+1.0D0)/1.36D0)
     NTERMS = XN - (XN*XLN+ALNEPS)/(XLN+1.36D0) + 1.5D0
     XBND = D1MACH(3)
  end if
  FIRST = .FALSE.
!
  ABSX = ABS(X)
  if (ABSX > 0.5D0) DEXPRL = (EXP(X)-1.0D0)/X
  if (ABSX > 0.5D0) RETURN
!
  DEXPRL = 1.0D0
  if (ABSX < XBND) RETURN
!
  DEXPRL = 0.0D0
  DO 20 I=1,NTERMS
    DEXPRL = 1.0D0 + DEXPRL*X/(NTERMS+2-I)
 20   CONTINUE
!
  return
end
