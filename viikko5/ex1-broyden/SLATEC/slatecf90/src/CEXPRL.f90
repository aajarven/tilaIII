FUNCTION CEXPRL (Z)
!
!! CEXPRL calculates the relative error exponential (EXP(X)-1)/X.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4B
!***TYPE      COMPLEX (EXPREL-S, DEXPRL-D, CEXPRL-C)
!***KEYWORDS  ELEMENTARY FUNCTIONS, EXPONENTIAL, FIRST ORDER, FNLIB
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate  (EXP(Z)-1)/Z .  For small ABS(Z), we use the Taylor
! series.  We could instead use the expression
!        CEXPRL(Z) = (EXP(X)*EXP(I*Y)-1)/Z
!                  = (X*EXPREL(X) * (1 - 2*SIN(Y/2)**2) - 2*SIN(Y/2)**2
!                                    + I*SIN(Y)*(1+X*EXPREL(X))) / Z
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH
!***REVISION HISTORY  (YYMMDD)
!   770801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CEXPRL
  COMPLEX CEXPRL
  COMPLEX Z
  LOGICAL FIRST
  SAVE NTERMS, RBND, FIRST
  DATA FIRST / .TRUE. /
!***FIRST EXECUTABLE STATEMENT  CEXPRL
  if (FIRST) THEN
     ALNEPS = LOG(R1MACH(3))
     XN = 3.72 - 0.3*ALNEPS
     XLN = LOG((XN+1.0)/1.36)
     NTERMS = XN - (XN*XLN+ALNEPS)/(XLN+1.36) + 1.5
     RBND = R1MACH(3)
  end if
  FIRST = .FALSE.
!
  R = ABS(Z)
  if (R > 0.5) CEXPRL = (EXP(Z) - 1.0) / Z
  if (R > 0.5) RETURN
!
  CEXPRL = (1.0, 0.0)
  if (R < RBND) RETURN
!
  CEXPRL = (0.0, 0.0)
  DO 20 I=1,NTERMS
    CEXPRL = 1.0 + CEXPRL*Z/(NTERMS+2-I)
 20   CONTINUE
!
  return
end
