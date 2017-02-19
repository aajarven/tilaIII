  DOUBLE PRECISION FUNCTION DLBETA (A, B)
!
!! DLBETA computes the natural logarithm of the complete Beta function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7B
!***TYPE      DOUBLE PRECISION (ALBETA-S, DLBETA-D, CLBETA-C)
!***KEYWORDS  FNLIB, LOGARITHM OF THE COMPLETE BETA FUNCTION,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DLBETA(A,B) calculates the double precision natural logarithm of
! the complete beta function for double precision arguments
! A and B.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D9LGMC, DGAMMA, DLNGAM, DLNREL, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  DLBETA
  DOUBLE PRECISION A, B, P, Q, CORR, SQ2PIL, D9LGMC, DGAMMA, DLNGAM, &
    DLNREL
  EXTERNAL DGAMMA
  SAVE SQ2PIL
  DATA SQ2PIL / 0.91893853320467274178032973640562D0 /
!***FIRST EXECUTABLE STATEMENT  DLBETA
  P = MIN (A, B)
  Q = MAX (A, B)
!
  if (P  <=  0.D0) call XERMSG ('SLATEC', 'DLBETA', &
     'BOTH ARGUMENTS MUST BE GT ZERO', 1, 2)
!
  if (P >= 10.D0) go to 30
  if (Q >= 10.D0) go to 20
!
! P AND Q ARE SMALL.
!
  DLBETA = LOG (DGAMMA(P) * (DGAMMA(Q)/DGAMMA(P+Q)) )
  return
!
! P IS SMALL, BUT Q IS BIG.
!
 20   CORR = D9LGMC(Q) - D9LGMC(P+Q)
  DLBETA = DLNGAM(P) + CORR + P - P*LOG(P+Q) &
    + (Q-0.5D0)*DLNREL(-P/(P+Q))
  return
!
! P AND Q ARE BIG.
!
 30   CORR = D9LGMC(P) + D9LGMC(Q) - D9LGMC(P+Q)
  DLBETA = -0.5D0*LOG(Q) + SQ2PIL + CORR + (P-0.5D0)*LOG(P/(P+Q)) &
    + Q*DLNREL(-P/(P+Q))
  return
!
end
