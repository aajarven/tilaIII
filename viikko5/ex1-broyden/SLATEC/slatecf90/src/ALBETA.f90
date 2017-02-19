function ALBETA (A, B)
!
!! ALBETA computes the natural logarithm of the complete Beta function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7B
!***TYPE      SINGLE PRECISION (ALBETA-S, DLBETA-D, CLBETA-C)
!***KEYWORDS  FNLIB, LOGARITHM OF THE COMPLETE BETA FUNCTION,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! ALBETA computes the natural log of the complete beta function.
!
! Input Parameters:
!       A   real and positive
!       B   real and positive
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  ALNGAM, ALNREL, GAMMA, R9LGMC, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  ALBETA
  EXTERNAL GAMMA
  SAVE SQ2PIL
  DATA SQ2PIL / 0.91893853320467274E0 /
!***FIRST EXECUTABLE STATEMENT  ALBETA
  P = MIN (A, B)
  Q = MAX (A, B)
!
  if (P  <=  0.0) call XERMSG ('SLATEC', 'ALBETA', &
     'BOTH ARGUMENTS MUST BE GT ZERO', 1, 2)
  if (P >= 10.0) go to 30
  if (Q >= 10.0) go to 20
!
! P AND Q ARE SMALL.
!
  ALBETA = LOG(GAMMA(P) * (GAMMA(Q)/GAMMA(P+Q)) )
  return
!
! P IS SMALL, BUT Q IS BIG.
!
 20   CORR = R9LGMC(Q) - R9LGMC(P+Q)
  ALBETA = ALNGAM(P) + CORR + P - P*LOG(P+Q) + &
    (Q-0.5)*ALNREL(-P/(P+Q))
  return
!
! P AND Q ARE BIG.
!
 30   CORR = R9LGMC(P) + R9LGMC(Q) - R9LGMC(P+Q)
  ALBETA = -0.5*LOG(Q) + SQ2PIL + CORR + (P-0.5)*LOG(P/(P+Q)) &
    + Q*ALNREL(-P/(P+Q))
  return
!
end
