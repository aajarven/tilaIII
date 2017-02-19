function BINOM (N, M)
!
!! BINOM computes the binomial coefficients.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C1
!***TYPE      SINGLE PRECISION (BINOM-S, DBINOM-D)
!***KEYWORDS  BINOMIAL COEFFICIENTS, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BINOM(N,M) calculates the binomial coefficient (N!)/((M!)*(N-M)!).
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  ALNREL, R1MACH, R9LGMC, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  BINOM
  LOGICAL FIRST
  SAVE SQ2PIL, BILNMX, FINTMX, FIRST
  DATA SQ2PIL / 0.91893853320467274E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BINOM
  if (FIRST) THEN
     BILNMX = LOG (R1MACH(2))
     FINTMX = 0.9/R1MACH(3)
  end if
  FIRST = .FALSE.
!
  if (N  <  0 .OR. M  <  0) call XERMSG ('SLATEC', 'BINOM', &
     'N OR M LT ZERO', 1, 2)
  if (N  <  M) call XERMSG ('SLATEC', 'BINOM', 'N LT M', 2, 2)
!
  K = MIN (M, N-M)
  if (K > 20) go to 30
  if (K*LOG(AMAX0(N,1)) > BILNMX) go to 30
!
  BINOM = 1.
  if (K == 0) RETURN
!
  DO 20 I=1,K
    BINOM = BINOM * REAL(N-I+1)/I
 20   CONTINUE
!
  if (BINOM < FINTMX) BINOM = AINT (BINOM+0.5)
  return
!
! if K < 9, APPROX IS NOT VALID AND ANSWER IS CLOSE TO THE OVERFLOW LIM
 30   if (K  <  9) call XERMSG ('SLATEC', 'BINOM', &
     'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG', 3, 2)
!
  XN = N + 1
  XK = K + 1
  XNK = N - K + 1
!
  CORR = R9LGMC(XN) - R9LGMC(XK) - R9LGMC(XNK)
  BINOM = XK*LOG(XNK/XK) - XN*ALNREL(-(XK-1.)/XN) &
    - 0.5*LOG(XN*XNK/XK) + 1.0 - SQ2PIL + CORR
!
  if (BINOM  >  BILNMX) call XERMSG ('SLATEC', 'BINOM', &
     'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG', 3, 2)
!
  BINOM = EXP (BINOM)
  if (BINOM < FINTMX) BINOM = AINT (BINOM+0.5)
!
  return
end
