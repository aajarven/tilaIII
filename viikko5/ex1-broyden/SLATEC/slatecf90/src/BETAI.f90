FUNCTION BETAI (X, PIN, QIN)
!
!! BETAI calculates the incomplete Beta function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7F
!***TYPE      SINGLE PRECISION (BETAI-S, DBETAI-D)
!***KEYWORDS  FNLIB, INCOMPLETE BETA FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!   BETAI calculates the REAL incomplete beta function.
!
!   The incomplete beta function ratio is the probability that a
!   random variable from a beta distribution having parameters PIN and
!   QIN will be less than or equal to X.
!
!     -- Input Arguments -- All arguments are REAL.
!   X      upper limit of integration.  X must be in (0,1) inclusive.
!   PIN    first beta distribution parameter.  PIN must be  >  0.0.
!   QIN    second beta distribution parameter.  QIN must be  >  0.0.
!
!***REFERENCES  Nancy E. Bosten and E. L. Battiste, Remark on Algorithm
!                 179, Communications of the ACM 17, 3 (March 1974),
!                 pp. 156.
!***ROUTINES CALLED  ALBETA, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
!***END PROLOGUE  BETAI
  REAL BETAI
  LOGICAL FIRST
  SAVE EPS, ALNEPS, SML, ALNSML, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  BETAI
  if (FIRST) THEN
     EPS = R1MACH(3)
     ALNEPS = LOG(EPS)
     SML = R1MACH(1)
     ALNSML = LOG(SML)
  end if
  FIRST = .FALSE.
!
  if (X  <  0. .OR. X  >  1.0) call XERMSG ('SLATEC', 'BETAI', &
     'X IS NOT IN THE RANGE (0,1)', 1, 2)
  if (PIN  <=  0. .OR. QIN  <=  0.) call XERMSG ('SLATEC', 'BETAI', &
     'P AND/OR Q IS LE ZERO', 2, 2)
!
  Y = X
  P = PIN
  Q = QIN
  if (Q <= P .AND. X < 0.8) go to 20
  if (X < 0.2) go to 20
  Y = 1.0 - Y
  P = QIN
  Q = PIN
!
 20   if ((P+Q)*Y/(P+1.) < EPS) go to 80
!
! EVALUATE THE INFINITE SUM FIRST.
! TERM WILL EQUAL Y**P/BETA(PS,P) * (1.-PS)I * Y**I / FAC(I)
!
  PS = Q - AINT(Q)
  if (PS == 0.) PS = 1.0
  XB = P*LOG(Y) -  ALBETA(PS, P) - LOG(P)
  BETAI = 0.0
  if (XB < ALNSML) go to 40
!
  BETAI = EXP (XB)
  TERM = BETAI*P
  if (PS == 1.0) go to 40
!
  N = MAX (ALNEPS/LOG(Y), 4.0E0)
  DO 30 I=1,N
    TERM = TERM*(I-PS)*Y/I
    BETAI = BETAI + TERM/(P+I)
 30   CONTINUE
!
! NOW EVALUATE THE FINITE SUM, MAYBE.
!
 40   if (Q <= 1.0) go to 70
!
  XB = P*LOG(Y) + Q*LOG(1.0-Y) - ALBETA(P,Q) - LOG(Q)
  IB = MAX (XB/ALNSML, 0.0E0)
  TERM = EXP (XB - IB*ALNSML)
  C = 1.0/(1.0-Y)
  P1 = Q*C/(P+Q-1.)
!
  FINSUM = 0.0
  N = Q
  if (Q == REAL(N)) N = N - 1
  DO 50 I=1,N
    if (P1 <= 1.0 .AND. TERM/EPS <= FINSUM) go to 60
    TERM = (Q-I+1)*C*TERM/(P+Q-I)
!
    if (TERM > 1.0) IB = IB - 1
    if (TERM > 1.0) TERM = TERM*SML
!
    if (IB == 0) FINSUM = FINSUM + TERM
 50   CONTINUE
!
 60   BETAI = BETAI + FINSUM
 70   if (Y /= X .OR. P /= PIN) BETAI = 1.0 - BETAI
  BETAI = MAX (MIN (BETAI, 1.0), 0.0)
  return
!
 80   BETAI = 0.0
  XB = P*LOG(MAX(Y,SML)) - LOG(P) - ALBETA(P,Q)
  if (XB > ALNSML .AND. Y /= 0.) BETAI = EXP (XB)
  if (Y /= X .OR. P /= PIN) BETAI = 1.0 - BETAI
  return

end
