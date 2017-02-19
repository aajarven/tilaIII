FUNCTION DBETAI (X, PIN, QIN)
!
!! DBETAI calculates the incomplete Beta function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7F
!***TYPE      DOUBLE PRECISION (BETAI-S, DBETAI-D)
!***KEYWORDS  FNLIB, INCOMPLETE BETA FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!   DBETAI calculates the DOUBLE PRECISION incomplete beta function.
!
!   The incomplete beta function ratio is the probability that a
!   random variable from a beta distribution having parameters PIN and
!   QIN will be less than or equal to X.
!
!     -- Input Arguments -- All arguments are DOUBLE PRECISION.
!   X      upper limit of integration.  X must be in (0,1) inclusive.
!   PIN    first beta distribution parameter.  PIN must be  >  0.0.
!   QIN    second beta distribution parameter.  QIN must be  >  0.0.
!
!***REFERENCES  Nancy E. Bosten and E. L. Battiste, Remark on Algorithm
!                 179, Communications of the ACM 17, 3 (March 1974),
!                 pp. 156.
!***ROUTINES CALLED  D1MACH, DLBETA, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
!***END PROLOGUE  DBETAI
  DOUBLE PRECISION DBETAI
  DOUBLE PRECISION X, PIN, QIN, ALNEPS, ALNSML, C, EPS, FINSUM, P, &
    PS, Q, SML, TERM, XB, XI, Y, D1MACH, DLBETA, P1
  LOGICAL FIRST
  SAVE EPS, ALNEPS, SML, ALNSML, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBETAI
  if (FIRST) THEN
     EPS = D1MACH(3)
     ALNEPS = LOG (EPS)
     SML = D1MACH(1)
     ALNSML = LOG (SML)
  end if
  FIRST = .FALSE.
!
  if (X  <  0.D0 .OR. X  >  1.D0) call XERMSG ('SLATEC', 'DBETAI', &
     'X IS NOT IN THE RANGE (0,1)', 1, 2)
  if (PIN  <=  0.D0 .OR. QIN  <=  0.D0) call XERMSG ('SLATEC', &
     'DBETAI', 'P AND/OR Q IS LE ZERO', 2, 2)
!
  Y = X
  P = PIN
  Q = QIN
  if (Q <= P .AND. X < 0.8D0) go to 20
  if (X < 0.2D0) go to 20
  Y = 1.0D0 - Y
  P = QIN
  Q = PIN
!
 20   if ((P+Q)*Y/(P+1.D0) < EPS) go to 80
!
! EVALUATE THE INFINITE SUM FIRST.  TERM WILL EQUAL
! Y**P/BETA(PS,P) * (1.-PS)-SUB-I * Y**I / FAC(I) .
!
  PS = Q - AINT(Q)
  if (PS == 0.D0) PS = 1.0D0
  XB = P*LOG(Y) - DLBETA(PS,P) - LOG(P)
  DBETAI = 0.0D0
  if (XB < ALNSML) go to 40
!
  DBETAI = EXP (XB)
  TERM = DBETAI*P
  if (PS == 1.0D0) go to 40
  N = MAX (ALNEPS/LOG(Y), 4.0D0)
  DO 30 I=1,N
    XI = I
    TERM = TERM * (XI-PS)*Y/XI
    DBETAI = DBETAI + TERM/(P+XI)
 30   CONTINUE
!
! NOW EVALUATE THE FINITE SUM, MAYBE.
!
 40   if (Q <= 1.0D0) go to 70
!
  XB = P*LOG(Y) + Q*LOG(1.0D0-Y) - DLBETA(P,Q) - LOG(Q)
  IB = MAX (XB/ALNSML, 0.0D0)
  TERM = EXP(XB - IB*ALNSML)
  C = 1.0D0/(1.D0-Y)
  P1 = Q*C/(P+Q-1.D0)
!
  FINSUM = 0.0D0
  N = Q
  if (Q == DBLE(N)) N = N - 1
  DO 50 I=1,N
    if (P1 <= 1.0D0 .AND. TERM/EPS <= FINSUM) go to 60
    XI = I
    TERM = (Q-XI+1.0D0)*C*TERM/(P+Q-XI)
!
    if (TERM > 1.0D0) IB = IB - 1
    if (TERM > 1.0D0) TERM = TERM*SML
!
    if (IB == 0) FINSUM = FINSUM + TERM
 50   CONTINUE
!
 60   DBETAI = DBETAI + FINSUM
 70   if (Y /= X .OR. P /= PIN) DBETAI = 1.0D0 - DBETAI
  DBETAI = MAX (MIN (DBETAI, 1.0D0), 0.0D0)
  return
!
 80   DBETAI = 0.0D0
  XB = P*LOG(MAX(Y,SML)) - LOG(P) - DLBETA(P,Q)
  if (XB > ALNSML .AND. Y /= 0.0D0) DBETAI = EXP(XB)
  if (Y /= X .OR. P /= PIN) DBETAI = 1.0D0 - DBETAI
!
  return
end
