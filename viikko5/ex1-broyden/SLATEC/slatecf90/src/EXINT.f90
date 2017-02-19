subroutine EXINT (X, N, KODE, M, TOL, EN, NZ, IERR)
!
!! EXINT computes an M member sequence of exponential integrals ...
!            E(N+K,X), K=0,1,...,M-1 for N  >=  1 and X  >=  0.
!***LIBRARY   SLATEC
!***CATEGORY  C5
!***TYPE      SINGLE PRECISION (EXINT-S, DEXINT-D)
!***KEYWORDS  EXPONENTIAL INTEGRAL, SPECIAL FUNCTIONS
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!         EXINT computes M member sequences of exponential integrals
!         E(N+K,X), K=0,1,...,M-1 for N  >=  1 and X  >=  0.  The
!         exponential integral is defined by
!
!         E(N,X)=integral on (1,infinity) of EXP(-XT)/T**N
!
!         where X=0.0 and N=1 cannot occur simultaneously.  Formulas
!         and notation are found in the NBS Handbook of Mathematical
!         Functions (ref. 1).
!
!         The power series is implemented for X  <=  XCUT and the
!         confluent hypergeometric representation
!
!                     E(A,X) = EXP(-X)*(X**(A-1))*U(A,A,X)
!
!         is computed for X  >  XCUT.  Since sequences are computed in
!         a stable fashion by recurring away from X, A is selected as
!         the integer closest to X within the constraint N  <=  A  <=
!         N+M-1.  For the U computation, A is further modified to be the
!         nearest even integer.  Indices are carried forward or
!         backward by the two term recursion relation
!
!                     K*E(K+1,X) + X*E(K,X) = EXP(-X)
!
!         once E(A,X) is computed.  The U function is computed by means
!         of the backward recursive Miller algorithm applied to the
!         three term contiguous relation for U(A+K,A,X), K=0,1,...
!         This produces accurate ratios and determines U(A+K,A,X), and
!         hence E(A,X), to within a multiplicative constant C.
!         Another contiguous relation applied to C*U(A,A,X) and
!         C*U(A+1,A,X) gets C*U(A+1,A+1,X), a quantity proportional to
!         E(A+1,X).  The normalizing constant C is obtained from the
!         two term recursion relation above with K=A.
!
!     Description of Arguments
!
!         Input
!           X       X  >  0.0 for N=1 and  X  >=  0.0 for N  >=  2
!           N       order of the first member of the sequence, N  >=  1
!                   (X=0.0 and N=1 is an error)
!           KODE    a selection parameter for scaled values
!                   KODE=1   returns        E(N+K,X), K=0,1,...,M-1.
!                       =2   returns EXP(X)*E(N+K,X), K=0,1,...,M-1.
!           M       number of exponential integrals in the sequence,
!                   M  >=  1
!           TOL     relative accuracy wanted, ETOL  <=  TOL  <=  0.1
!                   ETOL = single precision unit roundoff = R1MACH(4)
!
!         Output
!           EN      a vector of dimension at least M containing values
!                   EN(K) = E(N+K-1,X) or EXP(X)*E(N+K-1,X), K=1,M
!                   depending on KODE
!           NZ      underflow indicator
!                   NZ=0   a normal return
!                   NZ=M   X exceeds XLIM and an underflow occurs.
!                          EN(K)=0.0E0 , K=1,M returned on KODE=1
!           IERR    error flag
!                   IERR=0, normal return, computation completed
!                   IERR=1, input error,   no computation
!                   IERR=2, error,         no computation
!                           algorithm termination condition not met
!
!***REFERENCES  M. Abramowitz and I. A. Stegun, Handbook of
!                 Mathematical Functions, NBS AMS Series 55, U.S. Dept.
!                 of Commerce, 1955.
!               D. E. Amos, Computation of exponential integrals, ACM
!                 Transactions on Mathematical Software 6, (1980),
!                 pp. 365-377 and pp. 420-428.
!***ROUTINES CALLED  I1MACH, PSIXN, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   800501  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   910408  Updated the REFERENCES section.  (WRB)
!   920207  Updated with code with a revision date of 880811 from
!           D. Amos.  Included correction of argument list.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  EXINT
  REAL             A,AA,AAMS,AH,AK,AT,B,BK,BT,CC,CNORM,CT,EM,EMX,EN, &
                   ETOL,FNM,FX,PT,P1,P2,S,TOL,TX,X,XCUT,XLIM,XTOL,Y, &
                   YT,Y1,Y2
  REAL             R1MACH,PSIXN
  INTEGER I,IC,ICASE,ICT,IERR,IK,IND,IX,I1M,JSET,K,KK,KN,KODE,KS,M, &
          ML,MU,N,ND,NM,NZ
  INTEGER I1MACH
  DIMENSION EN(*), A(99), B(99), Y(2)
!***FIRST EXECUTABLE STATEMENT  EXINT
  IERR = 0
  NZ = 0
  ETOL = MAX(R1MACH(4),0.5E-18)
  if (X < 0.0E0) IERR = 1
  if (N < 1) IERR = 1
  if (KODE < 1 .OR. KODE > 2) IERR = 1
  if (M < 1) IERR = 1
  if (TOL < ETOL .OR. TOL > 0.1E0) IERR = 1
  if (X == 0.0E0 .AND. N == 1) IERR = 1
  if (IERR /= 0) RETURN
  I1M = -I1MACH(12)
  PT = 2.3026E0*R1MACH(5)*I1M
  XLIM = PT - 6.907755E0
  BT = PT + (N+M-1)
  if (BT > 1000.0E0) XLIM = PT - LOG(BT)
!
  XCUT = 2.0E0
  if (ETOL > 2.0E-7) XCUT = 1.0E0
  if (X > XCUT) go to 100
  if (X == 0.0E0 .AND. N > 1) go to 80
!-----------------------------------------------------------------------
!     SERIES FOR E(N,X) FOR X <= XCUT
!-----------------------------------------------------------------------
  TX = X + 0.5E0
  IX = TX
!-----------------------------------------------------------------------
!     ICASE=1 MEANS INTEGER CLOSEST TO X IS 2 AND N=1
!     ICASE=2 MEANS INTEGER CLOSEST TO X IS 0,1, OR 2 AND N >= 2
!-----------------------------------------------------------------------
  ICASE = 2
  if (IX > N) ICASE = 1
  NM = N - ICASE + 1
  ND = NM + 1
  IND = 3 - ICASE
  MU = M - IND
  ML = 1
  KS = ND
  FNM = NM
  S = 0.0E0
  XTOL = 3.0E0*TOL
  if (ND == 1) go to 10
  XTOL = 0.3333E0*TOL
  S = 1.0E0/FNM
   10 CONTINUE
  AA = 1.0E0
  AK = 1.0E0
  IC = 35
  if (X < ETOL) IC = 1
  DO 50 I=1,IC
    AA = -AA*X/AK
    if (I == NM) go to 30
    S = S - AA/(AK-FNM)
    if (ABS(AA) <= XTOL*ABS(S)) go to 20
    AK = AK + 1.0E0
    go to 50
   20   CONTINUE
    if (I < 2) go to 40
    if (ND-2 > I .OR. I > ND-1) go to 60
    AK = AK + 1.0E0
    go to 50
   30   S = S + AA*(-LOG(X)+PSIXN(ND))
    XTOL = 3.0E0*TOL
   40   AK = AK + 1.0E0
   50 CONTINUE
  if (IC /= 1) go to 340
   60 if (ND == 1) S = S + (-LOG(X)+PSIXN(1))
  if (KODE == 2) S = S*EXP(X)
  EN(1) = S
  EMX = 1.0E0
  if (M == 1) go to 70
  EN(IND) = S
  AA = KS
  if (KODE == 1) EMX = EXP(-X)
  go to (220, 240), ICASE
   70 if (ICASE == 2) RETURN
  if (KODE == 1) EMX = EXP(-X)
  EN(1) = (EMX-S)/X
  return
   80 CONTINUE
  DO 90 I=1,M
    EN(I) = 1.0E0/(N+I-2)
   90 CONTINUE
  return
!-----------------------------------------------------------------------
!     BACKWARD RECURSIVE MILLER ALGORITHM FOR
!              E(N,X)=EXP(-X)*(X**(N-1))*U(N,N,X)
!     WITH RECURSION AWAY FROM N=INTEGER CLOSEST TO X.
!     U(A,B,X) IS THE SECOND CONFLUENT HYPERGEOMETRIC FUNCTION
!-----------------------------------------------------------------------
  100 CONTINUE
  EMX = 1.0E0
  if (KODE == 2) go to 130
  if (X <= XLIM) go to 120
  NZ = M
  DO 110 I=1,M
    EN(I) = 0.0E0
  110 CONTINUE
  return
  120 EMX = EXP(-X)
  130 CONTINUE
  IX = X+0.5E0
  KN = N + M - 1
  if (KN <= IX) go to 140
  if (N < IX .AND. IX < KN) go to 170
  if (N >= IX) go to 160
  go to 340
  140 ICASE = 1
  KS = KN
  ML = M - 1
  MU = -1
  IND = M
  if (KN > 1) go to 180
  150 KS = 2
  ICASE = 3
  go to 180
  160 ICASE = 2
  IND = 1
  KS = N
  MU = M - 1
  if (N > 1) go to 180
  if (KN == 1) go to 150
  IX = 2
  170 ICASE = 1
  KS = IX
  ML = IX - N
  IND = ML + 1
  MU = KN - IX
  180 CONTINUE
  IK = KS/2
  AH = IK
  JSET = 1 + KS - (IK+IK)
!-----------------------------------------------------------------------
!     START COMPUTATION FOR
!              EN(IND) = C*U( A , A ,X)    JSET=1
!              EN(IND) = C*U(A+1,A+1,X)    JSET=2
!     FOR AN EVEN INTEGER A.
!-----------------------------------------------------------------------
  IC = 0
  AA = AH + AH
  AAMS = AA - 1.0E0
  AAMS = AAMS*AAMS
  TX = X + X
  FX = TX + TX
  AK = AH
  XTOL = TOL
  if (TOL <= 1.0E-3) XTOL = 20.0E0*TOL
  CT = AAMS + FX*AH
  EM = (AH+1.0E0)/((X+AA)*XTOL*SQRT(CT))
  BK = AA
  CC = AH*AH
!-----------------------------------------------------------------------
!     FORWARD RECURSION FOR P(IC),P(IC+1) AND INDEX IC FOR BACKWARD
!     RECURSION
!-----------------------------------------------------------------------
  P1 = 0.0E0
  P2 = 1.0E0
  190 CONTINUE
  if (IC == 99) go to 340
  IC = IC + 1
  AK = AK + 1.0E0
  AT = BK/(BK+AK+CC+IC)
  BK = BK + AK + AK
  A(IC) = AT
  BT = (AK+AK+X)/(AK+1.0E0)
  B(IC) = BT
  PT = P2
  P2 = BT*P2 - AT*P1
  P1 = PT
  CT = CT + FX
  EM = EM*AT*(1.0E0-TX/CT)
  if (EM*(AK+1.0E0) > P1*P1) go to 190
  ICT = IC
  KK = IC + 1
  BT = TX/(CT+FX)
  Y2 = (BK/(BK+CC+KK))*(P1/P2)*(1.0E0-BT+0.375E0*BT*BT)
  Y1 = 1.0E0
!-----------------------------------------------------------------------
!     BACKWARD RECURRENCE FOR
!              Y1=             C*U( A ,A,X)
!              Y2= C*(A/(1+A/2))*U(A+1,A,X)
!-----------------------------------------------------------------------
  DO 200 K=1,ICT
    KK = KK - 1
    YT = Y1
    Y1 = (B(KK)*Y1-Y2)/A(KK)
    Y2 = YT
  200 CONTINUE
!-----------------------------------------------------------------------
!     THE CONTIGUOUS RELATION
!              X*U(B,C+1,X)=(C-B)*U(B,C,X)+U(B-1,C,X)
!     WITH  B=A+1 , C=A IS USED FOR
!              Y(2) = C * U(A+1,A+1,X)
!     X IS INCORPORATED INTO THE NORMALIZING RELATION
!-----------------------------------------------------------------------
  PT = Y2/Y1
  CNORM = 1.0E0 - PT*(AH+1.0E0)/AA
  Y(1) = 1.0E0/(CNORM*AA+X)
  Y(2) = CNORM*Y(1)
  if (ICASE == 3) go to 210
  EN(IND) = EMX*Y(JSET)
  if (M == 1) RETURN
  AA = KS
  go to (220, 240), ICASE
!-----------------------------------------------------------------------
!     RECURSION SECTION  N*E(N+1,X) + X*E(N,X)=EMX
!-----------------------------------------------------------------------
  210 EN(1) = EMX*(1.0E0-Y(1))/X
  return
  220 K = IND - 1
  DO 230 I=1,ML
    AA = AA - 1.0E0
    EN(K) = (EMX-AA*EN(K+1))/X
    K = K - 1
  230 CONTINUE
  if (MU <= 0) RETURN
  AA = KS
  240 K = IND
  DO 250 I=1,MU
    EN(K+1) = (EMX-X*EN(K))/AA
    AA = AA + 1.0E0
    K = K + 1
  250 CONTINUE
  return
  340 CONTINUE
  IERR = 2
  return
end
