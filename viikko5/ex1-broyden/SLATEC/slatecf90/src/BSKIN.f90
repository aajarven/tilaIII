subroutine BSKIN (X, N, KODE, M, Y, NZ, IERR)
!
!! BSKIN computes repeated integrals of the K-zero Bessel function.
!
!***LIBRARY   SLATEC
!***CATEGORY  C10F
!***TYPE      SINGLE PRECISION (BSKIN-S, DBSKIN-D)
!***KEYWORDS  BICKLEY FUNCTIONS, EXPONENTIAL INTEGRAL,
!             INTEGRALS OF BESSEL FUNCTIONS, K-ZERO BESSEL FUNCTION
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!         The following definitions are used in BSKIN:
!
!   Definition 1
!         KI(0,X) = K-zero Bessel function.
!
!   Definition 2
!         KI(N,X) = Bickley Function
!                 =  integral from X to infinity of KI(N-1,t)dt
!                     for X .ge. 0 and N = 1,2,...
!   ____________________________________________________________________
!      BSKIN computes sequences of Bickley functions (repeated integrals
!      of the K0 Bessel function); i.e. for fixed X and N and K=1,...,
!      BSKIN computes the M-member sequence
!
!                     Y(K) =        KI(N+K-1,X) for KODE=1
!      or
!                     Y(K) = EXP(X)*KI(N+K-1,X) for KODE=2,
!
!      for N.ge.0 and X.ge.0 (N and X cannot be zero simultaneously).
!
!      INPUT
!        X      - Argument, X .ge. 0.0E0
!        N      - Order of first member of the sequence N .ge. 0
!        KODE   - Selection parameter
!                 KODE = 1 returns Y(K)=       KI(N+K-1,X), K=1,M
!                      = 2 returns Y(K)=EXP(X)*KI(N+K-1,X), K=1,M
!        M      - Number of members in the sequence, M.ge.1
!
!      OUTPUT
!        Y      - A vector of dimension at least M containing the
!                 sequence selected by KODE.
!        NZ     - Underflow flag
!                 NZ = 0 means computation completed
!                    = M means an exponential underflow occurred on
!                        KODE=1.  Y(K)=0.0E0, K=1,...,M is returned
!        IERR   - Error flag
!                 IERR = 0, Normal return, computation completed.
!                      = 1, Input error,   no computation.
!                      = 2, Error,         no computation.  The
!                           termination condition was not met.
!
!      The nominal computational accuracy is the maximum of unit
!      roundoff (=R1MACH(4)) and 1.0e-18 since critical constants
!      are given to only 18 digits.
!
!      DBSKIN is the double precision version of BSKIN.
!
! *Long Description:
!
!         Numerical recurrence on
!
!      (L-1)*KI(L,X) = X(KI(L-3,X) - KI(L-1,X)) + (L-2)*KI(L-2,X)
!
!         is stable where recurrence is carried forward or backward
!         away from INT(X+0.5).  The power series for indices 0,1 and 2
!         on 0.le.X.le. 2 starts a stable recurrence for indices
!         greater than 2.  If N is sufficiently large (N.gt.NLIM), the
!         uniform asymptotic expansion for N to INFINITY is more
!         economical.  On X.gt.2 the recursion is started by evaluating
!         the uniform expansion for the three members whose indices are
!         closest to INT(X+0.5) within the set N,...,N+M-1.  Forward
!         recurrence, backward recurrence or both, complete the
!         sequence depending on the relation of INT(X+0.5) to the
!         indices N,...,N+M-1.
!
!***REFERENCES  D. E. Amos, Uniform asymptotic expansions for
!                 exponential integrals E(N,X) and Bickley functions
!                 KI(N,X), ACM Transactions on Mathematical Software,
!                 1983.
!               D. E. Amos, A portable Fortran subroutine for the
!                 Bickley functions KI(N,X), Algorithm 609, ACM
!                 Transactions on Mathematical Software, 1983.
!***ROUTINES CALLED  BKIAS, BKISR, EXINT, GAMRN, I1MACH, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   820601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891009  Removed unreferenced statement label.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BSKIN
  INTEGER I, ICASE, IERR, IL, I1M, K, KK, KODE, KTRMS, M, &
   M3, N, NE, NFLG, NL, NLIM, NN, NP, NS, NT, NZ
  INTEGER I1MACH
  REAL A, ENLIM, EXI, FN, GR, H, HN, HRTPI, SS, TOL, T1, T2, W, X, &
   XLIM, XNLIM, XP, Y, YS, YSS
  REAL GAMRN, R1MACH
  DIMENSION EXI(102), A(50), YS(3), YSS(3), H(31), Y(*)
  SAVE A, HRTPI
!-----------------------------------------------------------------------
!             COEFFICIENTS IN SERIES OF EXPONENTIAL INTEGRALS
!-----------------------------------------------------------------------
  DATA A(1), A(2), A(3), A(4), A(5), A(6), A(7), A(8), A(9), A(10), &
   A(11), A(12), A(13), A(14), A(15), A(16), A(17), A(18), A(19), &
   A(20), A(21), A(22), A(23), A(24) /1.00000000000000000E+00, &
   5.00000000000000000E-01,3.75000000000000000E-01, &
   3.12500000000000000E-01,2.73437500000000000E-01, &
   2.46093750000000000E-01,2.25585937500000000E-01, &
   2.09472656250000000E-01,1.96380615234375000E-01, &
   1.85470581054687500E-01,1.76197052001953125E-01, &
   1.68188095092773438E-01,1.61180257797241211E-01, &
   1.54981017112731934E-01,1.49445980787277222E-01, &
   1.44464448094367981E-01,1.39949934091418982E-01, &
   1.35833759559318423E-01,1.32060599571559578E-01, &
   1.28585320635465905E-01,1.25370687619579257E-01, &
   1.22385671247684513E-01,1.19604178719328047E-01, &
   1.17004087877603524E-01/
  DATA A(25), A(26), A(27), A(28), A(29), A(30), A(31), A(32), &
   A(33), A(34), A(35), A(36), A(37), A(38), A(39), A(40), A(41), &
   A(42), A(43), A(44), A(45), A(46), A(47), A(48) &
   /1.14566502713486784E-01,1.12275172659217048E-01, &
   1.10116034723462874E-01,1.08076848895250599E-01, &
   1.06146905164978267E-01,1.04316786110409676E-01, &
   1.02578173008569515E-01,1.00923686347140974E-01, &
   9.93467537479668965E-02,9.78414999033007314E-02, &
   9.64026543164874854E-02,9.50254735405376642E-02, &
   9.37056752969190855E-02,9.24393823875012600E-02, &
   9.12230747245078224E-02,9.00535481254756708E-02, &
   8.89278787739072249E-02,8.78433924473961612E-02, &
   8.67976377754033498E-02,8.57883629175498224E-02, &
   8.48134951571231199E-02,8.38711229887106408E-02, &
   8.29594803475290034E-02,8.20769326842574183E-02/
  DATA A(49), A(50) /8.12219646354630702E-02,8.03931690779583449E-02 &
   /
!-----------------------------------------------------------------------
!             SQRT(PI)/2
!-----------------------------------------------------------------------
  DATA HRTPI /8.86226925452758014E-01/
!
!***FIRST EXECUTABLE STATEMENT  BSKIN
  IERR = 0
  NZ=0
  if (X < 0.0E0) IERR=1
  if (N < 0) IERR=1
  if (KODE < 1 .OR. KODE > 2) IERR=1
  if (M < 1) IERR=1
  if (X == 0.0E0 .AND. N == 0) IERR=1
  if (IERR /= 0) RETURN
  if (X == 0.0E0) go to 300
  I1M = -I1MACH(12)
  T1 = 2.3026E0*R1MACH(5)*I1M
  XLIM = T1 - 3.228086E0
  T2 = T1 + N + M - 1
  if (T2 > 1000.0E0) XLIM = T1 - 0.5E0*(LOG(T2)-0.451583E0)
  if (X > XLIM .AND. KODE == 1) go to 320
  TOL = MAX(R1MACH(4),1.0E-18)
  I1M = I1MACH(11)
!-----------------------------------------------------------------------
!     LN(NLIM) = 0.125*LN(EPS),   NLIM = 2*KTRMS+N
!-----------------------------------------------------------------------
  XNLIM = 0.287823E0*(I1M-1)*R1MACH(5)
  ENLIM = EXP(XNLIM)
  NLIM = INT(ENLIM) + 2
  NLIM = MIN(100,NLIM)
  NLIM = MAX(20,NLIM)
  M3 = MIN(M,3)
  NL = N + M - 1
  if (X > 2.0E0) go to 130
  if (N > NLIM) go to 280
!-----------------------------------------------------------------------
!     COMPUTATION BY SERIES FOR 0 <= X <= 2
!-----------------------------------------------------------------------
  NFLG = 0
  NN = N
  if (NL <= 2) go to 60
  M3 = 3
  NN = 0
  NFLG = 1
   60 CONTINUE
  XP = 1.0E0
  if (KODE == 2) XP = EXP(X)
  DO 80 I=1,M3
    call BKISR(X, NN, W, IERR)
  if ( IERR /= 0) RETURN
    W = W*XP
    if (NN < N) go to 70
    KK = NN - N + 1
    Y(KK) = W
   70   CONTINUE
    YS(I) = W
    NN = NN + 1
   80 CONTINUE
  if (NFLG == 0) RETURN
  NS = NN
  XP = 1.0E0
   90 CONTINUE
!-----------------------------------------------------------------------
!     FORWARD RECURSION SCALED BY EXP(X) ON ICASE=0,1,2
!-----------------------------------------------------------------------
  FN = NS - 1
  IL = NL - NS + 1
  if (IL <= 0) RETURN
  DO 110 I=1,IL
    T1 = YS(2)
    T2 = YS(3)
    YS(3) = (X*(YS(1)-YS(3))+(FN-1.0E0)*YS(2))/FN
    YS(2) = T2
    YS(1) = T1
    FN = FN + 1.0E0
    if (NS < N) go to 100
    KK = NS - N + 1
    Y(KK) = YS(3)*XP
  100   CONTINUE
    NS = NS + 1
  110 CONTINUE
  return
!-----------------------------------------------------------------------
!     COMPUTATION BY ASYMPTOTIC EXPANSION FOR X > 2
!-----------------------------------------------------------------------
  130 CONTINUE
  W = X + 0.5E0
  NT = INT(W)
  if (NL > NT) go to 270
!-----------------------------------------------------------------------
!     CASE NL <= NT, ICASE=0
!-----------------------------------------------------------------------
  ICASE = 0
  NN = NL
  NFLG = MIN(M-M3,1)
  140 CONTINUE
  KK = (NLIM-NN)/2
  KTRMS = MAX(0,KK)
  NS = NN + 1
  NP = NN - M3 + 1
  XP = 1.0E0
  if (KODE == 1) XP = EXP(-X)
  DO 150 I=1,M3
    KK = I
    call BKIAS(X, NP, KTRMS, A, W, KK, NE, GR, H, IERR)
  if ( IERR /= 0) RETURN
    YS(I) = W
    NP = NP + 1
  150 CONTINUE
!-----------------------------------------------------------------------
!     SUM SERIES OF EXPONENTIAL INTEGRALS BACKWARD
!-----------------------------------------------------------------------
  if (KTRMS == 0) go to 160
  NE = KTRMS + KTRMS + 1
  NP = NN - M3 + 2
  call EXINT(X, NP, 2, NE, TOL, EXI, NZ, IERR)
  if ( NZ /= 0) go to 320
  if ( IERR == 2) RETURN
  160 CONTINUE
  DO 190 I=1,M3
    SS = 0.0E0
    if (KTRMS == 0) go to 180
    KK = I + KTRMS + KTRMS - 2
    IL = KTRMS
    DO 170 K=1,KTRMS
      SS = SS + A(IL)*EXI(KK)
      KK = KK - 2
      IL = IL - 1
  170   CONTINUE
  180   CONTINUE
    YS(I) = YS(I) + SS
  190 CONTINUE
  if (ICASE == 1) go to 200
  if (NFLG /= 0) go to 220
  200 CONTINUE
  DO 210 I=1,M3
    Y(I) = YS(I)*XP
  210 CONTINUE
  if (ICASE == 1 .AND. NFLG == 1) go to 90
  return
  220 CONTINUE
!-----------------------------------------------------------------------
!     BACKWARD RECURSION SCALED BY EXP(X) ICASE=0,2
!-----------------------------------------------------------------------
  KK = NN - N + 1
  K = M3
  DO 230 I=1,M3
    Y(KK) = YS(K)*XP
    YSS(I) = YS(I)
    KK = KK - 1
    K = K - 1
  230 CONTINUE
  IL = KK
  if (IL <= 0) go to 250
  FN = NN - 3
  DO 240 I=1,IL
    T1 = YS(2)
    T2 = YS(1)
    YS(1) = YS(2) + ((FN+2.0E0)*YS(3)-(FN+1.0E0)*YS(1))/X
    YS(2) = T2
    YS(3) = T1
    Y(KK) = YS(1)*XP
    KK = KK - 1
    FN = FN - 1.0E0
  240 CONTINUE
  250 CONTINUE
  if (ICASE /= 2) RETURN
  DO 260 I=1,M3
    YS(I) = YSS(I)
  260 CONTINUE
  go to 90
  270 CONTINUE
  if (N < NT) go to 290
!-----------------------------------------------------------------------
!     ICASE=1, NT <= N <= NL WITH FORWARD RECURSION
!-----------------------------------------------------------------------
  280 CONTINUE
  NN = N + M3 - 1
  NFLG = MIN(M-M3,1)
  ICASE = 1
  go to 140
!-----------------------------------------------------------------------
!     ICASE=2, N < NT < NL WITH BOTH FORWARD AND BACKWARD RECURSION
!-----------------------------------------------------------------------
  290 CONTINUE
  NN = NT + 1
  NFLG = MIN(M-M3,1)
  ICASE = 2
  go to 140
!-----------------------------------------------------------------------
!     X=0 CASE
!-----------------------------------------------------------------------
  300 CONTINUE
  FN = N
  HN = 0.5E0*FN
  GR = GAMRN(HN)
  Y(1) = HRTPI*GR
  if (M == 1) RETURN
  Y(2) = HRTPI/(HN*GR)
  if (M == 2) RETURN
  DO 310 K=3,M
    Y(K) = FN*Y(K-2)/(FN+1.0E0)
    FN = FN + 1.0E0
  310 CONTINUE
  return
!-----------------------------------------------------------------------
!     UNDERFLOW ON KODE=1, X > XLIM
!-----------------------------------------------------------------------
  320 CONTINUE
  NZ=M
  DO 330 I=1,M
    Y(I) = 0.0E0
  330 CONTINUE
  return
end
