subroutine DPSIFN (X, N, KODE, M, ANS, NZ, IERR)
!
!! DPSIFN computes derivatives of the Psi function.
!
!***LIBRARY   SLATEC
!***CATEGORY  C7C
!***TYPE      DOUBLE PRECISION (PSIFN-S, DPSIFN-D)
!***KEYWORDS  DERIVATIVES OF THE GAMMA FUNCTION, POLYGAMMA FUNCTION,
!             PSI FUNCTION
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!         The following definitions are used in DPSIFN:
!
!      Definition 1
!         PSI(X) = d/dx (ln(GAMMA(X)), the first derivative of
!                  the log GAMMA function.
!      Definition 2
!                     K   K
!         PSI(K,X) = d /dx (PSI(X)), the K-th derivative of PSI(X).
!   ___________________________________________________________________
!      DPSIFN computes a sequence of SCALED derivatives of
!      the PSI function; i.e. for fixed X and M it computes
!      the M-member sequence
!
!                    ((-1)**(K+1)/GAMMA(K+1))*PSI(K,X)
!                       for K = N,...,N+M-1
!
!      where PSI(K,X) is as defined above.   For KODE=1, DPSIFN returns
!      the scaled derivatives as described.  KODE=2 is operative only
!      when K=0 and in that case DPSIFN returns -PSI(X) + LN(X).  That
!      is, the logarithmic behavior for large X is removed when KODE=2
!      and K=0.  When sums or differences of PSI functions are computed
!      the logarithmic terms can be combined analytically and computed
!      separately to help retain significant digits.
!
!         Note that call DPSIFN(X,0,1,1,ANS) results in
!                   ANS = -PSI(X)
!
!     Input      X is DOUBLE PRECISION
!           X      - Argument, X .gt. 0.0D0
!           N      - First member of the sequence, 0 .le. N .le. 100
!                    N=0 gives ANS(1) = -PSI(X)       for KODE=1
!                                       -PSI(X)+LN(X) for KODE=2
!           KODE   - Selection parameter
!                    KODE=1 returns scaled derivatives of the PSI
!                    function.
!                    KODE=2 returns scaled derivatives of the PSI
!                    function EXCEPT when N=0. In this case,
!                    ANS(1) = -PSI(X) + LN(X) is returned.
!           M      - Number of members of the sequence, M.ge.1
!
!    Output     ANS is DOUBLE PRECISION
!           ANS    - A vector of length at least M whose first M
!                    components contain the sequence of derivatives
!                    scaled according to KODE.
!           NZ     - Underflow flag
!                    NZ.eq.0, A normal return
!                    NZ.ne.0, Underflow, last NZ components of ANS are
!                             set to zero, ANS(M-K+1)=0.0, K=1,...,NZ
!           IERR   - Error flag
!                    IERR=0, A normal return, computation completed
!                    IERR=1, Input error,     no computation
!                    IERR=2, Overflow,        X too small or N+M-1 too
!                            large or both
!                    IERR=3, Error,           N too large. Dimensioned
!                            array TRMR(NMAX) is not large enough for N
!
!         The nominal computational accuracy is the maximum of unit
!         roundoff (=D1MACH(4)) and 1.0D-18 since critical constants
!         are given to only 18 digits.
!
!         PSIFN is the single precision version of DPSIFN.
!
! *Long Description:
!
!         The basic method of evaluation is the asymptotic expansion
!         for large X.ge.XMIN followed by backward recursion on a two
!         term recursion relation
!
!                  W(X+1) + X**(-N-1) = W(X).
!
!         This is supplemented by a series
!
!                  SUM( (X+K)**(-N-1) , K=0,1,2,... )
!
!         which converges rapidly for large N. Both XMIN and the
!         number of terms of the series are calculated from the unit
!         roundoff of the machine environment.
!
!***REFERENCES  Handbook of Mathematical Functions, National Bureau
!                 of Standards Applied Mathematics Series 55, edited
!                 by M. Abramowitz and I. A. Stegun, equations 6.3.5,
!                 6.3.18, 6.4.6, 6.4.9 and 6.4.10, pp.258-260, 1964.
!               D. E. Amos, A portable Fortran subroutine for
!                 derivatives of the Psi function, Algorithm 610, ACM
!                 Transactions on Mathematical Software 9, 4 (1983),
!                 pp. 494-502.
!***ROUTINES CALLED  D1MACH, I1MACH
!***REVISION HISTORY  (YYMMDD)
!   820601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DPSIFN
  INTEGER I, IERR, J, K, KODE, M, MM, MX, N, NMAX, NN, NP, NX, NZ, &
    FN
  INTEGER I1MACH
  DOUBLE PRECISION ANS, ARG, B, DEN, ELIM, EPS, FLN, &
   FX, RLN, RXSQ, R1M4, R1M5, S, SLOPE, T, TA, TK, TOL, TOLS, TRM, &
   TRMR, TSS, TST, TT, T1, T2, WDTOL, X, XDMLN, XDMY, XINC, XLN, &
   XM, XMIN, XQ, YINT
  DOUBLE PRECISION D1MACH
  DIMENSION B(22), TRM(22), TRMR(100), ANS(*)
  SAVE NMAX, B
  DATA NMAX /100/
!-----------------------------------------------------------------------
!             BERNOULLI NUMBERS
!-----------------------------------------------------------------------
  DATA B(1), B(2), B(3), B(4), B(5), B(6), B(7), B(8), B(9), B(10), &
   B(11), B(12), B(13), B(14), B(15), B(16), B(17), B(18), B(19), &
   B(20), B(21), B(22) /1.00000000000000000D+00, &
   -5.00000000000000000D-01,1.66666666666666667D-01, &
   -3.33333333333333333D-02,2.38095238095238095D-02, &
   -3.33333333333333333D-02,7.57575757575757576D-02, &
   -2.53113553113553114D-01,1.16666666666666667D+00, &
   -7.09215686274509804D+00,5.49711779448621554D+01, &
   -5.29124242424242424D+02,6.19212318840579710D+03, &
   -8.65802531135531136D+04,1.42551716666666667D+06, &
   -2.72982310678160920D+07,6.01580873900642368D+08, &
   -1.51163157670921569D+10,4.29614643061166667D+11, &
   -1.37116552050883328D+13,4.88332318973593167D+14, &
   -1.92965793419400681D+16/
!
!***FIRST EXECUTABLE STATEMENT  DPSIFN
  IERR = 0
  NZ=0
  if (X <= 0.0D0) IERR=1
  if (N < 0) IERR=1
  if (KODE < 1 .OR. KODE > 2) IERR=1
  if (M < 1) IERR=1
  if (IERR /= 0) RETURN
  MM=M
  NX = MIN(-I1MACH(15),I1MACH(16))
  R1M5 = D1MACH(5)
  R1M4 = D1MACH(4)*0.5D0
  WDTOL = MAX(R1M4,0.5D-18)
!-----------------------------------------------------------------------
!     ELIM = APPROXIMATE EXPONENTIAL OVER AND UNDERFLOW LIMIT
!-----------------------------------------------------------------------
  ELIM = 2.302D0*(NX*R1M5-3.0D0)
  XLN = LOG(X)
   41 CONTINUE
  NN = N + MM - 1
  FN = NN
  T = (FN+1)*XLN
!-----------------------------------------------------------------------
!     OVERFLOW AND UNDERFLOW TEST FOR SMALL AND LARGE X
!-----------------------------------------------------------------------
  if (ABS(T) > ELIM) go to 290
  if (X < WDTOL) go to 260
!-----------------------------------------------------------------------
!     COMPUTE XMIN AND THE NUMBER OF TERMS OF THE SERIES, FLN+1
!-----------------------------------------------------------------------
  RLN = R1M5*I1MACH(14)
  RLN = MIN(RLN,18.06D0)
  FLN = MAX(RLN,3.0D0) - 3.0D0
  YINT = 3.50D0 + 0.40D0*FLN
  SLOPE = 0.21D0 + FLN*(0.0006038D0*FLN+0.008677D0)
  XM = YINT + SLOPE*FN
  MX = INT(XM) + 1
  XMIN = MX
  if (N == 0) go to 50
  XM = -2.302D0*RLN - MIN(0.0D0,XLN)
  ARG = XM/N
  ARG = MIN(0.0D0,ARG)
  EPS = EXP(ARG)
  XM = 1.0D0 - EPS
  if (ABS(ARG) < 1.0D-3) XM = -ARG
  FLN = X*XM/EPS
  XM = XMIN - X
  if (XM > 7.0D0 .AND. FLN < 15.0D0) go to 200
   50 CONTINUE
  XDMY = X
  XDMLN = XLN
  XINC = 0.0D0
  if (X >= XMIN) go to 60
  NX = INT(X)
  XINC = XMIN - NX
  XDMY = X + XINC
  XDMLN = LOG(XDMY)
   60 CONTINUE
!-----------------------------------------------------------------------
!     GENERATE W(N+MM-1,X) BY THE ASYMPTOTIC EXPANSION
!-----------------------------------------------------------------------
  T = FN*XDMLN
  T1 = XDMLN + XDMLN
  T2 = T + XDMLN
  TK = MAX(ABS(T),ABS(T1),ABS(T2))
  if (TK > ELIM) go to 380
  TSS = EXP(-T)
  TT = 0.5D0/XDMY
  T1 = TT
  TST = WDTOL*TT
  if (NN /= 0) T1 = TT + 1.0D0/FN
  RXSQ = 1.0D0/(XDMY*XDMY)
  TA = 0.5D0*RXSQ
  T = (FN+1)*TA
  S = T*B(3)
  if (ABS(S) < TST) go to 80
  TK = 2.0D0
  DO 70 K=4,22
    T = T*((TK+FN+1)/(TK+1.0D0))*((TK+FN)/(TK+2.0D0))*RXSQ
    TRM(K) = T*B(K)
    if (ABS(TRM(K)) < TST) go to 80
    S = S + TRM(K)
    TK = TK + 2.0D0
   70 CONTINUE
   80 CONTINUE
  S = (S+T1)*TSS
  if (XINC == 0.0D0) go to 100
!-----------------------------------------------------------------------
!     BACKWARD RECUR FROM XDMY TO X
!-----------------------------------------------------------------------
  NX = INT(XINC)
  NP = NN + 1
  if (NX > NMAX) go to 390
  if (NN == 0) go to 160
  XM = XINC - 1.0D0
  FX = X + XM
!-----------------------------------------------------------------------
!     THIS LOOP SHOULD NOT BE CHANGED. FX IS ACCURATE WHEN X IS SMALL
!-----------------------------------------------------------------------
  DO 90 I=1,NX
    TRMR(I) = FX**(-NP)
    S = S + TRMR(I)
    XM = XM - 1.0D0
    FX = X + XM
   90 CONTINUE
  100 CONTINUE
  ANS(MM) = S
  if (FN == 0) go to 180
!-----------------------------------------------------------------------
!     GENERATE LOWER DERIVATIVES, J < N+MM-1
!-----------------------------------------------------------------------
  if (MM == 1) RETURN
  DO 150 J=2,MM
    FN = FN - 1
    TSS = TSS*XDMY
    T1 = TT
    if (FN /= 0) T1 = TT + 1.0D0/FN
    T = (FN+1)*TA
    S = T*B(3)
    if (ABS(S) < TST) go to 120
    TK = 4 + FN
    DO 110 K=4,22
      TRM(K) = TRM(K)*(FN+1)/TK
      if (ABS(TRM(K)) < TST) go to 120
      S = S + TRM(K)
      TK = TK + 2.0D0
  110   CONTINUE
  120   CONTINUE
    S = (S+T1)*TSS
    if (XINC == 0.0D0) go to 140
    if (FN == 0) go to 160
    XM = XINC - 1.0D0
    FX = X + XM
    DO 130 I=1,NX
      TRMR(I) = TRMR(I)*FX
      S = S + TRMR(I)
      XM = XM - 1.0D0
      FX = X + XM
  130   CONTINUE
  140   CONTINUE
    MX = MM - J + 1
    ANS(MX) = S
    if (FN == 0) go to 180
  150 CONTINUE
  return
!-----------------------------------------------------------------------
!     RECURSION FOR N = 0
!-----------------------------------------------------------------------
  160 CONTINUE
  DO 170 I=1,NX
    S = S + 1.0D0/(X+NX-I)
  170 CONTINUE
  180 CONTINUE
  if (KODE == 2) go to 190
  ANS(1) = S - XDMLN
  return
  190 CONTINUE
  if (XDMY == X) RETURN
  XQ = XDMY/X
  ANS(1) = S - LOG(XQ)
  return
!-----------------------------------------------------------------------
!     COMPUTE BY SERIES (X+K)**(-(N+1)) , K=0,1,2,...
!-----------------------------------------------------------------------
  200 CONTINUE
  NN = INT(FLN) + 1
  NP = N + 1
  T1 = (N+1)*XLN
  T = EXP(-T1)
  S = T
  DEN = X
  DO 210 I=1,NN
    DEN = DEN + 1.0D0
    TRM(I) = DEN**(-NP)
    S = S + TRM(I)
  210 CONTINUE
  ANS(1) = S
  if (N /= 0) go to 220
  if (KODE == 2) ANS(1) = S + XLN
  220 CONTINUE
  if (MM == 1) RETURN
!-----------------------------------------------------------------------
!     GENERATE HIGHER DERIVATIVES, J > N
!-----------------------------------------------------------------------
  TOL = WDTOL/5.0D0
  DO 250 J=2,MM
    T = T/X
    S = T
    TOLS = T*TOL
    DEN = X
    DO 230 I=1,NN
      DEN = DEN + 1.0D0
      TRM(I) = TRM(I)/DEN
      S = S + TRM(I)
      if (TRM(I) < TOLS) go to 240
  230   CONTINUE
  240   CONTINUE
    ANS(J) = S
  250 CONTINUE
  return
!-----------------------------------------------------------------------
!     SMALL X < UNIT ROUND OFF
!-----------------------------------------------------------------------
  260 CONTINUE
  ANS(1) = X**(-N-1)
  if (MM == 1) go to 280
  K = 1
  DO 270 I=2,MM
    ANS(K+1) = ANS(K)/X
    K = K + 1
  270 CONTINUE
  280 CONTINUE
  if (N /= 0) RETURN
  if (KODE == 2) ANS(1) = ANS(1) + XLN
  return
  290 CONTINUE
  if (T > 0.0D0) go to 380
  NZ=0
  IERR=2
  return
  380 CONTINUE
  NZ=NZ+1
  ANS(MM)=0.0D0
  MM=MM-1
  if (MM == 0) RETURN
  go to 41
  390 CONTINUE
  NZ=0
  IERR=3
  return
end
