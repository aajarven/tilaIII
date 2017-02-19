subroutine DHKSEQ (X, M, H, IERR)
!
!! DHKSEQ is subsidiary to DBSKIN.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (HKSEQ-S, DHKSEQ-D)
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!   DHKSEQ is an adaptation of subroutine DPSIFN described in the
!   reference below.  DHKSEQ generates the sequence
!   H(K,X) = (-X)**(K+1)*(PSI(K,X) PSI(K,X+0.5))/GAMMA(K+1), for
!            K=0,...,M.
!
!***SEE ALSO  DBSKIN
!***REFERENCES  D. E. Amos, A portable Fortran subroutine for
!                 derivatives of the Psi function, Algorithm 610, ACM
!                 Transactions on Mathematical Software 9, 4 (1983),
!                 pp. 494-502.
!***ROUTINES CALLED  D1MACH, I1MACH
!***REVISION HISTORY  (YYMMDD)
!   820601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
!***END PROLOGUE  DHKSEQ
  INTEGER I, IERR, J, K, M, MX, NX
  INTEGER I1MACH
  DOUBLE PRECISION B, FK, FLN, FN, FNP, H, HRX, RLN, RXSQ, R1M5, S, &
   SLOPE, T, TK, TRM, TRMH, TRMR, TST, U, V, WDTOL, X, XDMY, XH, &
   XINC, XM, XMIN, YINT
  DOUBLE PRECISION D1MACH
  DIMENSION B(22), TRM(22), TRMR(25), TRMH(25), U(25), V(25), H(*)
  SAVE B
!-----------------------------------------------------------------------
!             SCALED BERNOULLI NUMBERS 2.0*B(2K)*(1-2**(-2K))
!-----------------------------------------------------------------------
  DATA B(1), B(2), B(3), B(4), B(5), B(6), B(7), B(8), B(9), B(10), &
   B(11), B(12), B(13), B(14), B(15), B(16), B(17), B(18), B(19), &
   B(20), B(21), B(22) /1.00000000000000000D+00, &
   -5.00000000000000000D-01,2.50000000000000000D-01, &
   -6.25000000000000000D-02,4.68750000000000000D-02, &
   -6.64062500000000000D-02,1.51367187500000000D-01, &
   -5.06103515625000000D-01,2.33319091796875000D+00, &
   -1.41840972900390625D+01,1.09941936492919922D+02, &
   -1.05824747562408447D+03,1.23842434241771698D+04, &
   -1.73160495905935764D+05,2.85103429084961116D+06, &
   -5.45964619322445132D+07,1.20316174668075304D+09, &
   -3.02326315271452307D+10,8.59229286072319606D+11, &
   -2.74233104097776039D+13,9.76664637943633248D+14, &
   -3.85931586838450360D+16/
!
!***FIRST EXECUTABLE STATEMENT  DHKSEQ
  IERR=0
  WDTOL = MAX(D1MACH(4),1.0D-18)
  FN = M - 1
  FNP = FN + 1.0D0
!-----------------------------------------------------------------------
!     COMPUTE XMIN
!-----------------------------------------------------------------------
  R1M5 = D1MACH(5)
  RLN = R1M5*I1MACH(14)
  RLN = MIN(RLN,18.06D0)
  FLN = MAX(RLN,3.0D0) - 3.0D0
  YINT = 3.50D0 + 0.40D0*FLN
  SLOPE = 0.21D0 + FLN*(0.0006038D0*FLN+0.008677D0)
  XM = YINT + SLOPE*FN
  MX = INT(XM) + 1
  XMIN = MX
!-----------------------------------------------------------------------
!     GENERATE H(M-1,XDMY)*XDMY**(M) BY THE ASYMPTOTIC EXPANSION
!-----------------------------------------------------------------------
  XDMY = X
  XINC = 0.0D0
  if (X >= XMIN) go to 10
  NX = INT(X)
  XINC = XMIN - NX
  XDMY = X + XINC
   10 CONTINUE
  RXSQ = 1.0D0/(XDMY*XDMY)
  HRX = 0.5D0/XDMY
  TST = 0.5D0*WDTOL
  T = FNP*HRX
!-----------------------------------------------------------------------
!     INITIALIZE COEFFICIENT ARRAY
!-----------------------------------------------------------------------
  S = T*B(3)
  if (ABS(S) < TST) go to 30
  TK = 2.0D0
  DO 20 K=4,22
    T = T*((TK+FN+1.0D0)/(TK+1.0D0))*((TK+FN)/(TK+2.0D0))*RXSQ
    TRM(K) = T*B(K)
    if (ABS(TRM(K)) < TST) go to 30
    S = S + TRM(K)
    TK = TK + 2.0D0
   20 CONTINUE
  go to 110
   30 CONTINUE
  H(M) = S + 0.5D0
  if (M == 1) go to 70
!-----------------------------------------------------------------------
!     GENERATE LOWER DERIVATIVES, I < M-1
!-----------------------------------------------------------------------
  DO 60 I=2,M
    FNP = FN
    FN = FN - 1.0D0
    S = FNP*HRX*B(3)
    if (ABS(S) < TST) go to 50
    FK = FNP + 3.0D0
    DO 40 K=4,22
      TRM(K) = TRM(K)*FNP/FK
      if (ABS(TRM(K)) < TST) go to 50
      S = S + TRM(K)
      FK = FK + 2.0D0
   40   CONTINUE
    go to 110
   50   CONTINUE
    MX = M - I + 1
    H(MX) = S + 0.5D0
   60 CONTINUE
   70 CONTINUE
  if (XINC == 0.0D0) RETURN
!-----------------------------------------------------------------------
!     RECUR BACKWARD FROM XDMY TO X
!-----------------------------------------------------------------------
  XH = X + 0.5D0
  S = 0.0D0
  NX = INT(XINC)
  DO 80 I=1,NX
    TRMR(I) = X/(X+NX-I)
    U(I) = TRMR(I)
    TRMH(I) = X/(XH+NX-I)
    V(I) = TRMH(I)
    S = S + U(I) - V(I)
   80 CONTINUE
  MX = NX + 1
  TRMR(MX) = X/XDMY
  U(MX) = TRMR(MX)
  H(1) = H(1)*TRMR(MX) + S
  if (M == 1) RETURN
  DO 100 J=2,M
    S = 0.0D0
    DO 90 I=1,NX
      TRMR(I) = TRMR(I)*U(I)
      TRMH(I) = TRMH(I)*V(I)
      S = S + TRMR(I) - TRMH(I)
   90   CONTINUE
    TRMR(MX) = TRMR(MX)*U(MX)
    H(J) = H(J)*TRMR(MX) + S
  100 CONTINUE
  return
  110 CONTINUE
  IERR=2
  return
end
