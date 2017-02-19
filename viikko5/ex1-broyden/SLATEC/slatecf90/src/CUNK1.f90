subroutine CUNK1 (Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
!
!! CUNK1 is subsidiary to CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNK1-A, ZUNK1-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
!     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
!     UNIFORM ASYMPTOTIC EXPANSION.
!     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
!     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
!
!***SEE ALSO  CBESK
!***ROUTINES CALLED  CS1S2, CUCHK, CUNIK, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUNK1
  COMPLEX CFN, CK, CONE, CRSC, CS, CSCL, CSGN, CSPN, CSR, CSS, &
   CWRK, CY, CZERO, C1, C2, PHI,  RZ, SUM,  S1, S2, Y, Z, &
   ZETA1,  ZETA2,  ZR, PHID, ZETA1D, ZETA2D, SUMD
  REAL ALIM, ANG, APHI, ASC, ASCLE, BRY, CPN, C2I, C2M, C2R, ELIM, &
   FMR, FN, FNF, FNU, PI, RS1, SGN, SPN, TOL, X, R1MACH
  INTEGER I, IB, IFLAG, IFN, IL, INIT, INU, IUF, K, KDFLG, KFLAG, &
   KK, KODE, MR, N, NW, NZ, J, IPARD, INITD, IC, M
  DIMENSION BRY(3), INIT(2), Y(N), SUM(2), PHI(2), ZETA1(2), &
   ZETA2(2), CY(2), CWRK(16,3), CSS(3), CSR(3)
  DATA CZERO, CONE / (0.0E0,0.0E0) , (1.0E0,0.0E0) /
  DATA PI / 3.14159265358979324E0 /
!***FIRST EXECUTABLE STATEMENT  CUNK1
  KDFLG = 1
  NZ = 0
!-----------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
!     THE UNDERFLOW LIMIT
!-----------------------------------------------------------------------
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CRSC = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CRSC
  CSR(1) = CRSC
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = 1.0E+3*R1MACH(1)/TOL
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = R1MACH(2)
  X = REAL(Z)
  ZR = Z
  if (X < 0.0E0) ZR = -Z
  J=2
  DO 70 I=1,N
!-----------------------------------------------------------------------
!     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
!-----------------------------------------------------------------------
    J = 3 - J
    FN = FNU + (I-1)
    INIT(J) = 0
    call CUNIK(ZR, FN, 2, 0, TOL, INIT(J), PHI(J), ZETA1(J), &
     ZETA2(J), SUM(J), CWRK(1,J))
    if (KODE == 1) go to 20
    CFN = CMPLX(FN,0.0E0)
    S1 = ZETA1(J) - CFN*(CFN/(ZR+ZETA2(J)))
    go to 30
   20   CONTINUE
    S1 = ZETA1(J) - ZETA2(J)
   30   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = REAL(S1)
    if (ABS(RS1) > ELIM) go to 60
    if (KDFLG == 1) KFLAG = 2
    if (ABS(RS1) < ALIM) go to 40
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    APHI = ABS(PHI(J))
    RS1 = RS1 + ALOG(APHI)
    if (ABS(RS1) > ELIM) go to 60
    if (KDFLG == 1) KFLAG = 1
    if (RS1 < 0.0E0) go to 40
    if (KDFLG == 1) KFLAG = 3
   40   CONTINUE
!-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!-----------------------------------------------------------------------
    S2 = PHI(J)*SUM(J)
    C2R = REAL(S1)
    C2I = AIMAG(S1)
    C2M = EXP(C2R)*REAL(CSS(KFLAG))
    S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
    S2 = S2*S1
    if (KFLAG /= 1) go to 50
    call CUCHK(S2, NW, BRY(1), TOL)
    if (NW /= 0) go to 60
   50   CONTINUE
    CY(KDFLG) = S2
    Y(I) = S2*CSR(KFLAG)
    if (KDFLG == 2) go to 75
    KDFLG = 2
    go to 70
   60   CONTINUE
    if (RS1 > 0.0E0) go to 290
!-----------------------------------------------------------------------
!     FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
    if (X < 0.0E0) go to 290
    KDFLG = 1
    Y(I) = CZERO
    NZ=NZ+1
    if (I == 1) go to 70
    if (Y(I-1) == CZERO) go to 70
    Y(I-1) = CZERO
    NZ=NZ+1
   70 CONTINUE
  I=N
   75 CONTINUE
  RZ = CMPLX(2.0E0,0.0E0)/ZR
  CK = CMPLX(FN,0.0E0)*RZ
  IB = I+1
  if (N < IB) go to 160
!-----------------------------------------------------------------------
!     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
!     ON UNDERFLOW
!-----------------------------------------------------------------------
  FN = FNU+(N-1)
  IPARD = 1
  if (MR /= 0) IPARD = 0
  INITD = 0
  call CUNIK(ZR,FN,2,IPARD,TOL,INITD,PHID,ZETA1D,ZETA2D,SUMD, &
  CWRK(1,3))
  if (KODE == 1) go to 80
  CFN=CMPLX(FN,0.0E0)
  S1=ZETA1D-CFN*(CFN/(ZR+ZETA2D))
  go to 90
   80 CONTINUE
  S1=ZETA1D-ZETA2D
   90 CONTINUE
  RS1=REAL(S1)
  if (ABS(RS1) > ELIM) go to 95
  if (ABS(RS1) < ALIM) go to 100
!-----------------------------------------------------------------------
!     REFINE ESTIMATE AND TEST
!-----------------------------------------------------------------------
  APHI=ABS(PHID)
  RS1=RS1+ALOG(APHI)
  if (ABS(RS1) < ELIM) go to 100
   95 CONTINUE
  if (RS1 > 0.0E0) go to 290
!-----------------------------------------------------------------------
!     FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
  if (X < 0.0E0) go to 290
  NZ=N
  DO 96 I=1,N
    Y(I) = CZERO
   96 CONTINUE
  return
  100 CONTINUE
!-----------------------------------------------------------------------
!     RECUR FORWARD FOR REMAINDER OF THE SEQUENCE
!-----------------------------------------------------------------------
  S1 = CY(1)
  S2 = CY(2)
  C1 = CSR(KFLAG)
  ASCLE = BRY(KFLAG)
  DO 120 I=IB,N
    C2 = S2
    S2 = CK*S2 + S1
    S1 = C2
    CK = CK + RZ
    C2 = S2*C1
    Y(I) = C2
    if (KFLAG >= 3) go to 120
    C2R = REAL(C2)
    C2I = AIMAG(C2)
    C2R = ABS(C2R)
    C2I = ABS(C2I)
    C2M = MAX(C2R,C2I)
    if (C2M <= ASCLE) go to 120
    KFLAG = KFLAG + 1
    ASCLE = BRY(KFLAG)
    S1 = S1*C1
    S2 = C2
    S1 = S1*CSS(KFLAG)
    S2 = S2*CSS(KFLAG)
    C1 = CSR(KFLAG)
  120 CONTINUE
  160 CONTINUE
  if (MR == 0) RETURN
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION FOR RE(Z) < 0.0E0
!-----------------------------------------------------------------------
  NZ = 0
  FMR = MR
  SGN = -SIGN(PI,FMR)
!-----------------------------------------------------------------------
!     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
!-----------------------------------------------------------------------
  CSGN = CMPLX(0.0E0,SGN)
  INU = FNU
  FNF = FNU - INU
  IFN = INU + N - 1
  ANG = FNF*SGN
  CPN = COS(ANG)
  SPN = SIN(ANG)
  CSPN = CMPLX(CPN,SPN)
  if (MOD(IFN,2) == 1) CSPN = -CSPN
  ASC = BRY(1)
  KK = N
  IUF = 0
  KDFLG = 1
  IB = IB-1
  IC = IB-1
  DO 260 K=1,N
    FN = FNU + (KK-1)
!-----------------------------------------------------------------------
!     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
!     FUNCTION ABOVE
!-----------------------------------------------------------------------
    M=3
    if (N > 2) go to 175
  170   CONTINUE
    INITD = INIT(J)
    PHID = PHI(J)
    ZETA1D = ZETA1(J)
    ZETA2D = ZETA2(J)
    SUMD = SUM(J)
    M = J
    J = 3 - J
    go to 180
  175   CONTINUE
    if ((KK == N).AND.(IB < N)) go to 180
    if ((KK == IB).OR.(KK == IC)) go to 170
    INITD = 0
  180   CONTINUE
    call CUNIK(ZR, FN, 1, 0, TOL, INITD, PHID, ZETA1D, &
     ZETA2D, SUMD, CWRK(1,M))
    if (KODE == 1) go to 190
    CFN = CMPLX(FN,0.0E0)
    S1 = -ZETA1D + CFN*(CFN/(ZR+ZETA2D))
    go to 200
  190   CONTINUE
    S1 = -ZETA1D + ZETA2D
  200   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = REAL(S1)
    if (ABS(RS1) > ELIM) go to 250
    if (KDFLG == 1) IFLAG = 2
    if (ABS(RS1) < ALIM) go to 210
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    APHI = ABS(PHID)
    RS1 = RS1 + ALOG(APHI)
    if (ABS(RS1) > ELIM) go to 250
    if (KDFLG == 1) IFLAG = 1
    if (RS1 < 0.0E0) go to 210
    if (KDFLG == 1) IFLAG = 3
  210   CONTINUE
    S2 = CSGN*PHID*SUMD
    C2R = REAL(S1)
    C2I = AIMAG(S1)
    C2M = EXP(C2R)*REAL(CSS(IFLAG))
    S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
    S2 = S2*S1
    if (IFLAG /= 1) go to 220
    call CUCHK(S2, NW, BRY(1), TOL)
    if (NW /= 0) S2 = CMPLX(0.0E0,0.0E0)
  220   CONTINUE
    CY(KDFLG) = S2
    C2 = S2
    S2 = S2*CSR(IFLAG)
!-----------------------------------------------------------------------
!     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
!-----------------------------------------------------------------------
    S1 = Y(KK)
    if (KODE == 1) go to 240
    call CS1S2(ZR, S1, S2, NW, ASC, ALIM, IUF)
    NZ = NZ + NW
  240   CONTINUE
    Y(KK) = S1*CSPN + S2
    KK = KK - 1
    CSPN = -CSPN
    if (C2 /= CZERO) go to 245
    KDFLG = 1
    go to 260
  245   CONTINUE
    if (KDFLG == 2) go to 265
    KDFLG = 2
    go to 260
  250   CONTINUE
    if (RS1 > 0.0E0) go to 290
    S2 = CZERO
    go to 220
  260 CONTINUE
  K = N
  265 CONTINUE
  IL = N - K
  if (IL == 0) RETURN
!-----------------------------------------------------------------------
!     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
!     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
!     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
!-----------------------------------------------------------------------
  S1 = CY(1)
  S2 = CY(2)
  CS = CSR(IFLAG)
  ASCLE = BRY(IFLAG)
  FN = (INU+IL)
  DO 280 I=1,IL
    C2 = S2
    S2 = S1 + CMPLX(FN+FNF,0.0E0)*RZ*S2
    S1 = C2
    FN = FN - 1.0E0
    C2 = S2*CS
    CK = C2
    C1 = Y(KK)
    if (KODE == 1) go to 270
    call CS1S2(ZR, C1, C2, NW, ASC, ALIM, IUF)
    NZ = NZ + NW
  270   CONTINUE
    Y(KK) = C1*CSPN + C2
    KK = KK - 1
    CSPN = -CSPN
    if (IFLAG >= 3) go to 280
    C2R = REAL(CK)
    C2I = AIMAG(CK)
    C2R = ABS(C2R)
    C2I = ABS(C2I)
    C2M = MAX(C2R,C2I)
    if (C2M <= ASCLE) go to 280
    IFLAG = IFLAG + 1
    ASCLE = BRY(IFLAG)
    S1 = S1*CS
    S2 = CK
    S1 = S1*CSS(IFLAG)
    S2 = S2*CSS(IFLAG)
    CS = CSR(IFLAG)
  280 CONTINUE
  return
  290 CONTINUE
  NZ = -1
  return
end
