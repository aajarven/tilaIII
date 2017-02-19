subroutine ZUNK1 (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, &
     ALIM)
!
!! ZUNK1 is subsidiary to ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNK1-A, ZUNK1-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
!     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
!     UNIFORM ASYMPTOTIC EXPANSION.
!     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
!     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
!
!***SEE ALSO  ZBESK
!***ROUTINES CALLED  D1MACH, ZABS, ZS1S2, ZUCHK, ZUNIK
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ZUNK1
!     COMPLEX CFN,CK,CONE,CRSC,CS,CSCL,CSGN,CSPN,CSR,CSS,CWRK,CY,CZERO,
!    *C1,C2,PHI,PHID,RZ,SUM,SUMD,S1,S2,Y,Z,ZETA1,ZETA1D,ZETA2,ZETA2D,ZR
  DOUBLE PRECISION ALIM, ANG, APHI, ASC, ASCLE, BRY, CKI, CKR, &
   CONER, CRSC, CSCL, CSGNI, CSPNI, CSPNR, CSR, CSRR, CSSR, &
   CWRKI, CWRKR, CYI, CYR, C1I, C1R, C2I, C2M, C2R, ELIM, FMR, FN, &
   FNF, FNU, PHIDI, PHIDR, PHII, PHIR, PI, RAST, RAZR, RS1, RZI, &
   RZR, SGN, STI, STR, SUMDI, SUMDR, SUMI, SUMR, S1I, S1R, S2I, &
   S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, &
   ZET1DI, ZET1DR, ZET2DI, ZET2DR, ZI, ZR, ZRI, ZRR, D1MACH, ZABS
  INTEGER I, IB, IFLAG, IFN, IL, INIT, INU, IUF, K, KDFLG, KFLAG, &
   KK, KODE, MR, N, NW, NZ, INITD, IC, IPARD, J, M
  DIMENSION BRY(3), INIT(2), YR(N), YI(N), SUMR(2), SUMI(2), &
   ZETA1R(2), ZETA1I(2), ZETA2R(2), ZETA2I(2), CYR(2), CYI(2), &
   CWRKR(16,3), CWRKI(16,3), CSSR(3), CSRR(3), PHIR(2), PHII(2)
  EXTERNAL ZABS
  DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
  DATA PI / 3.14159265358979324D0 /
!***FIRST EXECUTABLE STATEMENT  ZUNK1
  KDFLG = 1
  NZ = 0
!-----------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
!     THE UNDERFLOW LIMIT
!-----------------------------------------------------------------------
  CSCL = 1.0D0/TOL
  CRSC = TOL
  CSSR(1) = CSCL
  CSSR(2) = CONER
  CSSR(3) = CRSC
  CSRR(1) = CRSC
  CSRR(2) = CONER
  CSRR(3) = CSCL
  BRY(1) = 1.0D+3*D1MACH(1)/TOL
  BRY(2) = 1.0D0/BRY(1)
  BRY(3) = D1MACH(2)
  ZRR = ZR
  ZRI = ZI
  if (ZR >= 0.0D0) go to 10
  ZRR = -ZR
  ZRI = -ZI
   10 CONTINUE
  J = 2
  DO 70 I=1,N
!-----------------------------------------------------------------------
!     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
!-----------------------------------------------------------------------
    J = 3 - J
    FN = FNU + (I-1)
    INIT(J) = 0
    call ZUNIK(ZRR, ZRI, FN, 2, 0, TOL, INIT(J), PHIR(J), PHII(J), &
     ZETA1R(J), ZETA1I(J), ZETA2R(J), ZETA2I(J), SUMR(J), SUMI(J), &
     CWRKR(1,J), CWRKI(1,J))
    if (KODE == 1) go to 20
    STR = ZRR + ZETA2R(J)
    STI = ZRI + ZETA2I(J)
    RAST = FN/ZABS(STR,STI)
    STR = STR*RAST*RAST
    STI = -STI*RAST*RAST
    S1R = ZETA1R(J) - STR
    S1I = ZETA1I(J) - STI
    go to 30
   20   CONTINUE
    S1R = ZETA1R(J) - ZETA2R(J)
    S1I = ZETA1I(J) - ZETA2I(J)
   30   CONTINUE
    RS1 = S1R
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    if (ABS(RS1) > ELIM) go to 60
    if (KDFLG == 1) KFLAG = 2
    if (ABS(RS1) < ALIM) go to 40
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    APHI = ZABS(PHIR(J),PHII(J))
    RS1 = RS1 + LOG(APHI)
    if (ABS(RS1) > ELIM) go to 60
    if (KDFLG == 1) KFLAG = 1
    if (RS1 < 0.0D0) go to 40
    if (KDFLG == 1) KFLAG = 3
   40   CONTINUE
!-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!-----------------------------------------------------------------------
    S2R = PHIR(J)*SUMR(J) - PHII(J)*SUMI(J)
    S2I = PHIR(J)*SUMI(J) + PHII(J)*SUMR(J)
    STR = EXP(S1R)*CSSR(KFLAG)
    S1R = STR*COS(S1I)
    S1I = STR*SIN(S1I)
    STR = S2R*S1R - S2I*S1I
    S2I = S1R*S2I + S2R*S1I
    S2R = STR
    if (KFLAG /= 1) go to 50
    call ZUCHK(S2R, S2I, NW, BRY(1), TOL)
    if (NW /= 0) go to 60
   50   CONTINUE
    CYR(KDFLG) = S2R
    CYI(KDFLG) = S2I
    YR(I) = S2R*CSRR(KFLAG)
    YI(I) = S2I*CSRR(KFLAG)
    if (KDFLG == 2) go to 75
    KDFLG = 2
    go to 70
   60   CONTINUE
    if (RS1 > 0.0D0) go to 300
!-----------------------------------------------------------------------
!     FOR ZR < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
    if (ZR < 0.0D0) go to 300
    KDFLG = 1
    YR(I)=ZEROR
    YI(I)=ZEROI
    NZ=NZ+1
    if (I == 1) go to 70
    if ((YR(I-1) == ZEROR).AND.(YI(I-1) == ZEROI)) go to 70
    YR(I-1)=ZEROR
    YI(I-1)=ZEROI
    NZ=NZ+1
   70 CONTINUE
  I = N
   75 CONTINUE
  RAZR = 1.0D0/ZABS(ZRR,ZRI)
  STR = ZRR*RAZR
  STI = -ZRI*RAZR
  RZR = (STR+STR)*RAZR
  RZI = (STI+STI)*RAZR
  CKR = FN*RZR
  CKI = FN*RZI
  IB = I + 1
  if (N < IB) go to 160
!-----------------------------------------------------------------------
!     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
!     ON UNDERFLOW.
!-----------------------------------------------------------------------
  FN = FNU + (N-1)
  IPARD = 1
  if (MR /= 0) IPARD = 0
  INITD = 0
  call ZUNIK(ZRR, ZRI, FN, 2, IPARD, TOL, INITD, PHIDR, PHIDI, &
   ZET1DR, ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI, CWRKR(1,3), &
   CWRKI(1,3))
  if (KODE == 1) go to 80
  STR = ZRR + ZET2DR
  STI = ZRI + ZET2DI
  RAST = FN/ZABS(STR,STI)
  STR = STR*RAST*RAST
  STI = -STI*RAST*RAST
  S1R = ZET1DR - STR
  S1I = ZET1DI - STI
  go to 90
   80 CONTINUE
  S1R = ZET1DR - ZET2DR
  S1I = ZET1DI - ZET2DI
   90 CONTINUE
  RS1 = S1R
  if (ABS(RS1) > ELIM) go to 95
  if (ABS(RS1) < ALIM) go to 100
!-----------------------------------------------------------------------
!     REFINE ESTIMATE AND TEST
!-----------------------------------------------------------------------
  APHI = ZABS(PHIDR,PHIDI)
  RS1 = RS1+LOG(APHI)
  if (ABS(RS1) < ELIM) go to 100
   95 CONTINUE
  if (ABS(RS1) > 0.0D0) go to 300
!-----------------------------------------------------------------------
!     FOR ZR < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
  if (ZR < 0.0D0) go to 300
  NZ = N
  DO 96 I=1,N
    YR(I) = ZEROR
    YI(I) = ZEROI
   96 CONTINUE
  return
!-----------------------------------------------------------------------
!     FORWARD RECUR FOR REMAINDER OF THE SEQUENCE
!-----------------------------------------------------------------------
  100 CONTINUE
  S1R = CYR(1)
  S1I = CYI(1)
  S2R = CYR(2)
  S2I = CYI(2)
  C1R = CSRR(KFLAG)
  ASCLE = BRY(KFLAG)
  DO 120 I=IB,N
    C2R = S2R
    C2I = S2I
    S2R = CKR*C2R - CKI*C2I + S1R
    S2I = CKR*C2I + CKI*C2R + S1I
    S1R = C2R
    S1I = C2I
    CKR = CKR + RZR
    CKI = CKI + RZI
    C2R = S2R*C1R
    C2I = S2I*C1R
    YR(I) = C2R
    YI(I) = C2I
    if (KFLAG >= 3) go to 120
    STR = ABS(C2R)
    STI = ABS(C2I)
    C2M = MAX(STR,STI)
    if (C2M <= ASCLE) go to 120
    KFLAG = KFLAG + 1
    ASCLE = BRY(KFLAG)
    S1R = S1R*C1R
    S1I = S1I*C1R
    S2R = C2R
    S2I = C2I
    S1R = S1R*CSSR(KFLAG)
    S1I = S1I*CSSR(KFLAG)
    S2R = S2R*CSSR(KFLAG)
    S2I = S2I*CSSR(KFLAG)
    C1R = CSRR(KFLAG)
  120 CONTINUE
  160 CONTINUE
  if (MR == 0) RETURN
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION FOR RE(Z) < 0.0D0
!-----------------------------------------------------------------------
  NZ = 0
  FMR = MR
  SGN = -DSIGN(PI,FMR)
!-----------------------------------------------------------------------
!     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
!-----------------------------------------------------------------------
  CSGNI = SGN
  INU = FNU
  FNF = FNU - INU
  IFN = INU + N - 1
  ANG = FNF*SGN
  CSPNR = COS(ANG)
  CSPNI = SIN(ANG)
  if (MOD(IFN,2) == 0) go to 170
  CSPNR = -CSPNR
  CSPNI = -CSPNI
  170 CONTINUE
  ASC = BRY(1)
  IUF = 0
  KK = N
  KDFLG = 1
  IB = IB - 1
  IC = IB - 1
  DO 270 K=1,N
    FN = FNU + (KK-1)
!-----------------------------------------------------------------------
!     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
!     FUNCTION ABOVE
!-----------------------------------------------------------------------
    M=3
    if (N > 2) go to 175
  172   CONTINUE
    INITD = INIT(J)
    PHIDR = PHIR(J)
    PHIDI = PHII(J)
    ZET1DR = ZETA1R(J)
    ZET1DI = ZETA1I(J)
    ZET2DR = ZETA2R(J)
    ZET2DI = ZETA2I(J)
    SUMDR = SUMR(J)
    SUMDI = SUMI(J)
    M = J
    J = 3 - J
    go to 180
  175   CONTINUE
    if ((KK == N).AND.(IB < N)) go to 180
    if ((KK == IB).OR.(KK == IC)) go to 172
    INITD = 0
  180   CONTINUE
    call ZUNIK(ZRR, ZRI, FN, 1, 0, TOL, INITD, PHIDR, PHIDI, &
     ZET1DR, ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI, &
     CWRKR(1,M), CWRKI(1,M))
    if (KODE == 1) go to 200
    STR = ZRR + ZET2DR
    STI = ZRI + ZET2DI
    RAST = FN/ZABS(STR,STI)
    STR = STR*RAST*RAST
    STI = -STI*RAST*RAST
    S1R = -ZET1DR + STR
    S1I = -ZET1DI + STI
    go to 210
  200   CONTINUE
    S1R = -ZET1DR + ZET2DR
    S1I = -ZET1DI + ZET2DI
  210   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = S1R
    if (ABS(RS1) > ELIM) go to 260
    if (KDFLG == 1) IFLAG = 2
    if (ABS(RS1) < ALIM) go to 220
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    APHI = ZABS(PHIDR,PHIDI)
    RS1 = RS1 + LOG(APHI)
    if (ABS(RS1) > ELIM) go to 260
    if (KDFLG == 1) IFLAG = 1
    if (RS1 < 0.0D0) go to 220
    if (KDFLG == 1) IFLAG = 3
  220   CONTINUE
    STR = PHIDR*SUMDR - PHIDI*SUMDI
    STI = PHIDR*SUMDI + PHIDI*SUMDR
    S2R = -CSGNI*STI
    S2I = CSGNI*STR
    STR = EXP(S1R)*CSSR(IFLAG)
    S1R = STR*COS(S1I)
    S1I = STR*SIN(S1I)
    STR = S2R*S1R - S2I*S1I
    S2I = S2R*S1I + S2I*S1R
    S2R = STR
    if (IFLAG /= 1) go to 230
    call ZUCHK(S2R, S2I, NW, BRY(1), TOL)
    if (NW == 0) go to 230
    S2R = ZEROR
    S2I = ZEROI
  230   CONTINUE
    CYR(KDFLG) = S2R
    CYI(KDFLG) = S2I
    C2R = S2R
    C2I = S2I
    S2R = S2R*CSRR(IFLAG)
    S2I = S2I*CSRR(IFLAG)
!-----------------------------------------------------------------------
!     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
!-----------------------------------------------------------------------
    S1R = YR(KK)
    S1I = YI(KK)
    if (KODE == 1) go to 250
    call ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF)
    NZ = NZ + NW
  250   CONTINUE
    YR(KK) = S1R*CSPNR - S1I*CSPNI + S2R
    YI(KK) = CSPNR*S1I + CSPNI*S1R + S2I
    KK = KK - 1
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    if (C2R /= 0.0D0 .OR. C2I /= 0.0D0) go to 255
    KDFLG = 1
    go to 270
  255   CONTINUE
    if (KDFLG == 2) go to 275
    KDFLG = 2
    go to 270
  260   CONTINUE
    if (RS1 > 0.0D0) go to 300
    S2R = ZEROR
    S2I = ZEROI
    go to 230
  270 CONTINUE
  K = N
  275 CONTINUE
  IL = N - K
  if (IL == 0) RETURN
!-----------------------------------------------------------------------
!     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
!     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
!     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
!-----------------------------------------------------------------------
  S1R = CYR(1)
  S1I = CYI(1)
  S2R = CYR(2)
  S2I = CYI(2)
  CSR = CSRR(IFLAG)
  ASCLE = BRY(IFLAG)
  FN = INU+IL
  DO 290 I=1,IL
    C2R = S2R
    C2I = S2I
    S2R = S1R + (FN+FNF)*(RZR*C2R-RZI*C2I)
    S2I = S1I + (FN+FNF)*(RZR*C2I+RZI*C2R)
    S1R = C2R
    S1I = C2I
    FN = FN - 1.0D0
    C2R = S2R*CSR
    C2I = S2I*CSR
    CKR = C2R
    CKI = C2I
    C1R = YR(KK)
    C1I = YI(KK)
    if (KODE == 1) go to 280
    call ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF)
    NZ = NZ + NW
  280   CONTINUE
    YR(KK) = C1R*CSPNR - C1I*CSPNI + C2R
    YI(KK) = C1R*CSPNI + C1I*CSPNR + C2I
    KK = KK - 1
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    if (IFLAG >= 3) go to 290
    C2R = ABS(CKR)
    C2I = ABS(CKI)
    C2M = MAX(C2R,C2I)
    if (C2M <= ASCLE) go to 290
    IFLAG = IFLAG + 1
    ASCLE = BRY(IFLAG)
    S1R = S1R*CSR
    S1I = S1I*CSR
    S2R = CKR
    S2I = CKI
    S1R = S1R*CSSR(IFLAG)
    S1I = S1I*CSSR(IFLAG)
    S2R = S2R*CSSR(IFLAG)
    S2I = S2I*CSSR(IFLAG)
    CSR = CSRR(IFLAG)
  290 CONTINUE
  return
  300 CONTINUE
  NZ = -1
  return
end
