subroutine CUNK2 (Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
!
!! CUNK2 is subsidiary to CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNK2-A, ZUNK2-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
!     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
!     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
!     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
!     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z if Z IS IN THE RIGHT
!     HALF PLANE OR ZR=-Z if Z IS IN THE LEFT HALF PLANE. MR INDIC-
!     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
!     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
!
!***SEE ALSO  CBESK
!***ROUTINES CALLED  CAIRY, CS1S2, CUCHK, CUNHJ, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUNK2
  COMPLEX AI, ARG, ASUM, BSUM, CFN, CI, CIP, &
   CK, CONE, CRSC, CR1, CR2, CS, CSCL, CSGN, CSPN, CSR, CSS, CY, &
   CZERO, C1, C2, DAI, PHI,  RZ, S1, S2, Y, Z, ZB, ZETA1, &
   ZETA2, ZN, ZR, PHID, ARGD, ZETA1D, ZETA2D, ASUMD, BSUMD
  REAL AARG, AIC, ALIM, ANG, APHI, ASC, ASCLE, BRY, CAR, CPN, C2I, &
   C2M, C2R, ELIM, FMR, FN, FNF, FNU, HPI, PI, RS1, SAR, SGN, SPN, &
   TOL, X, YY, R1MACH
  INTEGER I, IB, IFLAG, IFN, IL, IN, INU, IUF, K, KDFLG, KFLAG, KK, &
   KODE, MR, N, NAI, NDAI, NW, NZ, IDUM, J, IPARD, IC
  DIMENSION BRY(3), Y(N), ASUM(2), BSUM(2), PHI(2), ARG(2), &
   ZETA1(2), ZETA2(2), CY(2), CIP(4), CSS(3), CSR(3)
  DATA CZERO, CONE, CI, CR1, CR2 / &
           (0.0E0,0.0E0),(1.0E0,0.0E0),(0.0E0,1.0E0), &
  (1.0E0,1.73205080756887729E0),(-0.5E0,-8.66025403784438647E-01)/
  DATA HPI, PI, AIC / &
       1.57079632679489662E+00,     3.14159265358979324E+00, &
       1.26551212348464539E+00/
  DATA CIP(1),CIP(2),CIP(3),CIP(4)/ &
   (1.0E0,0.0E0), (0.0E0,-1.0E0), (-1.0E0,0.0E0), (0.0E0,1.0E0)/
!***FIRST EXECUTABLE STATEMENT  CUNK2
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
  YY = AIMAG(ZR)
  ZN = -ZR*CI
  ZB = ZR
  INU = FNU
  FNF = FNU - INU
  ANG = -HPI*FNF
  CAR = COS(ANG)
  SAR = SIN(ANG)
  CPN = -HPI*CAR
  SPN = -HPI*SAR
  C2 = CMPLX(-SPN,CPN)
  KK = MOD(INU,4) + 1
  CS = CR1*C2*CIP(KK)
  if (YY > 0.0E0) go to 10
  ZN = CONJG(-ZN)
  ZB = CONJG(ZB)
   10 CONTINUE
!-----------------------------------------------------------------------
!     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
!     QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE COMPUTED BY
!     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
!-----------------------------------------------------------------------
  J = 2
  DO 70 I=1,N
!-----------------------------------------------------------------------
!     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
!-----------------------------------------------------------------------
    J = 3 - J
    FN = FNU + (I-1)
    call CUNHJ(ZN, FN, 0, TOL, PHI(J), ARG(J), ZETA1(J), ZETA2(J), &
     ASUM(J), BSUM(J))
    if (KODE == 1) go to 20
    CFN = CMPLX(FN,0.0E0)
    S1 = ZETA1(J) - CFN*(CFN/(ZB+ZETA2(J)))
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
    AARG = ABS(ARG(J))
    RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
    if (ABS(RS1) > ELIM) go to 60
    if (KDFLG == 1) KFLAG = 1
    if (RS1 < 0.0E0) go to 40
    if (KDFLG == 1) KFLAG = 3
   40   CONTINUE
!-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!-----------------------------------------------------------------------
    C2 = ARG(J)*CR2
    call CAIRY(C2, 0, 2, AI, NAI, IDUM)
    call CAIRY(C2, 1, 2, DAI, NDAI, IDUM)
    S2 = CS*PHI(J)*(AI*ASUM(J)+CR2*DAI*BSUM(J))
    C2R = REAL(S1)
    C2I = AIMAG(S1)
    C2M = EXP(C2R)*REAL(CSS(KFLAG))
    S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
    S2 = S2*S1
    if (KFLAG /= 1) go to 50
    call CUCHK(S2, NW, BRY(1), TOL)
    if (NW /= 0) go to 60
   50   CONTINUE
    if (YY <= 0.0E0) S2 = CONJG(S2)
    CY(KDFLG) = S2
    Y(I) = S2*CSR(KFLAG)
    CS = -CI*CS
    if (KDFLG == 2) go to 75
    KDFLG = 2
    go to 70
   60   CONTINUE
    if (RS1 > 0.0E0) go to 300
!-----------------------------------------------------------------------
!     FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
    if (X < 0.0E0) go to 300
    KDFLG = 1
    Y(I) = CZERO
    CS = -CI*CS
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
  IB = I + 1
  if (N < IB) go to 170
!-----------------------------------------------------------------------
!     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
!     ON UNDERFLOW
!-----------------------------------------------------------------------
  FN = FNU+(N-1)
  IPARD = 1
  if (MR /= 0) IPARD = 0
  call CUNHJ(ZN,FN,IPARD,TOL,PHID,ARGD,ZETA1D,ZETA2D,ASUMD,BSUMD)
  if (KODE == 1) go to 80
  CFN=CMPLX(FN,0.0E0)
  S1=ZETA1D-CFN*(CFN/(ZB+ZETA2D))
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
  AARG = ABS(ARGD)
  RS1=RS1+ALOG(APHI)-0.25E0*ALOG(AARG)-AIC
  if (ABS(RS1) < ELIM) go to 100
   95 CONTINUE
  if (RS1 > 0.0E0) go to 300
!-----------------------------------------------------------------------
!     FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
  if (X < 0.0E0) go to 300
  NZ=N
  DO 96 I=1,N
    Y(I) = CZERO
   96 CONTINUE
  return
  100 CONTINUE
!-----------------------------------------------------------------------
!     SCALED FORWARD RECURRENCE FOR REMAINDER OF THE SEQUENCE
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
  170 CONTINUE
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
  if (YY <= 0.0E0) CSGN = CONJG(CSGN)
  IFN = INU + N - 1
  ANG = FNF*SGN
  CPN = COS(ANG)
  SPN = SIN(ANG)
  CSPN = CMPLX(CPN,SPN)
  if (MOD(IFN,2) == 1) CSPN = -CSPN
!-----------------------------------------------------------------------
!     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
!     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
!     QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE COMPUTED BY
!     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
!-----------------------------------------------------------------------
  CS = CMPLX(CAR,-SAR)*CSGN
  IN = MOD(IFN,4) + 1
  C2 = CIP(IN)
  CS = CS*CONJG(C2)
  ASC = BRY(1)
  KK = N
  KDFLG = 1
  IB = IB-1
  IC = IB-1
  IUF = 0
  DO 270 K=1,N
!-----------------------------------------------------------------------
!     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
!     FUNCTION ABOVE
!-----------------------------------------------------------------------
    FN = FNU+(KK-1)
    if (N > 2) go to 180
  175   CONTINUE
    PHID = PHI(J)
    ARGD = ARG(J)
    ZETA1D = ZETA1(J)
    ZETA2D = ZETA2(J)
    ASUMD = ASUM(J)
    BSUMD = BSUM(J)
    J = 3 - J
    go to 190
  180   CONTINUE
    if ((KK == N).AND.(IB < N)) go to 190
    if ((KK == IB).OR.(KK == IC)) go to 175
    call CUNHJ(ZN, FN, 0, TOL, PHID, ARGD, ZETA1D, ZETA2D, &
     ASUMD, BSUMD)
  190   CONTINUE
    if (KODE == 1) go to 200
    CFN = CMPLX(FN,0.0E0)
    S1 = -ZETA1D + CFN*(CFN/(ZB+ZETA2D))
    go to 210
  200   CONTINUE
    S1 = -ZETA1D + ZETA2D
  210   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = REAL(S1)
    if (ABS(RS1) > ELIM) go to 260
    if (KDFLG == 1) IFLAG = 2
    if (ABS(RS1) < ALIM) go to 220
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    APHI = ABS(PHID)
    AARG = ABS(ARGD)
    RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
    if (ABS(RS1) > ELIM) go to 260
    if (KDFLG == 1) IFLAG = 1
    if (RS1 < 0.0E0) go to 220
    if (KDFLG == 1) IFLAG = 3
  220   CONTINUE
    call CAIRY(ARGD, 0, 2, AI, NAI, IDUM)
    call CAIRY(ARGD, 1, 2, DAI, NDAI, IDUM)
    S2 = CS*PHID*(AI*ASUMD+DAI*BSUMD)
    C2R = REAL(S1)
    C2I = AIMAG(S1)
    C2M = EXP(C2R)*REAL(CSS(IFLAG))
    S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
    S2 = S2*S1
    if (IFLAG /= 1) go to 230
    call CUCHK(S2, NW, BRY(1), TOL)
    if (NW /= 0) S2 = CMPLX(0.0E0,0.0E0)
  230   CONTINUE
    if (YY <= 0.0E0) S2 = CONJG(S2)
    CY(KDFLG) = S2
    C2 = S2
    S2 = S2*CSR(IFLAG)
!-----------------------------------------------------------------------
!     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
!-----------------------------------------------------------------------
    S1 = Y(KK)
    if (KODE == 1) go to 250
    call CS1S2(ZR, S1, S2, NW, ASC, ALIM, IUF)
    NZ = NZ + NW
  250   CONTINUE
    Y(KK) = S1*CSPN + S2
    KK = KK - 1
    CSPN = -CSPN
    CS = -CS*CI
    if (C2 /= CZERO) go to 255
    KDFLG = 1
    go to 270
  255   CONTINUE
    if (KDFLG == 2) go to 275
    KDFLG = 2
    go to 270
  260   CONTINUE
    if (RS1 > 0.0E0) go to 300
    S2 = CZERO
    go to 230
  270 CONTINUE
  K = N
  275 CONTINUE
  IL = N-K
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
  FN = INU+IL
  DO 290 I=1,IL
    C2 = S2
    S2 = S1 + CMPLX(FN+FNF,0.0E0)*RZ*S2
    S1 = C2
    FN = FN - 1.0E0
    C2 = S2*CS
    CK = C2
    C1 = Y(KK)
    if (KODE == 1) go to 280
    call CS1S2(ZR, C1, C2, NW, ASC, ALIM, IUF)
    NZ = NZ + NW
  280   CONTINUE
    Y(KK) = C1*CSPN + C2
    KK = KK - 1
    CSPN = -CSPN
    if (IFLAG >= 3) go to 290
    C2R = REAL(CK)
    C2I = AIMAG(CK)
    C2R = ABS(C2R)
    C2I = ABS(C2I)
    C2M = MAX(C2R,C2I)
    if (C2M <= ASCLE) go to 290
    IFLAG = IFLAG + 1
    ASCLE = BRY(IFLAG)
    S1 = S1*CS
    S2 = CK
    S1 = S1*CSS(IFLAG)
    S2 = S2*CSS(IFLAG)
    CS = CSR(IFLAG)
  290 CONTINUE
  return
  300 CONTINUE
  NZ = -1
  return
end
