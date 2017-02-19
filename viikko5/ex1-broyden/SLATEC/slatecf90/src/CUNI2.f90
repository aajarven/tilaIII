subroutine CUNI2 (Z, FNU, KODE, N, Y, NZ, NLAST, FNUL, TOL, ELIM, &
     ALIM)
!
!! CUNI2 is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNI2-A, ZUNI2-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
!     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
!     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
!
!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST /= 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
!     Y(I)=CZERO FOR I=NLAST+1,N
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  CAIRY, CUCHK, CUNHJ, CUOIK, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUNI2
  COMPLEX AI, ARG, ASUM, BSUM, CFN, CI, CID, CIP, CONE, CRSC, CSCL, &
   CSR, CSS, CY, CZERO, C1, C2, DAI, PHI, RZ, S1, S2, Y, Z, ZB, &
   ZETA1, ZETA2, ZN, ZAR
  REAL AARG, AIC, ALIM, ANG, APHI, ASCLE, AY, BRY, CAR, C2I, C2M, &
   C2R, ELIM, FN, FNU, FNUL, HPI, RS1, SAR, TOL, YY, R1MACH
  INTEGER I, IFLAG, IN, INU, J, K, KODE, N, NAI, ND, NDAI, NLAST, &
   NN, NUF, NW, NZ, IDUM
  DIMENSION BRY(3), Y(N), CIP(4), CSS(3), CSR(3), CY(2)
  DATA CZERO,CONE,CI/(0.0E0,0.0E0),(1.0E0,0.0E0),(0.0E0,1.0E0)/
  DATA CIP(1),CIP(2),CIP(3),CIP(4)/ &
   (1.0E0,0.0E0), (0.0E0,1.0E0), (-1.0E0,0.0E0), (0.0E0,-1.0E0)/
  DATA HPI, AIC  / &
        1.57079632679489662E+00,     1.265512123484645396E+00/
!***FIRST EXECUTABLE STATEMENT  CUNI2
  NZ = 0
  ND = N
  NLAST = 0
!-----------------------------------------------------------------------
!     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
!     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
!     EXP(ALIM)=EXP(ELIM)*TOL
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
  YY = AIMAG(Z)
!-----------------------------------------------------------------------
!     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
!-----------------------------------------------------------------------
  ZN = -Z*CI
  ZB = Z
  CID = -CI
  INU = FNU
  ANG = HPI*(FNU-INU)
  CAR = COS(ANG)
  SAR = SIN(ANG)
  C2 = CMPLX(CAR,SAR)
  ZAR = C2
  IN = INU + N - 1
  IN = MOD(IN,4)
  C2 = C2*CIP(IN+1)
  if (YY > 0.0E0) go to 10
  ZN = CONJG(-ZN)
  ZB = CONJG(ZB)
  CID = -CID
  C2 = CONJG(C2)
   10 CONTINUE
!-----------------------------------------------------------------------
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!-----------------------------------------------------------------------
  FN = MAX(FNU,1.0E0)
  call CUNHJ(ZN, FN, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
  if (KODE == 1) go to 20
  CFN = CMPLX(FNU,0.0E0)
  S1 = -ZETA1 + CFN*(CFN/(ZB+ZETA2))
  go to 30
   20 CONTINUE
  S1 = -ZETA1 + ZETA2
   30 CONTINUE
  RS1 = REAL(S1)
  if (ABS(RS1) > ELIM) go to 150
   40 CONTINUE
  NN = MIN(2,ND)
  DO 90 I=1,NN
    FN = FNU + (ND-I)
    call CUNHJ(ZN, FN, 0, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
    if (KODE == 1) go to 50
    CFN = CMPLX(FN,0.0E0)
    AY = ABS(YY)
    S1 = -ZETA1 + CFN*(CFN/(ZB+ZETA2)) + CMPLX(0.0E0,AY)
    go to 60
   50   CONTINUE
    S1 = -ZETA1 + ZETA2
   60   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = REAL(S1)
    if (ABS(RS1) > ELIM) go to 120
    if (I == 1) IFLAG = 2
    if (ABS(RS1) < ALIM) go to 70
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    APHI = ABS(PHI)
    AARG = ABS(ARG)
    RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
    if (ABS(RS1) > ELIM) go to 120
    if (I == 1) IFLAG = 1
    if (RS1 < 0.0E0) go to 70
    if (I == 1) IFLAG = 3
   70   CONTINUE
!-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!-----------------------------------------------------------------------
    call CAIRY(ARG, 0, 2, AI, NAI, IDUM)
    call CAIRY(ARG, 1, 2, DAI, NDAI, IDUM)
    S2 = PHI*(AI*ASUM+DAI*BSUM)
    C2R = REAL(S1)
    C2I = AIMAG(S1)
    C2M = EXP(C2R)*REAL(CSS(IFLAG))
    S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
    S2 = S2*S1
    if (IFLAG /= 1) go to 80
    call CUCHK(S2, NW, BRY(1), TOL)
    if (NW /= 0) go to 120
   80   CONTINUE
    if (YY <= 0.0E0) S2 = CONJG(S2)
    J = ND - I + 1
    S2 = S2*C2
    CY(I) = S2
    Y(J) = S2*CSR(IFLAG)
    C2 = C2*CID
   90 CONTINUE
  if (ND <= 2) go to 110
  RZ = CMPLX(2.0E0,0.0E0)/Z
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = R1MACH(2)
  S1 = CY(1)
  S2 = CY(2)
  C1 = CSR(IFLAG)
  ASCLE = BRY(IFLAG)
  K = ND - 2
  FN = K
  DO 100 I=3,ND
    C2 = S2
    S2 = S1 + CMPLX(FNU+FN,0.0E0)*RZ*S2
    S1 = C2
    C2 = S2*C1
    Y(K) = C2
    K = K - 1
    FN = FN - 1.0E0
    if (IFLAG >= 3) go to 100
    C2R = REAL(C2)
    C2I = AIMAG(C2)
    C2R = ABS(C2R)
    C2I = ABS(C2I)
    C2M = MAX(C2R,C2I)
    if (C2M <= ASCLE) go to 100
    IFLAG = IFLAG + 1
    ASCLE = BRY(IFLAG)
    S1 = S1*C1
    S2 = C2
    S1 = S1*CSS(IFLAG)
    S2 = S2*CSS(IFLAG)
    C1 = CSR(IFLAG)
  100 CONTINUE
  110 CONTINUE
  return
  120 CONTINUE
  if (RS1 > 0.0E0) go to 140
!-----------------------------------------------------------------------
!     SET UNDERFLOW AND UPDATE PARAMETERS
!-----------------------------------------------------------------------
  Y(ND) = CZERO
  NZ = NZ + 1
  ND = ND - 1
  if (ND == 0) go to 110
  call CUOIK(Z, FNU, KODE, 1, ND, Y, NUF, TOL, ELIM, ALIM)
  if (NUF < 0) go to 140
  ND = ND - NUF
  NZ = NZ + NUF
  if (ND == 0) go to 110
  FN = FNU + (ND-1)
  if (FN < FNUL) go to 130
!      FN = AIMAG(CID)
!      J = NUF + 1
!      K = MOD(J,4) + 1
!      S1 = CIP(K)
!      if (FN < 0.0E0) S1 = CONJG(S1)
!      C2 = C2*S1
  IN = INU + ND - 1
  IN = MOD(IN,4) + 1
  C2 = ZAR*CIP(IN)
  if (YY <= 0.0E0)C2=CONJG(C2)
  go to 40
  130 CONTINUE
  NLAST = ND
  return
  140 CONTINUE
  NZ = -1
  return
  150 CONTINUE
  if (RS1 > 0.0E0) go to 140
  NZ = N
  DO 160 I=1,N
    Y(I) = CZERO
  160 CONTINUE
  return
end
