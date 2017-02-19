subroutine ZUNI2 (ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL, &
     TOL, ELIM, ALIM)
!
!! ZUNI2 is subsidiary to ZBESI and ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNI2-A, ZUNI2-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
!     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
!     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
!
!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST /= 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
!     Y(I)=CZERO FOR I=NLAST+1,N
!
!***SEE ALSO  ZBESI, ZBESK
!***ROUTINES CALLED  D1MACH, ZABS, ZAIRY, ZUCHK, ZUNHJ, ZUOIK
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ZUNI2
!     COMPLEX AI,ARG,ASUM,BSUM,CFN,CI,CID,CIP,CONE,CRSC,CSCL,CSR,CSS,
!    *CZERO,C1,C2,DAI,PHI,RZ,S1,S2,Y,Z,ZB,ZETA1,ZETA2,ZN
  DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGI, &
   ARGR, ASCLE, ASUMI, ASUMR, BRY, BSUMI, BSUMR, CIDI, CIPI, CIPR, &
   CONER, CRSC, CSCL, CSRR, CSSR, C1R, C2I, C2M, C2R, DAII, &
   DAIR, ELIM, FN, FNU, FNUL, HPI, PHII, PHIR, RAST, RAZ, RS1, RZI, &
   RZR, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZBI, ZBR, ZEROI, &
   ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI, ZNI, ZNR, ZR, CYR, &
   CYI, D1MACH, ZABS, CAR, SAR
  INTEGER I, IFLAG, IN, INU, J, K, KODE, N, NAI, ND, NDAI, NLAST, &
   NN, NUF, NW, NZ, IDUM
  DIMENSION BRY(3), YR(N), YI(N), CIPR(4), CIPI(4), CSSR(3), &
   CSRR(3), CYR(2), CYI(2)
  EXTERNAL ZABS
  DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
  DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4), &
   CIPI(4)/ 1.0D0,0.0D0, 0.0D0,1.0D0, -1.0D0,0.0D0, 0.0D0,-1.0D0/
  DATA HPI, AIC  / &
        1.57079632679489662D+00,     1.265512123484645396D+00/
!***FIRST EXECUTABLE STATEMENT  ZUNI2
  NZ = 0
  ND = N
  NLAST = 0
!-----------------------------------------------------------------------
!     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
!     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
!     EXP(ALIM)=EXP(ELIM)*TOL
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
!-----------------------------------------------------------------------
!     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
!-----------------------------------------------------------------------
  ZNR = ZI
  ZNI = -ZR
  ZBR = ZR
  ZBI = ZI
  CIDI = -CONER
  INU = FNU
  ANG = HPI*(FNU-INU)
  C2R = COS(ANG)
  C2I = SIN(ANG)
  CAR = C2R
  SAR = C2I
  IN = INU + N - 1
  IN = MOD(IN,4) + 1
  STR = C2R*CIPR(IN) - C2I*CIPI(IN)
  C2I = C2R*CIPI(IN) + C2I*CIPR(IN)
  C2R = STR
  if (ZI > 0.0D0) go to 10
  ZNR = -ZNR
  ZBI = -ZBI
  CIDI = -CIDI
  C2I = -C2I
   10 CONTINUE
!-----------------------------------------------------------------------
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!-----------------------------------------------------------------------
  FN = MAX(FNU,1.0D0)
  call ZUNHJ(ZNR, ZNI, FN, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, &
   ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
  if (KODE == 1) go to 20
  STR = ZBR + ZETA2R
  STI = ZBI + ZETA2I
  RAST = FN/ZABS(STR,STI)
  STR = STR*RAST*RAST
  STI = -STI*RAST*RAST
  S1R = -ZETA1R + STR
  S1I = -ZETA1I + STI
  go to 30
   20 CONTINUE
  S1R = -ZETA1R + ZETA2R
  S1I = -ZETA1I + ZETA2I
   30 CONTINUE
  RS1 = S1R
  if (ABS(RS1) > ELIM) go to 150
   40 CONTINUE
  NN = MIN(2,ND)
  DO 90 I=1,NN
    FN = FNU + (ND-I)
    call ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR, PHII, ARGR, ARGI, &
     ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
    if (KODE == 1) go to 50
    STR = ZBR + ZETA2R
    STI = ZBI + ZETA2I
    RAST = FN/ZABS(STR,STI)
    STR = STR*RAST*RAST
    STI = -STI*RAST*RAST
    S1R = -ZETA1R + STR
    S1I = -ZETA1I + STI + ABS(ZI)
    go to 60
   50   CONTINUE
    S1R = -ZETA1R + ZETA2R
    S1I = -ZETA1I + ZETA2I
   60   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = S1R
    if (ABS(RS1) > ELIM) go to 120
    if (I == 1) IFLAG = 2
    if (ABS(RS1) < ALIM) go to 70
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    APHI = ZABS(PHIR,PHII)
    AARG = ZABS(ARGR,ARGI)
    RS1 = RS1 + LOG(APHI) - 0.25D0*LOG(AARG) - AIC
    if (ABS(RS1) > ELIM) go to 120
    if (I == 1) IFLAG = 1
    if (RS1 < 0.0D0) go to 70
    if (I == 1) IFLAG = 3
   70   CONTINUE
!-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!-----------------------------------------------------------------------
    call ZAIRY(ARGR, ARGI, 0, 2, AIR, AII, NAI, IDUM)
    call ZAIRY(ARGR, ARGI, 1, 2, DAIR, DAII, NDAI, IDUM)
    STR = DAIR*BSUMR - DAII*BSUMI
    STI = DAIR*BSUMI + DAII*BSUMR
    STR = STR + (AIR*ASUMR-AII*ASUMI)
    STI = STI + (AIR*ASUMI+AII*ASUMR)
    S2R = PHIR*STR - PHII*STI
    S2I = PHIR*STI + PHII*STR
    STR = EXP(S1R)*CSSR(IFLAG)
    S1R = STR*COS(S1I)
    S1I = STR*SIN(S1I)
    STR = S2R*S1R - S2I*S1I
    S2I = S2R*S1I + S2I*S1R
    S2R = STR
    if (IFLAG /= 1) go to 80
    call ZUCHK(S2R, S2I, NW, BRY(1), TOL)
    if (NW /= 0) go to 120
   80   CONTINUE
    if (ZI <= 0.0D0) S2I = -S2I
    STR = S2R*C2R - S2I*C2I
    S2I = S2R*C2I + S2I*C2R
    S2R = STR
    CYR(I) = S2R
    CYI(I) = S2I
    J = ND - I + 1
    YR(J) = S2R*CSRR(IFLAG)
    YI(J) = S2I*CSRR(IFLAG)
    STR = -C2I*CIDI
    C2I = C2R*CIDI
    C2R = STR
   90 CONTINUE
  if (ND <= 2) go to 110
  RAZ = 1.0D0/ZABS(ZR,ZI)
  STR = ZR*RAZ
  STI = -ZI*RAZ
  RZR = (STR+STR)*RAZ
  RZI = (STI+STI)*RAZ
  BRY(2) = 1.0D0/BRY(1)
  BRY(3) = D1MACH(2)
  S1R = CYR(1)
  S1I = CYI(1)
  S2R = CYR(2)
  S2I = CYI(2)
  C1R = CSRR(IFLAG)
  ASCLE = BRY(IFLAG)
  K = ND - 2
  FN = K
  DO 100 I=3,ND
    C2R = S2R
    C2I = S2I
    S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
    S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
    S1R = C2R
    S1I = C2I
    C2R = S2R*C1R
    C2I = S2I*C1R
    YR(K) = C2R
    YI(K) = C2I
    K = K - 1
    FN = FN - 1.0D0
    if (IFLAG >= 3) go to 100
    STR = ABS(C2R)
    STI = ABS(C2I)
    C2M = MAX(STR,STI)
    if (C2M <= ASCLE) go to 100
    IFLAG = IFLAG + 1
    ASCLE = BRY(IFLAG)
    S1R = S1R*C1R
    S1I = S1I*C1R
    S2R = C2R
    S2I = C2I
    S1R = S1R*CSSR(IFLAG)
    S1I = S1I*CSSR(IFLAG)
    S2R = S2R*CSSR(IFLAG)
    S2I = S2I*CSSR(IFLAG)
    C1R = CSRR(IFLAG)
  100 CONTINUE
  110 CONTINUE
  return
  120 CONTINUE
  if (RS1 > 0.0D0) go to 140
!-----------------------------------------------------------------------
!     SET UNDERFLOW AND UPDATE PARAMETERS
!-----------------------------------------------------------------------
  YR(ND) = ZEROR
  YI(ND) = ZEROI
  NZ = NZ + 1
  ND = ND - 1
  if (ND == 0) go to 110
  call ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
  if (NUF < 0) go to 140
  ND = ND - NUF
  NZ = NZ + NUF
  if (ND == 0) go to 110
  FN = FNU + (ND-1)
  if (FN < FNUL) go to 130
!      FN = CIDI
!      J = NUF + 1
!      K = MOD(J,4) + 1
!      S1R = CIPR(K)
!      S1I = CIPI(K)
!      if (FN < 0.0D0) S1I = -S1I
!      STR = C2R*S1R - C2I*S1I
!      C2I = C2R*S1I + C2I*S1R
!      C2R = STR
  IN = INU + ND - 1
  IN = MOD(IN,4) + 1
  C2R = CAR*CIPR(IN) - SAR*CIPI(IN)
  C2I = CAR*CIPI(IN) + SAR*CIPR(IN)
  if (ZI <= 0.0D0) C2I = -C2I
  go to 40
  130 CONTINUE
  NLAST = ND
  return
  140 CONTINUE
  NZ = -1
  return
  150 CONTINUE
  if (RS1 > 0.0D0) go to 140
  NZ = N
  DO 160 I=1,N
    YR(I) = ZEROR
    YI(I) = ZEROI
  160 CONTINUE
  return
end
