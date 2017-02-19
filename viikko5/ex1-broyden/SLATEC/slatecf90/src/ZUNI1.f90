subroutine ZUNI1 (ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL, &
     TOL, ELIM, ALIM)
!
!! ZUNI1 is subsidiary to ZBESI and ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNI1-A, ZUNI1-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
!     EXPANSION FOR I(FNU,Z) IN -PI/3 <= ARG Z <= PI/3.
!
!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST /= 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
!     Y(I)=CZERO FOR I=NLAST+1,N
!
!***SEE ALSO  ZBESI, ZBESK
!***ROUTINES CALLED  D1MACH, ZABS, ZUCHK, ZUNIK, ZUOIK
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ZUNI1
!     COMPLEX CFN,CONE,CRSC,CSCL,CSR,CSS,CWRK,CZERO,C1,C2,PHI,RZ,SUM,S1,
!    *S2,Y,Z,ZETA1,ZETA2
  DOUBLE PRECISION ALIM, APHI, ASCLE, BRY, CONER, CRSC, &
   CSCL, CSRR, CSSR, CWRKI, CWRKR, C1R, C2I, C2M, C2R, ELIM, FN, &
   FNU, FNUL, PHII, PHIR, RAST, RS1, RZI, RZR, STI, STR, SUMI, &
   SUMR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I, &
   ZETA1R, ZETA2I, ZETA2R, ZI, ZR, CYR, CYI, D1MACH, ZABS
  INTEGER I, IFLAG, INIT, K, KODE, M, N, ND, NLAST, NN, NUF, NW, NZ
  DIMENSION BRY(3), YR(N), YI(N), CWRKR(16), CWRKI(16), CSSR(3), &
   CSRR(3), CYR(2), CYI(2)
  EXTERNAL ZABS
  DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
!***FIRST EXECUTABLE STATEMENT  ZUNI1
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
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!-----------------------------------------------------------------------
  FN = MAX(FNU,1.0D0)
  INIT = 0
  call ZUNIK(ZR, ZI, FN, 1, 1, TOL, INIT, PHIR, PHII, ZETA1R, &
   ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
  if (KODE == 1) go to 10
  STR = ZR + ZETA2R
  STI = ZI + ZETA2I
  RAST = FN/ZABS(STR,STI)
  STR = STR*RAST*RAST
  STI = -STI*RAST*RAST
  S1R = -ZETA1R + STR
  S1I = -ZETA1I + STI
  go to 20
   10 CONTINUE
  S1R = -ZETA1R + ZETA2R
  S1I = -ZETA1I + ZETA2I
   20 CONTINUE
  RS1 = S1R
  if (ABS(RS1) > ELIM) go to 130
   30 CONTINUE
  NN = MIN(2,ND)
  DO 80 I=1,NN
    FN = FNU + (ND-I)
    INIT = 0
    call ZUNIK(ZR, ZI, FN, 1, 0, TOL, INIT, PHIR, PHII, ZETA1R, &
     ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
    if (KODE == 1) go to 40
    STR = ZR + ZETA2R
    STI = ZI + ZETA2I
    RAST = FN/ZABS(STR,STI)
    STR = STR*RAST*RAST
    STI = -STI*RAST*RAST
    S1R = -ZETA1R + STR
    S1I = -ZETA1I + STI + ZI
    go to 50
   40   CONTINUE
    S1R = -ZETA1R + ZETA2R
    S1I = -ZETA1I + ZETA2I
   50   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = S1R
    if (ABS(RS1) > ELIM) go to 110
    if (I == 1) IFLAG = 2
    if (ABS(RS1) < ALIM) go to 60
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    APHI = ZABS(PHIR,PHII)
    RS1 = RS1 + LOG(APHI)
    if (ABS(RS1) > ELIM) go to 110
    if (I == 1) IFLAG = 1
    if (RS1 < 0.0D0) go to 60
    if (I == 1) IFLAG = 3
   60   CONTINUE
!-----------------------------------------------------------------------
!     SCALE S1 if ABS(S1) < ASCLE
!-----------------------------------------------------------------------
    S2R = PHIR*SUMR - PHII*SUMI
    S2I = PHIR*SUMI + PHII*SUMR
    STR = EXP(S1R)*CSSR(IFLAG)
    S1R = STR*COS(S1I)
    S1I = STR*SIN(S1I)
    STR = S2R*S1R - S2I*S1I
    S2I = S2R*S1I + S2I*S1R
    S2R = STR
    if (IFLAG /= 1) go to 70
    call ZUCHK(S2R, S2I, NW, BRY(1), TOL)
    if (NW /= 0) go to 110
   70   CONTINUE
    CYR(I) = S2R
    CYI(I) = S2I
    M = ND - I + 1
    YR(M) = S2R*CSRR(IFLAG)
    YI(M) = S2I*CSRR(IFLAG)
   80 CONTINUE
  if (ND <= 2) go to 100
  RAST = 1.0D0/ZABS(ZR,ZI)
  STR = ZR*RAST
  STI = -ZI*RAST
  RZR = (STR+STR)*RAST
  RZI = (STI+STI)*RAST
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
  DO 90 I=3,ND
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
    if (IFLAG >= 3) go to 90
    STR = ABS(C2R)
    STI = ABS(C2I)
    C2M = MAX(STR,STI)
    if (C2M <= ASCLE) go to 90
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
   90 CONTINUE
  100 CONTINUE
  return
!-----------------------------------------------------------------------
!     SET UNDERFLOW AND UPDATE PARAMETERS
!-----------------------------------------------------------------------
  110 CONTINUE
  if (RS1 > 0.0D0) go to 120
  YR(ND) = ZEROR
  YI(ND) = ZEROI
  NZ = NZ + 1
  ND = ND - 1
  if (ND == 0) go to 100
  call ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
  if (NUF < 0) go to 120
  ND = ND - NUF
  NZ = NZ + NUF
  if (ND == 0) go to 100
  FN = FNU + (ND-1)
  if (FN >= FNUL) go to 30
  NLAST = ND
  return
  120 CONTINUE
  NZ = -1
  return
  130 CONTINUE
  if (RS1 > 0.0D0) go to 120
  NZ = N
  DO 140 I=1,N
    YR(I) = ZEROR
    YI(I) = ZEROI
  140 CONTINUE
  return
end
