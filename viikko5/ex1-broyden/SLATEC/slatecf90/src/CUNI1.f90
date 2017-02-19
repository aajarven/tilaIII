subroutine CUNI1 (Z, FNU, KODE, N, Y, NZ, NLAST, FNUL, TOL, ELIM, &
     ALIM)
!
!! CUNI1 is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNI1-A, ZUNI1-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
!     EXPANSION FOR I(FNU,Z) IN -PI/3 <= ARG Z <= PI/3.
!
!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST /= 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
!     Y(I)=CZERO FOR I=NLAST+1,N
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  CUCHK, CUNIK, CUOIK, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUNI1
  COMPLEX CFN, CONE, CRSC, CSCL, CSR, CSS, CWRK, CZERO, C1, C2, &
   PHI, RZ, SUM, S1, S2, Y, Z, ZETA1, ZETA2, CY
  REAL ALIM, APHI, ASCLE, BRY, C2I, C2M, C2R, ELIM, FN, FNU, FNUL, &
   RS1, TOL, YY, R1MACH
  INTEGER I, IFLAG, INIT, K, KODE, M, N, ND, NLAST, NN, NUF, NW, NZ
  DIMENSION BRY(3), Y(N), CWRK(16), CSS(3), CSR(3), CY(2)
  DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CUNI1
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
!-----------------------------------------------------------------------
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!-----------------------------------------------------------------------
  FN = MAX(FNU,1.0E0)
  INIT = 0
  call CUNIK(Z, FN, 1, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM, CWRK)
  if (KODE == 1) go to 10
  CFN = CMPLX(FN,0.0E0)
  S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2))
  go to 20
   10 CONTINUE
  S1 = -ZETA1 + ZETA2
   20 CONTINUE
  RS1 = REAL(S1)
  if (ABS(RS1) > ELIM) go to 130
   30 CONTINUE
  NN = MIN(2,ND)
  DO 80 I=1,NN
    FN = FNU + (ND-I)
    INIT = 0
    call CUNIK(Z, FN, 1, 0, TOL, INIT, PHI, ZETA1, ZETA2, SUM, CWRK)
    if (KODE == 1) go to 40
    CFN = CMPLX(FN,0.0E0)
    YY = AIMAG(Z)
    S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2)) + CMPLX(0.0E0,YY)
    go to 50
   40   CONTINUE
    S1 = -ZETA1 + ZETA2
   50   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = REAL(S1)
    if (ABS(RS1) > ELIM) go to 110
    if (I == 1) IFLAG = 2
    if (ABS(RS1) < ALIM) go to 60
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    APHI = ABS(PHI)
    RS1 = RS1 + ALOG(APHI)
    if (ABS(RS1) > ELIM) go to 110
    if (I == 1) IFLAG = 1
    if (RS1 < 0.0E0) go to 60
    if (I == 1) IFLAG = 3
   60   CONTINUE
!-----------------------------------------------------------------------
!     SCALE S1 if ABS(S1) < ASCLE
!-----------------------------------------------------------------------
    S2 = PHI*SUM
    C2R = REAL(S1)
    C2I = AIMAG(S1)
    C2M = EXP(C2R)*REAL(CSS(IFLAG))
    S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
    S2 = S2*S1
    if (IFLAG /= 1) go to 70
    call CUCHK(S2, NW, BRY(1), TOL)
    if (NW /= 0) go to 110
   70   CONTINUE
    M = ND - I + 1
    CY(I) = S2
    Y(M) = S2*CSR(IFLAG)
   80 CONTINUE
  if (ND <= 2) go to 100
  RZ = CMPLX(2.0E0,0.0E0)/Z
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = R1MACH(2)
  S1 = CY(1)
  S2 = CY(2)
  C1 = CSR(IFLAG)
  ASCLE = BRY(IFLAG)
  K = ND - 2
  FN = K
  DO 90 I=3,ND
    C2 = S2
    S2 = S1 + CMPLX(FNU+FN,0.0E0)*RZ*S2
    S1 = C2
    C2 = S2*C1
    Y(K) = C2
    K = K - 1
    FN = FN - 1.0E0
    if (IFLAG >= 3) go to 90
    C2R = REAL(C2)
    C2I = AIMAG(C2)
    C2R = ABS(C2R)
    C2I = ABS(C2I)
    C2M = MAX(C2R,C2I)
    if (C2M <= ASCLE) go to 90
    IFLAG = IFLAG + 1
    ASCLE = BRY(IFLAG)
    S1 = S1*C1
    S2 = C2
    S1 = S1*CSS(IFLAG)
    S2 = S2*CSS(IFLAG)
    C1 = CSR(IFLAG)
   90 CONTINUE
  100 CONTINUE
  return
!-----------------------------------------------------------------------
!     SET UNDERFLOW AND UPDATE PARAMETERS
!-----------------------------------------------------------------------
  110 CONTINUE
  if (RS1 > 0.0E0) go to 120
  Y(ND) = CZERO
  NZ = NZ + 1
  ND = ND - 1
  if (ND == 0) go to 100
  call CUOIK(Z, FNU, KODE, 1, ND, Y, NUF, TOL, ELIM, ALIM)
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
  if (RS1 > 0.0E0) go to 120
  NZ = N
  DO 140 I=1,N
    Y(I) = CZERO
  140 CONTINUE
  return
end
