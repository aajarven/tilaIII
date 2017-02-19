subroutine CSERI (Z, FNU, KODE, N, Y, NZ, TOL, ELIM, ALIM)
!
!! CSERI is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CSERI-A, ZSERI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z) >= 0.0 BY
!     MEANS OF THE POWER SERIES FOR LARGE ABS(Z) IN THE
!     REGION ABS(Z) <= 2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
!     NZ > 0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
!     DUE TO UNDERFLOW. NZ < 0 MEANS UNDERFLOW OCCURRED, BUT THE
!     CONDITION ABS(Z) <= 2*SQRT(FNU+1) WAS VIOLATED AND THE
!     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  CUCHK, GAMLN, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CSERI
  COMPLEX AK1, CK, COEF, CONE, CRSC, CZ, CZERO, HZ, RZ, S1, S2, W, &
   Y, Z
  REAL AA, ACZ, AK, ALIM, ARM, ASCLE, ATOL, AZ, DFNU, ELIM, FNU, &
   FNUP, RAK1, RS, RTR1, S, SS, TOL, X, GAMLN, R1MACH
  INTEGER I, IB, IDUM, IFLAG, IL, K, KODE, L, M, N, NN, NW, NZ
  DIMENSION Y(N), W(2)
  DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CSERI
  NZ = 0
  AZ = ABS(Z)
  if (AZ == 0.0E0) go to 150
  X = REAL(Z)
  ARM = 1.0E+3*R1MACH(1)
  RTR1 = SQRT(ARM)
  CRSC = CMPLX(1.0E0,0.0E0)
  IFLAG = 0
  if (AZ < ARM) go to 140
  HZ = Z*CMPLX(0.5E0,0.0E0)
  CZ = CZERO
  if (AZ > RTR1) CZ = HZ*HZ
  ACZ = ABS(CZ)
  NN = N
  CK = CLOG(HZ)
   10 CONTINUE
  DFNU = FNU + (NN-1)
  FNUP = DFNU + 1.0E0
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
  AK1 = CK*CMPLX(DFNU,0.0E0)
  AK = GAMLN(FNUP,IDUM)
  AK1 = AK1 - CMPLX(AK,0.0E0)
  if (KODE == 2) AK1 = AK1 - CMPLX(X,0.0E0)
  RAK1 = REAL(AK1)
  if (RAK1 > (-ELIM)) go to 30
   20 CONTINUE
  NZ = NZ + 1
  Y(NN) = CZERO
  if (ACZ > DFNU) go to 170
  NN = NN - 1
  if (NN == 0) RETURN
  go to 10
   30 CONTINUE
  if (RAK1 > (-ALIM)) go to 40
  IFLAG = 1
  SS = 1.0E0/TOL
  CRSC = CMPLX(TOL,0.0E0)
  ASCLE = ARM*SS
   40 CONTINUE
  AK = AIMAG(AK1)
  AA = EXP(RAK1)
  if (IFLAG == 1) AA = AA*SS
  COEF = CMPLX(AA,0.0E0)*CMPLX(COS(AK),SIN(AK))
  ATOL = TOL*ACZ/FNUP
  IL = MIN(2,NN)
  DO 80 I=1,IL
    DFNU = FNU + (NN-I)
    FNUP = DFNU + 1.0E0
    S1 = CONE
    if (ACZ < TOL*FNUP) go to 60
    AK1 = CONE
    AK = FNUP + 2.0E0
    S = FNUP
    AA = 2.0E0
   50   CONTINUE
    RS = 1.0E0/S
    AK1 = AK1*CZ*CMPLX(RS,0.0E0)
    S1 = S1 + AK1
    S = S + AK
    AK = AK + 2.0E0
    AA = AA*ACZ*RS
    if (AA > ATOL) go to 50
   60   CONTINUE
    M = NN - I + 1
    S2 = S1*COEF
    W(I) = S2
    if (IFLAG == 0) go to 70
    call CUCHK(S2, NW, ASCLE, TOL)
    if (NW /= 0) go to 20
   70   CONTINUE
    Y(M) = S2*CRSC
    if (I /= IL) COEF = COEF*CMPLX(DFNU,0.0E0)/HZ
   80 CONTINUE
  if (NN <= 2) RETURN
  K = NN - 2
  AK = K
  RZ = (CONE+CONE)/Z
  if (IFLAG == 1) go to 110
  IB = 3
   90 CONTINUE
  DO 100 I=IB,NN
    Y(K) = CMPLX(AK+FNU,0.0E0)*RZ*Y(K+1) + Y(K+2)
    AK = AK - 1.0E0
    K = K - 1
  100 CONTINUE
  return
!-----------------------------------------------------------------------
!     RECUR BACKWARD WITH SCALED VALUES
!-----------------------------------------------------------------------
  110 CONTINUE
!-----------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
!     UNDERFLOW LIMIT = ASCLE = R1MACH(1)*CSCL*1.0E+3
!-----------------------------------------------------------------------
  S1 = W(1)
  S2 = W(2)
  DO 120 L=3,NN
    CK = S2
    S2 = S1 + CMPLX(AK+FNU,0.0E0)*RZ*S2
    S1 = CK
    CK = S2*CRSC
    Y(K) = CK
    AK = AK - 1.0E0
    K = K - 1
    if (ABS(CK) > ASCLE) go to 130
  120 CONTINUE
  return
  130 CONTINUE
  IB = L + 1
  if (IB > NN) RETURN
  go to 90
  140 CONTINUE
  NZ = N
  if (FNU == 0.0E0) NZ = NZ - 1
  150 CONTINUE
  Y(1) = CZERO
  if (FNU == 0.0E0) Y(1) = CONE
  if (N == 1) RETURN
  DO 160 I=2,N
    Y(I) = CZERO
  160 CONTINUE
  return
!-----------------------------------------------------------------------
!     return WITH NZ < 0 if ABS(Z*Z/4) > FNU+N-NZ-1 COMPLETE
!     THE CALCULATION IN CBINU WITH N=N-ABS(NZ)
!-----------------------------------------------------------------------
  170 CONTINUE
  NZ = -NZ
  return
end
