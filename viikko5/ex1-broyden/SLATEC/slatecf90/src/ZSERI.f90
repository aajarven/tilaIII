subroutine ZSERI (ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM)
!
!! ZSERI is subsidiary to ZBESI and ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CSERI-A, ZSERI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z) >= 0.0 BY
!     MEANS OF THE POWER SERIES FOR LARGE ABS(Z) IN THE
!     REGION ABS(Z) <= 2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
!     NZ > 0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
!     DUE TO UNDERFLOW. NZ < 0 MEANS UNDERFLOW OCCURRED, BUT THE
!     CONDITION ABS(Z) <= 2*SQRT(FNU+1) WAS VIOLATED AND THE
!     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
!
!***SEE ALSO  ZBESI, ZBESK
!***ROUTINES CALLED  D1MACH, DGAMLN, ZABS, ZDIV, ZLOG, ZMLT, ZUCHK
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   930122  Added ZLOG to EXTERNAL statement.  (RWC)
!***END PROLOGUE  ZSERI
!     COMPLEX AK1,CK,COEF,CONE,CRSC,CSCL,CZ,CZERO,HZ,RZ,S1,S2,Y,Z
  DOUBLE PRECISION AA, ACZ, AK, AK1I, AK1R, ALIM, ARM, ASCLE, ATOL, &
   AZ, CKI, CKR, COEFI, COEFR, CONEI, CONER, CRSCR, CZI, CZR, DFNU, &
   ELIM, FNU, FNUP, HZI, HZR, RAZ, RS, RTR1, RZI, RZR, S, SS, STI, &
   STR, S1I, S1R, S2I, S2R, TOL, YI, YR, WI, WR, ZEROI, ZEROR, ZI, &
   ZR, DGAMLN, D1MACH, ZABS
  INTEGER I, IB, IDUM, IFLAG, IL, K, KODE, L, M, N, NN, NZ, NW
  DIMENSION YR(N), YI(N), WR(2), WI(2)
  EXTERNAL ZABS, ZLOG
  DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
!***FIRST EXECUTABLE STATEMENT  ZSERI
  NZ = 0
  AZ = ZABS(ZR,ZI)
  if (AZ == 0.0D0) go to 160
  ARM = 1.0D+3*D1MACH(1)
  RTR1 = SQRT(ARM)
  CRSCR = 1.0D0
  IFLAG = 0
  if (AZ < ARM) go to 150
  HZR = 0.5D0*ZR
  HZI = 0.5D0*ZI
  CZR = ZEROR
  CZI = ZEROI
  if (AZ <= RTR1) go to 10
  call ZMLT(HZR, HZI, HZR, HZI, CZR, CZI)
   10 CONTINUE
  ACZ = ZABS(CZR,CZI)
  NN = N
  call ZLOG(HZR, HZI, CKR, CKI, IDUM)
   20 CONTINUE
  DFNU = FNU + (NN-1)
  FNUP = DFNU + 1.0D0
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
  AK1R = CKR*DFNU
  AK1I = CKI*DFNU
  AK = DGAMLN(FNUP,IDUM)
  AK1R = AK1R - AK
  if (KODE == 2) AK1R = AK1R - ZR
  if (AK1R > (-ELIM)) go to 40
   30 CONTINUE
  NZ = NZ + 1
  YR(NN) = ZEROR
  YI(NN) = ZEROI
  if (ACZ > DFNU) go to 190
  NN = NN - 1
  if (NN == 0) RETURN
  go to 20
   40 CONTINUE
  if (AK1R > (-ALIM)) go to 50
  IFLAG = 1
  SS = 1.0D0/TOL
  CRSCR = TOL
  ASCLE = ARM*SS
   50 CONTINUE
  AA = EXP(AK1R)
  if (IFLAG == 1) AA = AA*SS
  COEFR = AA*COS(AK1I)
  COEFI = AA*SIN(AK1I)
  ATOL = TOL*ACZ/FNUP
  IL = MIN(2,NN)
  DO 90 I=1,IL
    DFNU = FNU + (NN-I)
    FNUP = DFNU + 1.0D0
    S1R = CONER
    S1I = CONEI
    if (ACZ < TOL*FNUP) go to 70
    AK1R = CONER
    AK1I = CONEI
    AK = FNUP + 2.0D0
    S = FNUP
    AA = 2.0D0
   60   CONTINUE
    RS = 1.0D0/S
    STR = AK1R*CZR - AK1I*CZI
    STI = AK1R*CZI + AK1I*CZR
    AK1R = STR*RS
    AK1I = STI*RS
    S1R = S1R + AK1R
    S1I = S1I + AK1I
    S = S + AK
    AK = AK + 2.0D0
    AA = AA*ACZ*RS
    if (AA > ATOL) go to 60
   70   CONTINUE
    S2R = S1R*COEFR - S1I*COEFI
    S2I = S1R*COEFI + S1I*COEFR
    WR(I) = S2R
    WI(I) = S2I
    if (IFLAG == 0) go to 80
    call ZUCHK(S2R, S2I, NW, ASCLE, TOL)
    if (NW /= 0) go to 30
   80   CONTINUE
    M = NN - I + 1
    YR(M) = S2R*CRSCR
    YI(M) = S2I*CRSCR
    if (I == IL) go to 90
    call ZDIV(COEFR, COEFI, HZR, HZI, STR, STI)
    COEFR = STR*DFNU
    COEFI = STI*DFNU
   90 CONTINUE
  if (NN <= 2) RETURN
  K = NN - 2
  AK = K
  RAZ = 1.0D0/AZ
  STR = ZR*RAZ
  STI = -ZI*RAZ
  RZR = (STR+STR)*RAZ
  RZI = (STI+STI)*RAZ
  if (IFLAG == 1) go to 120
  IB = 3
  100 CONTINUE
  DO 110 I=IB,NN
    YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2)
    YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2)
    AK = AK - 1.0D0
    K = K - 1
  110 CONTINUE
  return
!-----------------------------------------------------------------------
!     RECUR BACKWARD WITH SCALED VALUES
!-----------------------------------------------------------------------
  120 CONTINUE
!-----------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
!     UNDERFLOW LIMIT = ASCLE = D1MACH(1)*SS*1.0D+3
!-----------------------------------------------------------------------
  S1R = WR(1)
  S1I = WI(1)
  S2R = WR(2)
  S2I = WI(2)
  DO 130 L=3,NN
    CKR = S2R
    CKI = S2I
    S2R = S1R + (AK+FNU)*(RZR*CKR-RZI*CKI)
    S2I = S1I + (AK+FNU)*(RZR*CKI+RZI*CKR)
    S1R = CKR
    S1I = CKI
    CKR = S2R*CRSCR
    CKI = S2I*CRSCR
    YR(K) = CKR
    YI(K) = CKI
    AK = AK - 1.0D0
    K = K - 1
    if (ZABS(CKR,CKI) > ASCLE) go to 140
  130 CONTINUE
  return
  140 CONTINUE
  IB = L + 1
  if (IB > NN) RETURN
  go to 100
  150 CONTINUE
  NZ = N
  if (FNU == 0.0D0) NZ = NZ - 1
  160 CONTINUE
  YR(1) = ZEROR
  YI(1) = ZEROI
  if (FNU /= 0.0D0) go to 170
  YR(1) = CONER
  YI(1) = CONEI
  170 CONTINUE
  if (N == 1) RETURN
  DO 180 I=2,N
    YR(I) = ZEROR
    YI(I) = ZEROI
  180 CONTINUE
  return
!-----------------------------------------------------------------------
!     return WITH NZ < 0 if ABS(Z*Z/4) > FNU+N-NZ-1 COMPLETE
!     THE CALCULATION IN CBINU WITH N=N-ABS(NZ)
!-----------------------------------------------------------------------
  190 CONTINUE
  NZ = -NZ
  return
end
