subroutine ZASYI (ZR, ZI, FNU, KODE, N, YR, YI, NZ, RL, TOL, ELIM, &
     ALIM)
!
!! ZASYI is subsidiary to ZBESI and ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CASYI-A, ZASYI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z) >= 0.0 BY
!     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE ABS(Z) IN THE
!     REGION ABS(Z) > MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
!     NZ < 0 INDICATES AN OVERFLOW ON KODE=1.
!
!***SEE ALSO  ZBESI, ZBESK
!***ROUTINES CALLED  D1MACH, ZABS, ZDIV, ZEXP, ZMLT, ZSQRT
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   930122  Added ZEXP and ZSQRT to EXTERNAL statement.  (RWC)
!***END PROLOGUE  ZASYI
!     COMPLEX AK1,CK,CONE,CS1,CS2,CZ,CZERO,DK,EZ,P1,RZ,S2,Y,Z
  DOUBLE PRECISION AA, AEZ, AK, AK1I, AK1R, ALIM, ARG, ARM, ATOL, &
   AZ, BB, BK, CKI, CKR, CONEI, CONER, CS1I, CS1R, CS2I, CS2R, CZI, &
   CZR, DFNU, DKI, DKR, DNU2, ELIM, EZI, EZR, FDN, FNU, PI, P1I, &
   P1R, RAZ, RL, RTPI, RTR1, RZI, RZR, S, SGN, SQK, STI, STR, S2I, &
   S2R, TOL, TZI, TZR, YI, YR, ZEROI, ZEROR, ZI, ZR, D1MACH, ZABS
  INTEGER I, IB, IL, INU, J, JL, K, KODE, KODED, M, N, NN, NZ
  DIMENSION YR(N), YI(N)
  EXTERNAL ZABS, ZEXP, ZSQRT
  DATA PI, RTPI  /3.14159265358979324D0 , 0.159154943091895336D0 /
  DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
!***FIRST EXECUTABLE STATEMENT  ZASYI
  NZ = 0
  AZ = ZABS(ZR,ZI)
  ARM = 1.0D+3*D1MACH(1)
  RTR1 = SQRT(ARM)
  IL = MIN(2,N)
  DFNU = FNU + (N-IL)
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
  RAZ = 1.0D0/AZ
  STR = ZR*RAZ
  STI = -ZI*RAZ
  AK1R = RTPI*STR*RAZ
  AK1I = RTPI*STI*RAZ
  call ZSQRT(AK1R, AK1I, AK1R, AK1I)
  CZR = ZR
  CZI = ZI
  if (KODE /= 2) go to 10
  CZR = ZEROR
  CZI = ZI
   10 CONTINUE
  if (ABS(CZR) > ELIM) go to 100
  DNU2 = DFNU + DFNU
  KODED = 1
  if ((ABS(CZR) > ALIM) .AND. (N > 2)) go to 20
  KODED = 0
  call ZEXP(CZR, CZI, STR, STI)
  call ZMLT(AK1R, AK1I, STR, STI, AK1R, AK1I)
   20 CONTINUE
  FDN = 0.0D0
  if (DNU2 > RTR1) FDN = DNU2*DNU2
  EZR = ZR*8.0D0
  EZI = ZI*8.0D0
!-----------------------------------------------------------------------
!     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
!     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
!     EXPANSION FOR THE IMAGINARY PART.
!-----------------------------------------------------------------------
  AEZ = 8.0D0*AZ
  S = TOL/AEZ
  JL = RL+RL + 2
  P1R = ZEROR
  P1I = ZEROI
  if (ZI == 0.0D0) go to 30
!-----------------------------------------------------------------------
!     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
!     SIGNIFICANCE WHEN FNU OR N IS LARGE
!-----------------------------------------------------------------------
  INU = FNU
  ARG = (FNU-INU)*PI
  INU = INU + N - IL
  AK = -SIN(ARG)
  BK = COS(ARG)
  if (ZI < 0.0D0) BK = -BK
  P1R = AK
  P1I = BK
  if (MOD(INU,2) == 0) go to 30
  P1R = -P1R
  P1I = -P1I
   30 CONTINUE
  DO 70 K=1,IL
    SQK = FDN - 1.0D0
    ATOL = S*ABS(SQK)
    SGN = 1.0D0
    CS1R = CONER
    CS1I = CONEI
    CS2R = CONER
    CS2I = CONEI
    CKR = CONER
    CKI = CONEI
    AK = 0.0D0
    AA = 1.0D0
    BB = AEZ
    DKR = EZR
    DKI = EZI
    DO 40 J=1,JL
      call ZDIV(CKR, CKI, DKR, DKI, STR, STI)
      CKR = STR*SQK
      CKI = STI*SQK
      CS2R = CS2R + CKR
      CS2I = CS2I + CKI
      SGN = -SGN
      CS1R = CS1R + CKR*SGN
      CS1I = CS1I + CKI*SGN
      DKR = DKR + EZR
      DKI = DKI + EZI
      AA = AA*ABS(SQK)/BB
      BB = BB + AEZ
      AK = AK + 8.0D0
      SQK = SQK - AK
      if (AA <= ATOL) go to 50
   40   CONTINUE
    go to 110
   50   CONTINUE
    S2R = CS1R
    S2I = CS1I
    if (ZR+ZR >= ELIM) go to 60
    TZR = ZR + ZR
    TZI = ZI + ZI
    call ZEXP(-TZR, -TZI, STR, STI)
    call ZMLT(STR, STI, P1R, P1I, STR, STI)
    call ZMLT(STR, STI, CS2R, CS2I, STR, STI)
    S2R = S2R + STR
    S2I = S2I + STI
   60   CONTINUE
    FDN = FDN + 8.0D0*DFNU + 4.0D0
    P1R = -P1R
    P1I = -P1I
    M = N - IL + K
    YR(M) = S2R*AK1R - S2I*AK1I
    YI(M) = S2R*AK1I + S2I*AK1R
   70 CONTINUE
  if (N <= 2) RETURN
  NN = N
  K = NN - 2
  AK = K
  STR = ZR*RAZ
  STI = -ZI*RAZ
  RZR = (STR+STR)*RAZ
  RZI = (STI+STI)*RAZ
  IB = 3
  DO 80 I=IB,NN
    YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2)
    YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2)
    AK = AK - 1.0D0
    K = K - 1
   80 CONTINUE
  if (KODED == 0) RETURN
  call ZEXP(CZR, CZI, CKR, CKI)
  DO 90 I=1,NN
    STR = YR(I)*CKR - YI(I)*CKI
    YI(I) = YR(I)*CKI + YI(I)*CKR
    YR(I) = STR
   90 CONTINUE
  return
  100 CONTINUE
  NZ = -1
  return
  110 CONTINUE
  NZ=-2
  return
end
