subroutine CASYI (Z, FNU, KODE, N, Y, NZ, RL, TOL, ELIM, ALIM)
!
!! CASYI is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CASYI-A, ZASYI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z) >= 0.0 BY
!     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE ABS(Z) IN THE
!     REGION ABS(Z) > MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
!     NZ < 0 INDICATES AN OVERFLOW ON KODE=1.
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CASYI
  COMPLEX AK1, CK, CONE, CS1, CS2, CZ, CZERO, DK, EZ, P1, RZ, S2, &
   Y, Z
  REAL AA, ACZ, AEZ, AK, ALIM, ARG, ARM, ATOL, AZ, BB, BK, DFNU, &
   DNU2, ELIM, FDN, FNU, PI, RL, RTPI, RTR1, S, SGN, SQK, TOL, X, &
   YY, R1MACH
  INTEGER I, IB, IL, INU, J, JL, K, KODE, KODED, M, N, NN, NZ
  DIMENSION Y(N)
  DATA PI, RTPI  /3.14159265358979324E0 , 0.159154943091895336E0 /
  DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CASYI
  NZ = 0
  AZ = ABS(Z)
  X = REAL(Z)
  ARM = 1.0E+3*R1MACH(1)
  RTR1 = SQRT(ARM)
  IL = MIN(2,N)
  DFNU = FNU + (N-IL)
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
  AK1 = CMPLX(RTPI,0.0E0)/Z
  AK1 = CSQRT(AK1)
  CZ = Z
  if (KODE == 2) CZ = Z - CMPLX(X,0.0E0)
  ACZ = REAL(CZ)
  if (ABS(ACZ) > ELIM) go to 80
  DNU2 = DFNU + DFNU
  KODED = 1
  if ((ABS(ACZ) > ALIM) .AND. (N > 2)) go to 10
  KODED = 0
  AK1 = AK1*CEXP(CZ)
   10 CONTINUE
  FDN = 0.0E0
  if (DNU2 > RTR1) FDN = DNU2*DNU2
  EZ = Z*CMPLX(8.0E0,0.0E0)
!-----------------------------------------------------------------------
!     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
!     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
!     EXPANSION FOR THE IMAGINARY PART.
!-----------------------------------------------------------------------
  AEZ = 8.0E0*AZ
  S = TOL/AEZ
  JL = RL+RL + 2
  YY = AIMAG(Z)
  P1 = CZERO
  if (YY == 0.0E0) go to 20
!-----------------------------------------------------------------------
!     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
!     SIGNIFICANCE WHEN FNU OR N IS LARGE
!-----------------------------------------------------------------------
  INU = FNU
  ARG = (FNU-INU)*PI
  INU = INU + N - IL
  AK = -SIN(ARG)
  BK = COS(ARG)
  if (YY < 0.0E0) BK = -BK
  P1 = CMPLX(AK,BK)
  if (MOD(INU,2) == 1) P1 = -P1
   20 CONTINUE
  DO 50 K=1,IL
    SQK = FDN - 1.0E0
    ATOL = S*ABS(SQK)
    SGN = 1.0E0
    CS1 = CONE
    CS2 = CONE
    CK = CONE
    AK = 0.0E0
    AA = 1.0E0
    BB = AEZ
    DK = EZ
    DO 30 J=1,JL
      CK = CK*CMPLX(SQK,0.0E0)/DK
      CS2 = CS2 + CK
      SGN = -SGN
      CS1 = CS1 + CK*CMPLX(SGN,0.0E0)
      DK = DK + EZ
      AA = AA*ABS(SQK)/BB
      BB = BB + AEZ
      AK = AK + 8.0E0
      SQK = SQK - AK
      if (AA <= ATOL) go to 40
   30   CONTINUE
    go to 90
   40   CONTINUE
    S2 = CS1
    if (X+X < ELIM) S2 = S2 + P1*CS2*CEXP(-Z-Z)
    FDN = FDN + 8.0E0*DFNU + 4.0E0
    P1 = -P1
    M = N - IL + K
    Y(M) = S2*AK1
   50 CONTINUE
  if (N <= 2) RETURN
  NN = N
  K = NN - 2
  AK = K
  RZ = (CONE+CONE)/Z
  IB = 3
  DO 60 I=IB,NN
    Y(K) = CMPLX(AK+FNU,0.0E0)*RZ*Y(K+1) + Y(K+2)
    AK = AK - 1.0E0
    K = K - 1
   60 CONTINUE
  if (KODED == 0) RETURN
  CK = CEXP(CZ)
  DO 70 I=1,NN
    Y(I) = Y(I)*CK
   70 CONTINUE
  return
   80 CONTINUE
  NZ = -1
  return
   90 CONTINUE
  NZ=-2
  return
end
