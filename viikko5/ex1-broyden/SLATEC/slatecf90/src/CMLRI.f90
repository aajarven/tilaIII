subroutine CMLRI (Z, FNU, KODE, N, Y, NZ, TOL)
!
!! CMLRI is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CMLRI-A, ZMLRI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY THE
!     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  GAMLN, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CMLRI
  COMPLEX CK, CNORM, CONE, CTWO, CZERO, PT, P1, P2, RZ, SUM, Y, Z
  REAL ACK, AK, AP, AT, AZ, BK, FKAP, FKK, FLAM, FNF, FNU, RHO, &
   RHO2, SCLE, TFNF, TOL, TST, X, GAMLN, R1MACH
  INTEGER I, IAZ, IDUM, IFNU, INU, ITIME, K, KK, KM, KODE, M, N, NZ
  DIMENSION Y(N)
  DATA CZERO,CONE,CTWO /(0.0E0,0.0E0),(1.0E0,0.0E0),(2.0E0,0.0E0)/
  SCLE = 1.0E+3*R1MACH(1)/TOL
!***FIRST EXECUTABLE STATEMENT  CMLRI
  NZ=0
  AZ = ABS(Z)
  X = REAL(Z)
  IAZ = AZ
  IFNU = FNU
  INU = IFNU + N - 1
  AT = IAZ + 1.0E0
  CK = CMPLX(AT,0.0E0)/Z
  RZ = CTWO/Z
  P1 = CZERO
  P2 = CONE
  ACK = (AT+1.0E0)/AZ
  RHO = ACK + SQRT(ACK*ACK-1.0E0)
  RHO2 = RHO*RHO
  TST = (RHO2+RHO2)/((RHO2-1.0E0)*(RHO-1.0E0))
  TST = TST/TOL
!-----------------------------------------------------------------------
!     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
!-----------------------------------------------------------------------
  AK = AT
  DO 10 I=1,80
    PT = P2
    P2 = P1 - CK*P2
    P1 = PT
    CK = CK + RZ
    AP = ABS(P2)
    if (AP > TST*AK*AK) go to 20
    AK = AK + 1.0E0
   10 CONTINUE
  go to 110
   20 CONTINUE
  I = I + 1
  K = 0
  if (INU < IAZ) go to 40
!-----------------------------------------------------------------------
!     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
!-----------------------------------------------------------------------
  P1 = CZERO
  P2 = CONE
  AT = INU + 1.0E0
  CK = CMPLX(AT,0.0E0)/Z
  ACK = AT/AZ
  TST = SQRT(ACK/TOL)
  ITIME = 1
  DO 30 K=1,80
    PT = P2
    P2 = P1 - CK*P2
    P1 = PT
    CK = CK + RZ
    AP = ABS(P2)
    if (AP < TST) go to 30
    if (ITIME == 2) go to 40
    ACK = ABS(CK)
    FLAM = ACK + SQRT(ACK*ACK-1.0E0)
    FKAP = AP/ABS(P1)
    RHO = MIN(FLAM,FKAP)
    TST = TST*SQRT(RHO/(RHO*RHO-1.0E0))
    ITIME = 2
   30 CONTINUE
  go to 110
   40 CONTINUE
!-----------------------------------------------------------------------
!     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
!-----------------------------------------------------------------------
  K = K + 1
  KK = MAX(I+IAZ,K+INU)
  FKK = KK
  P1 = CZERO
!-----------------------------------------------------------------------
!     SCALE P2 AND SUM BY SCLE
!-----------------------------------------------------------------------
  P2 = CMPLX(SCLE,0.0E0)
  FNF = FNU - IFNU
  TFNF = FNF + FNF
  BK = GAMLN(FKK+TFNF+1.0E0,IDUM) - GAMLN(FKK+1.0E0,IDUM) &
       -GAMLN(TFNF+1.0E0,IDUM)
  BK = EXP(BK)
  SUM = CZERO
  KM = KK - INU
  DO 50 I=1,KM
    PT = P2
    P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
    P1 = PT
    AK = 1.0E0 - TFNF/(FKK+TFNF)
    ACK = BK*AK
    SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
    BK = ACK
    FKK = FKK - 1.0E0
   50 CONTINUE
  Y(N) = P2
  if (N == 1) go to 70
  DO 60 I=2,N
    PT = P2
    P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
    P1 = PT
    AK = 1.0E0 - TFNF/(FKK+TFNF)
    ACK = BK*AK
    SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
    BK = ACK
    FKK = FKK - 1.0E0
    M = N - I + 1
    Y(M) = P2
   60 CONTINUE
   70 CONTINUE
  if (IFNU <= 0) go to 90
  DO 80 I=1,IFNU
    PT = P2
    P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
    P1 = PT
    AK = 1.0E0 - TFNF/(FKK+TFNF)
    ACK = BK*AK
    SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
    BK = ACK
    FKK = FKK - 1.0E0
   80 CONTINUE
   90 CONTINUE
  PT = Z
  if (KODE == 2) PT = PT - CMPLX(X,0.0E0)
  P1 = -CMPLX(FNF,0.0E0)*CLOG(RZ) + PT
  AP = GAMLN(1.0E0+FNF,IDUM)
  PT = P1 - CMPLX(AP,0.0E0)
!-----------------------------------------------------------------------
!     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
!     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
!-----------------------------------------------------------------------
  P2 = P2 + SUM
  AP = ABS(P2)
  P1 = CMPLX(1.0E0/AP,0.0E0)
  CK = CEXP(PT)*P1
  PT = CONJG(P2)*P1
  CNORM = CK*PT
  DO 100 I=1,N
    Y(I) = Y(I)*CNORM
  100 CONTINUE
  return
  110 CONTINUE
  NZ=-2
  return
end
