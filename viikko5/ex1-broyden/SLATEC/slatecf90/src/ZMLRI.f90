subroutine ZMLRI (ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL)
!
!! ZMLRI is subsidiary to ZBESI and ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CMLRI-A, ZMLRI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY THE
!     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
!
!***SEE ALSO  ZBESI, ZBESK
!***ROUTINES CALLED  D1MACH, DGAMLN, ZABS, ZEXP, ZLOG, ZMLT
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC)
!***END PROLOGUE  ZMLRI
!     COMPLEX CK,CNORM,CONE,CTWO,CZERO,PT,P1,P2,RZ,SUM,Y,Z
  DOUBLE PRECISION ACK, AK, AP, AT, AZ, BK, CKI, CKR, CNORMI, &
   CNORMR, CONEI, CONER, FKAP, FKK, FLAM, FNF, FNU, PTI, PTR, P1I, &
   P1R, P2I, P2R, RAZ, RHO, RHO2, RZI, RZR, SCLE, STI, STR, SUMI, &
   SUMR, TFNF, TOL, TST, YI, YR, ZEROI, ZEROR, ZI, ZR, DGAMLN, &
   D1MACH, ZABS
  INTEGER I, IAZ, IDUM, IFNU, INU, ITIME, K, KK, KM, KODE, M, N, NZ
  DIMENSION YR(N), YI(N)
  EXTERNAL ZABS, ZEXP, ZLOG
  DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
!***FIRST EXECUTABLE STATEMENT  ZMLRI
  SCLE = D1MACH(1)/TOL
  NZ=0
  AZ = ZABS(ZR,ZI)
  IAZ = AZ
  IFNU = FNU
  INU = IFNU + N - 1
  AT = IAZ + 1.0D0
  RAZ = 1.0D0/AZ
  STR = ZR*RAZ
  STI = -ZI*RAZ
  CKR = STR*AT*RAZ
  CKI = STI*AT*RAZ
  RZR = (STR+STR)*RAZ
  RZI = (STI+STI)*RAZ
  P1R = ZEROR
  P1I = ZEROI
  P2R = CONER
  P2I = CONEI
  ACK = (AT+1.0D0)*RAZ
  RHO = ACK + SQRT(ACK*ACK-1.0D0)
  RHO2 = RHO*RHO
  TST = (RHO2+RHO2)/((RHO2-1.0D0)*(RHO-1.0D0))
  TST = TST/TOL
!-----------------------------------------------------------------------
!     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
!-----------------------------------------------------------------------
  AK = AT
  DO 10 I=1,80
    PTR = P2R
    PTI = P2I
    P2R = P1R - (CKR*PTR-CKI*PTI)
    P2I = P1I - (CKI*PTR+CKR*PTI)
    P1R = PTR
    P1I = PTI
    CKR = CKR + RZR
    CKI = CKI + RZI
    AP = ZABS(P2R,P2I)
    if (AP > TST*AK*AK) go to 20
    AK = AK + 1.0D0
   10 CONTINUE
  go to 110
   20 CONTINUE
  I = I + 1
  K = 0
  if (INU < IAZ) go to 40
!-----------------------------------------------------------------------
!     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
!-----------------------------------------------------------------------
  P1R = ZEROR
  P1I = ZEROI
  P2R = CONER
  P2I = CONEI
  AT = INU + 1.0D0
  STR = ZR*RAZ
  STI = -ZI*RAZ
  CKR = STR*AT*RAZ
  CKI = STI*AT*RAZ
  ACK = AT*RAZ
  TST = SQRT(ACK/TOL)
  ITIME = 1
  DO 30 K=1,80
    PTR = P2R
    PTI = P2I
    P2R = P1R - (CKR*PTR-CKI*PTI)
    P2I = P1I - (CKR*PTI+CKI*PTR)
    P1R = PTR
    P1I = PTI
    CKR = CKR + RZR
    CKI = CKI + RZI
    AP = ZABS(P2R,P2I)
    if (AP < TST) go to 30
    if (ITIME == 2) go to 40
    ACK = ZABS(CKR,CKI)
    FLAM = ACK + SQRT(ACK*ACK-1.0D0)
    FKAP = AP/ZABS(P1R,P1I)
    RHO = MIN(FLAM,FKAP)
    TST = TST*SQRT(RHO/(RHO*RHO-1.0D0))
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
  P1R = ZEROR
  P1I = ZEROI
!-----------------------------------------------------------------------
!     SCALE P2 AND SUM BY SCLE
!-----------------------------------------------------------------------
  P2R = SCLE
  P2I = ZEROI
  FNF = FNU - IFNU
  TFNF = FNF + FNF
  BK = DGAMLN(FKK+TFNF+1.0D0,IDUM) - DGAMLN(FKK+1.0D0,IDUM) - &
   DGAMLN(TFNF+1.0D0,IDUM)
  BK = EXP(BK)
  SUMR = ZEROR
  SUMI = ZEROI
  KM = KK - INU
  DO 50 I=1,KM
    PTR = P2R
    PTI = P2I
    P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
    P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
    P1R = PTR
    P1I = PTI
    AK = 1.0D0 - TFNF/(FKK+TFNF)
    ACK = BK*AK
    SUMR = SUMR + (ACK+BK)*P1R
    SUMI = SUMI + (ACK+BK)*P1I
    BK = ACK
    FKK = FKK - 1.0D0
   50 CONTINUE
  YR(N) = P2R
  YI(N) = P2I
  if (N == 1) go to 70
  DO 60 I=2,N
    PTR = P2R
    PTI = P2I
    P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
    P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
    P1R = PTR
    P1I = PTI
    AK = 1.0D0 - TFNF/(FKK+TFNF)
    ACK = BK*AK
    SUMR = SUMR + (ACK+BK)*P1R
    SUMI = SUMI + (ACK+BK)*P1I
    BK = ACK
    FKK = FKK - 1.0D0
    M = N - I + 1
    YR(M) = P2R
    YI(M) = P2I
   60 CONTINUE
   70 CONTINUE
  if (IFNU <= 0) go to 90
  DO 80 I=1,IFNU
    PTR = P2R
    PTI = P2I
    P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
    P2I = P1I + (FKK+FNF)*(RZR*PTI+RZI*PTR)
    P1R = PTR
    P1I = PTI
    AK = 1.0D0 - TFNF/(FKK+TFNF)
    ACK = BK*AK
    SUMR = SUMR + (ACK+BK)*P1R
    SUMI = SUMI + (ACK+BK)*P1I
    BK = ACK
    FKK = FKK - 1.0D0
   80 CONTINUE
   90 CONTINUE
  PTR = ZR
  PTI = ZI
  if (KODE == 2) PTR = ZEROR
  call ZLOG(RZR, RZI, STR, STI, IDUM)
  P1R = -FNF*STR + PTR
  P1I = -FNF*STI + PTI
  AP = DGAMLN(1.0D0+FNF,IDUM)
  PTR = P1R - AP
  PTI = P1I
!-----------------------------------------------------------------------
!     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
!     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
!-----------------------------------------------------------------------
  P2R = P2R + SUMR
  P2I = P2I + SUMI
  AP = ZABS(P2R,P2I)
  P1R = 1.0D0/AP
  call ZEXP(PTR, PTI, STR, STI)
  CKR = STR*P1R
  CKI = STI*P1R
  PTR = P2R*P1R
  PTI = -P2I*P1R
  call ZMLT(CKR, CKI, PTR, PTI, CNORMR, CNORMI)
  DO 100 I=1,N
    STR = YR(I)*CNORMR - YI(I)*CNORMI
    YI(I) = YR(I)*CNORMI + YI(I)*CNORMR
    YR(I) = STR
  100 CONTINUE
  return
  110 CONTINUE
  NZ=-2
  return
end
