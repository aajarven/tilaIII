subroutine ZBKNU (ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM)
!
!! ZBKNU is subsidiary to ZAIRY, ZBESH, ZBESI and ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CBKNU-A, ZBKNU-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE.
!
!***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESK
!***ROUTINES CALLED  D1MACH, DGAMLN, I1MACH, ZABS, ZDIV, ZEXP, ZKSCL,
!                    ZLOG, ZMLT, ZSHCH, ZSQRT, ZUCHK
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   930122  Added ZEXP, ZLOG and ZSQRT to EXTERNAL statement.  (RWC)
!***END PROLOGUE  ZBKNU
!
  DOUBLE PRECISION AA, AK, ALIM, ASCLE, A1, A2, BB, BK, BRY, CAZ, &
   CBI, CBR, CC, CCHI, CCHR, CKI, CKR, COEFI, COEFR, CONEI, CONER, &
   CRSCR, CSCLR, CSHI, CSHR, CSI, CSR, CSRR, CSSR, CTWOR, &
   CZEROI, CZEROR, CZI, CZR, DNU, DNU2, DPI, ELIM, ETEST, FC, FHS, &
   FI, FK, FKS, FMUI, FMUR, FNU, FPI, FR, G1, G2, HPI, PI, PR, PTI, &
   PTR, P1I, P1R, P2I, P2M, P2R, QI, QR, RAK, RCAZ, RTHPI, RZI, &
   RZR, R1, S, SMUI, SMUR, SPI, STI, STR, S1I, S1R, S2I, S2R, TM, &
   TOL, TTH, T1, T2, YI, YR, ZI, ZR, DGAMLN, D1MACH, ZABS, ELM, &
   CELMR, ZDR, ZDI, AS, ALAS, HELIM, CYR, CYI
  INTEGER I, IFLAG, INU, K, KFLAG, KK, KMAX, KODE, KODED, N, NZ, &
   IDUM, I1MACH, J, IC, INUB, NW
  DIMENSION YR(N), YI(N), CC(8), CSSR(3), CSRR(3), BRY(3), CYR(2), &
   CYI(2)
  EXTERNAL ZABS, ZEXP, ZLOG, ZSQRT
!     COMPLEX Z,Y,A,B,RZ,SMU,FU,FMU,F,FLRZ,CZ,S1,S2,CSH,CCH
!     COMPLEX CK,P,Q,COEF,P1,P2,CBK,PT,CZERO,CONE,CTWO,ST,EZ,CS,DK
!
  DATA KMAX / 30 /
  DATA CZEROR,CZEROI,CONER,CONEI,CTWOR,R1/ &
    0.0D0 , 0.0D0 , 1.0D0 , 0.0D0 , 2.0D0 , 2.0D0 /
  DATA DPI, RTHPI, SPI ,HPI, FPI, TTH / &
       3.14159265358979324D0,       1.25331413731550025D0, &
       1.90985931710274403D0,       1.57079632679489662D0, &
       1.89769999331517738D0,       6.66666666666666666D-01/
  DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)/ &
       5.77215664901532861D-01,    -4.20026350340952355D-02, &
      -4.21977345555443367D-02,     7.21894324666309954D-03, &
      -2.15241674114950973D-04,    -2.01348547807882387D-05, &
       1.13302723198169588D-06,     6.11609510448141582D-09/
!***FIRST EXECUTABLE STATEMENT  ZBKNU
  CAZ = ZABS(ZR,ZI)
  CSCLR = 1.0D0/TOL
  CRSCR = TOL
  CSSR(1) = CSCLR
  CSSR(2) = 1.0D0
  CSSR(3) = CRSCR
  CSRR(1) = CRSCR
  CSRR(2) = 1.0D0
  CSRR(3) = CSCLR
  BRY(1) = 1.0D+3*D1MACH(1)/TOL
  BRY(2) = 1.0D0/BRY(1)
  BRY(3) = D1MACH(2)
  NZ = 0
  IFLAG = 0
  KODED = KODE
  RCAZ = 1.0D0/CAZ
  STR = ZR*RCAZ
  STI = -ZI*RCAZ
  RZR = (STR+STR)*RCAZ
  RZI = (STI+STI)*RCAZ
  INU = FNU+0.5D0
  DNU = FNU - INU
  if (ABS(DNU) == 0.5D0) go to 110
  DNU2 = 0.0D0
  if (ABS(DNU) > TOL) DNU2 = DNU*DNU
  if (CAZ > R1) go to 110
!-----------------------------------------------------------------------
!     SERIES FOR ABS(Z) <= R1
!-----------------------------------------------------------------------
  FC = 1.0D0
  call ZLOG(RZR, RZI, SMUR, SMUI, IDUM)
  FMUR = SMUR*DNU
  FMUI = SMUI*DNU
  call ZSHCH(FMUR, FMUI, CSHR, CSHI, CCHR, CCHI)
  if (DNU == 0.0D0) go to 10
  FC = DNU*DPI
  FC = FC/SIN(FC)
  SMUR = CSHR/DNU
  SMUI = CSHI/DNU
   10 CONTINUE
  A2 = 1.0D0 + DNU
!-----------------------------------------------------------------------
!     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
!-----------------------------------------------------------------------
  T2 = EXP(-DGAMLN(A2,IDUM))
  T1 = 1.0D0/(T2*FC)
  if (ABS(DNU) > 0.1D0) go to 40
!-----------------------------------------------------------------------
!     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
!-----------------------------------------------------------------------
  AK = 1.0D0
  S = CC(1)
  DO 20 K=2,8
    AK = AK*DNU2
    TM = CC(K)*AK
    S = S + TM
    if (ABS(TM) < TOL) go to 30
   20 CONTINUE
   30 G1 = -S
  go to 50
   40 CONTINUE
  G1 = (T1-T2)/(DNU+DNU)
   50 CONTINUE
  G2 = (T1+T2)*0.5D0
  FR = FC*(CCHR*G1+SMUR*G2)
  FI = FC*(CCHI*G1+SMUI*G2)
  call ZEXP(FMUR, FMUI, STR, STI)
  PR = 0.5D0*STR/T2
  PI = 0.5D0*STI/T2
  call ZDIV(0.5D0, 0.0D0, STR, STI, PTR, PTI)
  QR = PTR/T1
  QI = PTI/T1
  S1R = FR
  S1I = FI
  S2R = PR
  S2I = PI
  AK = 1.0D0
  A1 = 1.0D0
  CKR = CONER
  CKI = CONEI
  BK = 1.0D0 - DNU2
  if (INU > 0 .OR. N > 1) go to 80
!-----------------------------------------------------------------------
!     GENERATE K(FNU,Z), 0.0D0  <=  FNU  <  0.5D0 AND N=1
!-----------------------------------------------------------------------
  if (CAZ < TOL) go to 70
  call ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
  CZR = 0.25D0*CZR
  CZI = 0.25D0*CZI
  T1 = 0.25D0*CAZ*CAZ
   60 CONTINUE
  FR = (FR*AK+PR+QR)/BK
  FI = (FI*AK+PI+QI)/BK
  STR = 1.0D0/(AK-DNU)
  PR = PR*STR
  PI = PI*STR
  STR = 1.0D0/(AK+DNU)
  QR = QR*STR
  QI = QI*STR
  STR = CKR*CZR - CKI*CZI
  RAK = 1.0D0/AK
  CKI = (CKR*CZI+CKI*CZR)*RAK
  CKR = STR*RAK
  S1R = CKR*FR - CKI*FI + S1R
  S1I = CKR*FI + CKI*FR + S1I
  A1 = A1*T1*RAK
  BK = BK + AK + AK + 1.0D0
  AK = AK + 1.0D0
  if (A1 > TOL) go to 60
   70 CONTINUE
  YR(1) = S1R
  YI(1) = S1I
  if (KODED == 1) RETURN
  call ZEXP(ZR, ZI, STR, STI)
  call ZMLT(S1R, S1I, STR, STI, YR(1), YI(1))
  return
!-----------------------------------------------------------------------
!     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
!-----------------------------------------------------------------------
   80 CONTINUE
  if (CAZ < TOL) go to 100
  call ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
  CZR = 0.25D0*CZR
  CZI = 0.25D0*CZI
  T1 = 0.25D0*CAZ*CAZ
   90 CONTINUE
  FR = (FR*AK+PR+QR)/BK
  FI = (FI*AK+PI+QI)/BK
  STR = 1.0D0/(AK-DNU)
  PR = PR*STR
  PI = PI*STR
  STR = 1.0D0/(AK+DNU)
  QR = QR*STR
  QI = QI*STR
  STR = CKR*CZR - CKI*CZI
  RAK = 1.0D0/AK
  CKI = (CKR*CZI+CKI*CZR)*RAK
  CKR = STR*RAK
  S1R = CKR*FR - CKI*FI + S1R
  S1I = CKR*FI + CKI*FR + S1I
  STR = PR - FR*AK
  STI = PI - FI*AK
  S2R = CKR*STR - CKI*STI + S2R
  S2I = CKR*STI + CKI*STR + S2I
  A1 = A1*T1*RAK
  BK = BK + AK + AK + 1.0D0
  AK = AK + 1.0D0
  if (A1 > TOL) go to 90
  100 CONTINUE
  KFLAG = 2
  A1 = FNU + 1.0D0
  AK = A1*ABS(SMUR)
  if (AK > ALIM) KFLAG = 3
  STR = CSSR(KFLAG)
  P2R = S2R*STR
  P2I = S2I*STR
  call ZMLT(P2R, P2I, RZR, RZI, S2R, S2I)
  S1R = S1R*STR
  S1I = S1I*STR
  if (KODED == 1) go to 210
  call ZEXP(ZR, ZI, FR, FI)
  call ZMLT(S1R, S1I, FR, FI, S1R, S1I)
  call ZMLT(S2R, S2I, FR, FI, S2R, S2I)
  go to 210
!-----------------------------------------------------------------------
!     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
!     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
!     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
!     RECURSION
!-----------------------------------------------------------------------
  110 CONTINUE
  call ZSQRT(ZR, ZI, STR, STI)
  call ZDIV(RTHPI, CZEROI, STR, STI, COEFR, COEFI)
  KFLAG = 2
  if (KODED == 2) go to 120
  if (ZR > ALIM) go to 290
!     BLANK LINE
  STR = EXP(-ZR)*CSSR(KFLAG)
  STI = -STR*SIN(ZI)
  STR = STR*COS(ZI)
  call ZMLT(COEFR, COEFI, STR, STI, COEFR, COEFI)
  120 CONTINUE
  if (ABS(DNU) == 0.5D0) go to 300
!-----------------------------------------------------------------------
!     MILLER ALGORITHM FOR ABS(Z) > R1
!-----------------------------------------------------------------------
  AK = COS(DPI*DNU)
  AK = ABS(AK)
  if (AK == CZEROR) go to 300
  FHS = ABS(0.25D0-DNU2)
  if (FHS == CZEROR) go to 300
!-----------------------------------------------------------------------
!     COMPUTE R2=F(E). if ABS(Z) >= R2, USE FORWARD RECURRENCE TO
!     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
!     12 <= E <= 60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(14))=
!     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
!-----------------------------------------------------------------------
  T1 = I1MACH(14)-1
  T1 = T1*D1MACH(5)*3.321928094D0
  T1 = MAX(T1,12.0D0)
  T1 = MIN(T1,60.0D0)
  T2 = TTH*T1 - 6.0D0
  if (ZR /= 0.0D0) go to 130
  T1 = HPI
  go to 140
  130 CONTINUE
  T1 = DATAN(ZI/ZR)
  T1 = ABS(T1)
  140 CONTINUE
  if (T2 > CAZ) go to 170
!-----------------------------------------------------------------------
!     FORWARD RECURRENCE LOOP WHEN ABS(Z) >= R2
!-----------------------------------------------------------------------
  ETEST = AK/(DPI*CAZ*TOL)
  FK = CONER
  if (ETEST < CONER) go to 180
  FKS = CTWOR
  CKR = CAZ + CAZ + CTWOR
  P1R = CZEROR
  P2R = CONER
  DO 150 I=1,KMAX
    AK = FHS/FKS
    CBR = CKR/(FK+CONER)
    PTR = P2R
    P2R = CBR*P2R - P1R*AK
    P1R = PTR
    CKR = CKR + CTWOR
    FKS = FKS + FK + FK + CTWOR
    FHS = FHS + FK + FK
    FK = FK + CONER
    STR = ABS(P2R)*FK
    if (ETEST < STR) go to 160
  150 CONTINUE
  go to 310
  160 CONTINUE
  FK = FK + SPI*T1*SQRT(T2/CAZ)
  FHS = ABS(0.25D0-DNU2)
  go to 180
  170 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE BACKWARD INDEX K FOR ABS(Z) < R2
!-----------------------------------------------------------------------
  A2 = SQRT(CAZ)
  AK = FPI*AK/(TOL*SQRT(A2))
  AA = 3.0D0*T1/(1.0D0+CAZ)
  BB = 14.7D0*T1/(28.0D0+CAZ)
  AK = (LOG(AK)+CAZ*COS(AA)/(1.0D0+0.008D0*CAZ))/COS(BB)
  FK = 0.12125D0*AK*AK/CAZ + 1.5D0
  180 CONTINUE
!-----------------------------------------------------------------------
!     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
!-----------------------------------------------------------------------
  K = FK
  FK = K
  FKS = FK*FK
  P1R = CZEROR
  P1I = CZEROI
  P2R = TOL
  P2I = CZEROI
  CSR = P2R
  CSI = P2I
  DO 190 I=1,K
    A1 = FKS - FK
    AK = (FKS+FK)/(A1+FHS)
    RAK = 2.0D0/(FK+CONER)
    CBR = (FK+ZR)*RAK
    CBI = ZI*RAK
    PTR = P2R
    PTI = P2I
    P2R = (PTR*CBR-PTI*CBI-P1R)*AK
    P2I = (PTI*CBR+PTR*CBI-P1I)*AK
    P1R = PTR
    P1I = PTI
    CSR = CSR + P2R
    CSI = CSI + P2I
    FKS = A1 - FK + CONER
    FK = FK - CONER
  190 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE (P2/CS)=(P2/ABS(CS))*(CONJG(CS)/ABS(CS)) FOR BETTER
!     SCALING
!-----------------------------------------------------------------------
  TM = ZABS(CSR,CSI)
  PTR = 1.0D0/TM
  S1R = P2R*PTR
  S1I = P2I*PTR
  CSR = CSR*PTR
  CSI = -CSI*PTR
  call ZMLT(COEFR, COEFI, S1R, S1I, STR, STI)
  call ZMLT(STR, STI, CSR, CSI, S1R, S1I)
  if (INU > 0 .OR. N > 1) go to 200
  ZDR = ZR
  ZDI = ZI
  if ( IFLAG == 1) go to 270
  go to 240
  200 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE P1/P2=(P1/ABS(P2)*CONJG(P2)/ABS(P2) FOR SCALING
!-----------------------------------------------------------------------
  TM = ZABS(P2R,P2I)
  PTR = 1.0D0/TM
  P1R = P1R*PTR
  P1I = P1I*PTR
  P2R = P2R*PTR
  P2I = -P2I*PTR
  call ZMLT(P1R, P1I, P2R, P2I, PTR, PTI)
  STR = DNU + 0.5D0 - PTR
  STI = -PTI
  call ZDIV(STR, STI, ZR, ZI, STR, STI)
  STR = STR + 1.0D0
  call ZMLT(STR, STI, S1R, S1I, S2R, S2I)
!-----------------------------------------------------------------------
!     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
!     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
!-----------------------------------------------------------------------
  210 CONTINUE
  STR = DNU + 1.0D0
  CKR = STR*RZR
  CKI = STR*RZI
  if (N == 1) INU = INU - 1
  if (INU > 0) go to 220
  if (N > 1) go to 215
  S1R = S2R
  S1I = S2I
  215 CONTINUE
  ZDR = ZR
  ZDI = ZI
  if ( IFLAG == 1) go to 270
  go to 240
  220 CONTINUE
  INUB = 1
  if ( IFLAG == 1) go to 261
  225 CONTINUE
  P1R = CSRR(KFLAG)
  ASCLE = BRY(KFLAG)
  DO 230 I=INUB,INU
    STR = S2R
    STI = S2I
    S2R = CKR*STR - CKI*STI + S1R
    S2I = CKR*STI + CKI*STR + S1I
    S1R = STR
    S1I = STI
    CKR = CKR + RZR
    CKI = CKI + RZI
    if (KFLAG >= 3) go to 230
    P2R = S2R*P1R
    P2I = S2I*P1R
    STR = ABS(P2R)
    STI = ABS(P2I)
    P2M = MAX(STR,STI)
    if (P2M <= ASCLE) go to 230
    KFLAG = KFLAG + 1
    ASCLE = BRY(KFLAG)
    S1R = S1R*P1R
    S1I = S1I*P1R
    S2R = P2R
    S2I = P2I
    STR = CSSR(KFLAG)
    S1R = S1R*STR
    S1I = S1I*STR
    S2R = S2R*STR
    S2I = S2I*STR
    P1R = CSRR(KFLAG)
  230 CONTINUE
  if (N /= 1) go to 240
  S1R = S2R
  S1I = S2I
  240 CONTINUE
  STR = CSRR(KFLAG)
  YR(1) = S1R*STR
  YI(1) = S1I*STR
  if (N == 1) RETURN
  YR(2) = S2R*STR
  YI(2) = S2I*STR
  if (N == 2) RETURN
  KK = 2
  250 CONTINUE
  KK = KK + 1
  if (KK > N) RETURN
  P1R = CSRR(KFLAG)
  ASCLE = BRY(KFLAG)
  DO 260 I=KK,N
    P2R = S2R
    P2I = S2I
    S2R = CKR*P2R - CKI*P2I + S1R
    S2I = CKI*P2R + CKR*P2I + S1I
    S1R = P2R
    S1I = P2I
    CKR = CKR + RZR
    CKI = CKI + RZI
    P2R = S2R*P1R
    P2I = S2I*P1R
    YR(I) = P2R
    YI(I) = P2I
    if (KFLAG >= 3) go to 260
    STR = ABS(P2R)
    STI = ABS(P2I)
    P2M = MAX(STR,STI)
    if (P2M <= ASCLE) go to 260
    KFLAG = KFLAG + 1
    ASCLE = BRY(KFLAG)
    S1R = S1R*P1R
    S1I = S1I*P1R
    S2R = P2R
    S2I = P2I
    STR = CSSR(KFLAG)
    S1R = S1R*STR
    S1I = S1I*STR
    S2R = S2R*STR
    S2I = S2I*STR
    P1R = CSRR(KFLAG)
  260 CONTINUE
  return
!-----------------------------------------------------------------------
!     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
!-----------------------------------------------------------------------
  261 CONTINUE
  HELIM = 0.5D0*ELIM
  ELM = EXP(-ELIM)
  CELMR = ELM
  ASCLE = BRY(1)
  ZDR = ZR
  ZDI = ZI
  IC = -1
  J = 2
  DO 262 I=1,INU
    STR = S2R
    STI = S2I
    S2R = STR*CKR-STI*CKI+S1R
    S2I = STI*CKR+STR*CKI+S1I
    S1R = STR
    S1I = STI
    CKR = CKR+RZR
    CKI = CKI+RZI
    AS = ZABS(S2R,S2I)
    ALAS = LOG(AS)
    P2R = -ZDR+ALAS
    if ( P2R < (-ELIM)) go to 263
    call ZLOG(S2R,S2I,STR,STI,IDUM)
    P2R = -ZDR+STR
    P2I = -ZDI+STI
    P2M = EXP(P2R)/TOL
    P1R = P2M*COS(P2I)
    P1I = P2M*SIN(P2I)
    call ZUCHK(P1R,P1I,NW,ASCLE,TOL)
    if ( NW /= 0) go to 263
    J = 3 - J
    CYR(J) = P1R
    CYI(J) = P1I
    if ( IC == (I-1)) go to 264
    IC = I
    go to 262
  263   CONTINUE
    if ( ALAS < HELIM) go to 262
    ZDR = ZDR-ELIM
    S1R = S1R*CELMR
    S1I = S1I*CELMR
    S2R = S2R*CELMR
    S2I = S2I*CELMR
  262 CONTINUE
  if ( N /= 1) go to 270
  S1R = S2R
  S1I = S2I
  go to 270
  264 CONTINUE
  KFLAG = 1
  INUB = I+1
  S2R = CYR(J)
  S2I = CYI(J)
  J = 3 - J
  S1R = CYR(J)
  S1I = CYI(J)
  if ( INUB <= INU) go to 225
  if ( N /= 1) go to 240
  S1R = S2R
  S1I = S2I
  go to 240
  270 CONTINUE
  YR(1) = S1R
  YI(1) = S1I
  if ( N == 1) go to 280
  YR(2) = S2R
  YI(2) = S2I
  280 CONTINUE
  ASCLE = BRY(1)
  call ZKSCL(ZDR,ZDI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM)
  INU = N - NZ
  if (INU <= 0) RETURN
  KK = NZ + 1
  S1R = YR(KK)
  S1I = YI(KK)
  YR(KK) = S1R*CSRR(1)
  YI(KK) = S1I*CSRR(1)
  if (INU == 1) RETURN
  KK = NZ + 2
  S2R = YR(KK)
  S2I = YI(KK)
  YR(KK) = S2R*CSRR(1)
  YI(KK) = S2I*CSRR(1)
  if (INU == 2) RETURN
  T2 = FNU + (KK-1)
  CKR = T2*RZR
  CKI = T2*RZI
  KFLAG = 1
  go to 250
  290 CONTINUE
!-----------------------------------------------------------------------
!     SCALE BY EXP(Z), IFLAG = 1 CASES
!-----------------------------------------------------------------------
  KODED = 2
  IFLAG = 1
  KFLAG = 2
  go to 120
!-----------------------------------------------------------------------
!     FNU=HALF ODD INTEGER CASE, DNU=-0.5
!-----------------------------------------------------------------------
  300 CONTINUE
  S1R = COEFR
  S1I = COEFI
  S2R = COEFR
  S2I = COEFI
  go to 210
!
!
  310 CONTINUE
  NZ=-2
  return
end
