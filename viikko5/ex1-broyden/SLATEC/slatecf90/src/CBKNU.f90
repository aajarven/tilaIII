subroutine CBKNU (Z, FNU, KODE, N, Y, NZ, TOL, ELIM, ALIM)
!
!! CBKNU is subsidiary to CAIRY, CBESH, CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CBKNU-A, ZBKNU-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE
!
!***SEE ALSO  CAIRY, CBESH, CBESI, CBESK
!***ROUTINES CALLED  CKSCL, CSHCH, CUCHK, GAMLN, I1MACH, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CBKNU
!
  COMPLEX CCH, CK, COEF, CONE, CRSC, CS, CSCL, CSH, CSR, CSS, CTWO, &
   CZ, CZERO, F, FMU, P, PT, P1, P2, Q, RZ, SMU, ST, S1, S2, Y, Z, &
   ZD, CELM, CY
  REAL AA, AK, ALIM, ASCLE, A1, A2, BB, BK, BRY, CAZ, CC, DNU, &
   DNU2, ELIM, ETEST, FC, FHS, FK, FKS, FNU, FPI, G1, G2, HPI, PI, &
   P2I, P2M, P2R, RK, RTHPI, R1, S, SPI, TM, TOL, TTH, T1, T2, XX, &
   YY, GAMLN, R1MACH, HELIM, ELM, XD, YD, ALAS, AS
  INTEGER I, IDUM, IFLAG, INU, K, KFLAG, KK, KMAX, KODE, KODED, N, &
   NZ, I1MACH, NW, J, IC, INUB
  DIMENSION BRY(3), CC(8), CSS(3), CSR(3), Y(N), CY(2)
!
  DATA KMAX / 30 /
  DATA R1 / 2.0E0 /
  DATA CZERO,CONE,CTWO /(0.0E0,0.0E0),(1.0E0,0.0E0),(2.0E0,0.0E0)/
!
  DATA PI, RTHPI, SPI ,HPI, FPI, TTH / &
       3.14159265358979324E0,       1.25331413731550025E0, &
       1.90985931710274403E0,       1.57079632679489662E0, &
       1.89769999331517738E0,       6.66666666666666666E-01/
!
  DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)/ &
       5.77215664901532861E-01,    -4.20026350340952355E-02, &
      -4.21977345555443367E-02,     7.21894324666309954E-03, &
      -2.15241674114950973E-04,    -2.01348547807882387E-05, &
       1.13302723198169588E-06,     6.11609510448141582E-09/
!
!***FIRST EXECUTABLE STATEMENT  CBKNU
  XX = REAL(Z)
  YY = AIMAG(Z)
  CAZ = ABS(Z)
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CRSC = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CRSC
  CSR(1) = CRSC
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = 1.0E+3*R1MACH(1)/TOL
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = R1MACH(2)
  NZ = 0
  IFLAG = 0
  KODED = KODE
  RZ = CTWO/Z
  INU = FNU+0.5E0
  DNU = FNU - INU
  if (ABS(DNU) == 0.5E0) go to 110
  DNU2 = 0.0E0
  if (ABS(DNU) > TOL) DNU2 = DNU*DNU
  if (CAZ > R1) go to 110
!-----------------------------------------------------------------------
!     SERIES FOR ABS(Z) <= R1
!-----------------------------------------------------------------------
  FC = 1.0E0
  SMU = CLOG(RZ)
  FMU = SMU*CMPLX(DNU,0.0E0)
  call CSHCH(FMU, CSH, CCH)
  if (DNU == 0.0E0) go to 10
  FC = DNU*PI
  FC = FC/SIN(FC)
  SMU = CSH*CMPLX(1.0E0/DNU,0.0E0)
   10 CONTINUE
  A2 = 1.0E0 + DNU
!-----------------------------------------------------------------------
!     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
!-----------------------------------------------------------------------
  T2 = EXP(-GAMLN(A2,IDUM))
  T1 = 1.0E0/(T2*FC)
  if (ABS(DNU) > 0.1E0) go to 40
!-----------------------------------------------------------------------
!     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
!-----------------------------------------------------------------------
  AK = 1.0E0
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
  G2 = 0.5E0*(T1+T2)*FC
  G1 = G1*FC
  F = CMPLX(G1,0.0E0)*CCH + SMU*CMPLX(G2,0.0E0)
  PT = CEXP(FMU)
  P = CMPLX(0.5E0/T2,0.0E0)*PT
  Q = CMPLX(0.5E0/T1,0.0E0)/PT
  S1 = F
  S2 = P
  AK = 1.0E0
  A1 = 1.0E0
  CK = CONE
  BK = 1.0E0 - DNU2
  if (INU > 0 .OR. N > 1) go to 80
!-----------------------------------------------------------------------
!     GENERATE K(FNU,Z), 0.0D0  <=  FNU  <  0.5D0 AND N=1
!-----------------------------------------------------------------------
  if (CAZ < TOL) go to 70
  CZ = Z*Z*CMPLX(0.25E0,0.0E0)
  T1 = 0.25E0*CAZ*CAZ
   60 CONTINUE
  F = (F*CMPLX(AK,0.0E0)+P+Q)*CMPLX(1.0E0/BK,0.0E0)
  P = P*CMPLX(1.0E0/(AK-DNU),0.0E0)
  Q = Q*CMPLX(1.0E0/(AK+DNU),0.0E0)
  RK = 1.0E0/AK
  CK = CK*CZ*CMPLX(RK,0.0)
  S1 = S1 + CK*F
  A1 = A1*T1*RK
  BK = BK + AK + AK + 1.0E0
  AK = AK + 1.0E0
  if (A1 > TOL) go to 60
   70 CONTINUE
  Y(1) = S1
  if (KODED == 1) RETURN
  Y(1) = S1*CEXP(Z)
  return
!-----------------------------------------------------------------------
!     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
!-----------------------------------------------------------------------
   80 CONTINUE
  if (CAZ < TOL) go to 100
  CZ = Z*Z*CMPLX(0.25E0,0.0E0)
  T1 = 0.25E0*CAZ*CAZ
   90 CONTINUE
  F = (F*CMPLX(AK,0.0E0)+P+Q)*CMPLX(1.0E0/BK,0.0E0)
  P = P*CMPLX(1.0E0/(AK-DNU),0.0E0)
  Q = Q*CMPLX(1.0E0/(AK+DNU),0.0E0)
  RK = 1.0E0/AK
  CK = CK*CZ*CMPLX(RK,0.0E0)
  S1 = S1 + CK*F
  S2 = S2 + CK*(P-F*CMPLX(AK,0.0E0))
  A1 = A1*T1*RK
  BK = BK + AK + AK + 1.0E0
  AK = AK + 1.0E0
  if (A1 > TOL) go to 90
  100 CONTINUE
  KFLAG = 2
  BK = REAL(SMU)
  A1 = FNU + 1.0E0
  AK = A1*ABS(BK)
  if (AK > ALIM) KFLAG = 3
  P2 = S2*CSS(KFLAG)
  S2 = P2*RZ
  S1 = S1*CSS(KFLAG)
  if (KODED == 1) go to 210
  F = CEXP(Z)
  S1 = S1*F
  S2 = S2*F
  go to 210
!-----------------------------------------------------------------------
!     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
!     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
!     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
!     RECURSION
!-----------------------------------------------------------------------
  110 CONTINUE
  COEF = CMPLX(RTHPI,0.0E0)/CSQRT(Z)
  KFLAG = 2
  if (KODED == 2) go to 120
  if (XX > ALIM) go to 290
!     BLANK LINE
  A1 = EXP(-XX)*REAL(CSS(KFLAG))
  PT = CMPLX(A1,0.0E0)*CMPLX(COS(YY),-SIN(YY))
  COEF = COEF*PT
  120 CONTINUE
  if (ABS(DNU) == 0.5E0) go to 300
!-----------------------------------------------------------------------
!     MILLER ALGORITHM FOR ABS(Z) > R1
!-----------------------------------------------------------------------
  AK = COS(PI*DNU)
  AK = ABS(AK)
  if (AK == 0.0E0) go to 300
  FHS = ABS(0.25E0-DNU2)
  if (FHS == 0.0E0) go to 300
!-----------------------------------------------------------------------
!     COMPUTE R2=F(E). if ABS(Z) >= R2, USE FORWARD RECURRENCE TO
!     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
!     12 <= E <= 60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(11))=
!     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
!-----------------------------------------------------------------------
  T1 = (I1MACH(11)-1)*R1MACH(5)*3.321928094E0
  T1 = MAX(T1,12.0E0)
  T1 = MIN(T1,60.0E0)
  T2 = TTH*T1 - 6.0E0
  if (XX /= 0.0E0) go to 130
  T1 = HPI
  go to 140
  130 CONTINUE
  T1 = ATAN(YY/XX)
  T1 = ABS(T1)
  140 CONTINUE
  if (T2 > CAZ) go to 170
!-----------------------------------------------------------------------
!     FORWARD RECURRENCE LOOP WHEN ABS(Z) >= R2
!-----------------------------------------------------------------------
  ETEST = AK/(PI*CAZ*TOL)
  FK = 1.0E0
  if (ETEST < 1.0E0) go to 180
  FKS = 2.0E0
  RK = CAZ + CAZ + 2.0E0
  A1 = 0.0E0
  A2 = 1.0E0
  DO 150 I=1,KMAX
    AK = FHS/FKS
    BK = RK/(FK+1.0E0)
    TM = A2
    A2 = BK*A2 - AK*A1
    A1 = TM
    RK = RK + 2.0E0
    FKS = FKS + FK + FK + 2.0E0
    FHS = FHS + FK + FK
    FK = FK + 1.0E0
    TM = ABS(A2)*FK
    if (ETEST < TM) go to 160
  150 CONTINUE
  go to 310
  160 CONTINUE
  FK = FK + SPI*T1*SQRT(T2/CAZ)
  FHS = ABS(0.25E0-DNU2)
  go to 180
  170 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE BACKWARD INDEX K FOR ABS(Z) < R2
!-----------------------------------------------------------------------
  A2 = SQRT(CAZ)
  AK = FPI*AK/(TOL*SQRT(A2))
  AA = 3.0E0*T1/(1.0E0+CAZ)
  BB = 14.7E0*T1/(28.0E0+CAZ)
  AK = (ALOG(AK)+CAZ*COS(AA)/(1.0E0+0.008E0*CAZ))/COS(BB)
  FK = 0.12125E0*AK*AK/CAZ + 1.5E0
  180 CONTINUE
  K = FK
!-----------------------------------------------------------------------
!     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
!-----------------------------------------------------------------------
  FK = K
  FKS = FK*FK
  P1 = CZERO
  P2 = CMPLX(TOL,0.0E0)
  CS = P2
  DO 190 I=1,K
    A1 = FKS - FK
    A2 = (FKS+FK)/(A1+FHS)
    RK = 2.0E0/(FK+1.0E0)
    T1 = (FK+XX)*RK
    T2 = YY*RK
    PT = P2
    P2 = (P2*CMPLX(T1,T2)-P1)*CMPLX(A2,0.0E0)
    P1 = PT
    CS = CS + P2
    FKS = A1 - FK + 1.0E0
    FK = FK - 1.0E0
  190 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE (P2/CS)=(P2/ABS(CS))*(CONJG(CS)/ABS(CS)) FOR BETTER
!     SCALING
!-----------------------------------------------------------------------
  TM = ABS(CS)
  PT = CMPLX(1.0E0/TM,0.0E0)
  S1 = PT*P2
  CS = CONJG(CS)*PT
  S1 = COEF*S1*CS
  if (INU > 0 .OR. N > 1) go to 200
  ZD = Z
  if ( IFLAG == 1) go to 270
  go to 240
  200 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE P1/P2=(P1/ABS(P2)*CONJG(P2)/ABS(P2) FOR SCALING
!-----------------------------------------------------------------------
  TM = ABS(P2)
  PT = CMPLX(1.0E0/TM,0.0E0)
  P1 = PT*P1
  P2 = CONJG(P2)*PT
  PT = P1*P2
  S2 = S1*(CONE+(CMPLX(DNU+0.5E0,0.0E0)-PT)/Z)
!-----------------------------------------------------------------------
!     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION WITH
!     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
!-----------------------------------------------------------------------
  210 CONTINUE
  CK = CMPLX(DNU+1.0E0,0.0E0)*RZ
  if (N == 1) INU = INU - 1
  if (INU > 0) go to 220
  if (N == 1) S1=S2
  ZD = Z
  if ( IFLAG == 1) go to 270
  go to 240
  220 CONTINUE
  INUB = 1
  if (IFLAG == 1) go to 261
  225 CONTINUE
  P1 = CSR(KFLAG)
  ASCLE = BRY(KFLAG)
  DO 230 I=INUB,INU
    ST = S2
    S2 = CK*S2 + S1
    S1 = ST
    CK = CK + RZ
    if (KFLAG >= 3) go to 230
    P2 = S2*P1
    P2R = REAL(P2)
    P2I = AIMAG(P2)
    P2R = ABS(P2R)
    P2I = ABS(P2I)
    P2M = MAX(P2R,P2I)
    if (P2M <= ASCLE) go to 230
    KFLAG = KFLAG + 1
    ASCLE = BRY(KFLAG)
    S1 = S1*P1
    S2 = P2
    S1 = S1*CSS(KFLAG)
    S2 = S2*CSS(KFLAG)
    P1 = CSR(KFLAG)
  230 CONTINUE
  if (N == 1) S1 = S2
  240 CONTINUE
  Y(1) = S1*CSR(KFLAG)
  if (N == 1) RETURN
  Y(2) = S2*CSR(KFLAG)
  if (N == 2) RETURN
  KK = 2
  250 CONTINUE
  KK = KK + 1
  if (KK > N) RETURN
  P1 = CSR(KFLAG)
  ASCLE = BRY(KFLAG)
  DO 260 I=KK,N
    P2 = S2
    S2 = CK*S2 + S1
    S1 = P2
    CK = CK + RZ
    P2 = S2*P1
    Y(I) = P2
    if (KFLAG >= 3) go to 260
    P2R = REAL(P2)
    P2I = AIMAG(P2)
    P2R = ABS(P2R)
    P2I = ABS(P2I)
    P2M = MAX(P2R,P2I)
    if (P2M <= ASCLE) go to 260
    KFLAG = KFLAG + 1
    ASCLE = BRY(KFLAG)
    S1 = S1*P1
    S2 = P2
    S1 = S1*CSS(KFLAG)
    S2 = S2*CSS(KFLAG)
    P1 = CSR(KFLAG)
  260 CONTINUE
  return
!-----------------------------------------------------------------------
!     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
!-----------------------------------------------------------------------
  261 CONTINUE
  HELIM = 0.5E0*ELIM
  ELM = EXP(-ELIM)
  CELM = CMPLX(ELM,0.0)
  ASCLE = BRY(1)
  ZD = Z
  XD = XX
  YD = YY
  IC = -1
  J = 2
  DO 262 I=1,INU
    ST = S2
    S2 = CK*S2+S1
    S1 = ST
    CK = CK+RZ
    AS = ABS(S2)
    ALAS = ALOG(AS)
    P2R = -XD+ALAS
    if ( P2R < (-ELIM)) go to 263
    P2 = -ZD+CLOG(S2)
    P2R = REAL(P2)
    P2I = AIMAG(P2)
    P2M = EXP(P2R)/TOL
    P1 = CMPLX(P2M,0.0E0)*CMPLX(COS(P2I),SIN(P2I))
    call CUCHK(P1,NW,ASCLE,TOL)
    if ( NW /= 0) go to 263
    J=3-J
    CY(J) = P1
    if ( IC == (I-1)) go to 264
    IC = I
    go to 262
  263   CONTINUE
    if ( ALAS < HELIM) go to 262
    XD = XD-ELIM
    S1 = S1*CELM
    S2 = S2*CELM
    ZD = CMPLX(XD,YD)
  262 CONTINUE
  if ( N == 1) S1 = S2
  go to 270
  264 CONTINUE
  KFLAG = 1
  INUB = I+1
  S2 = CY(J)
  J = 3 - J
  S1 = CY(J)
  if ( INUB <= INU) go to 225
  if ( N == 1) S1 = S2
  go to 240
  270 CONTINUE
  Y(1) = S1
  if (N == 1) go to 280
  Y(2) = S2
  280 CONTINUE
  ASCLE = BRY(1)
  call CKSCL(ZD, FNU, N, Y, NZ, RZ, ASCLE, TOL, ELIM)
  INU = N - NZ
  if (INU <= 0) RETURN
  KK = NZ + 1
  S1 = Y(KK)
  Y(KK) = S1*CSR(1)
  if (INU == 1) RETURN
  KK = NZ + 2
  S2 = Y(KK)
  Y(KK) = S2*CSR(1)
  if (INU == 2) RETURN
  T2 = FNU + (KK-1)
  CK = CMPLX(T2,0.0E0)*RZ
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
  S1 = COEF
  S2 = COEF
  go to 210
  310 CONTINUE
  NZ=-2
  return
end
