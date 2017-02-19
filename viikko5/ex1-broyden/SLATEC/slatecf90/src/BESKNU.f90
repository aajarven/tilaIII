subroutine BESKNU (X, FNU, KODE, N, Y, NZ)
!
!! BESKNU is subsidiary to BESK.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (BESKNU-S, DBSKNU-D)
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!         BESKNU computes N member sequences of K Bessel functions
!         K/SUB(FNU+I-1)/(X), I=1,N for non-negative orders FNU and
!         positive X. Equations of the references are implemented on
!         small orders DNU for K/SUB(DNU)/(X) and K/SUB(DNU+1)/(X).
!         Forward recursion with the three term recursion relation
!         generates higher orders FNU+I-1, I=1,...,N. The parameter
!         KODE permits K/SUB(FNU+I-1)/(X) values or scaled values
!         EXP(X)*K/SUB(FNU+I-1)/(X), I=1,N to be returned.
!
!         To start the recursion FNU is normalized to the interval
!         -0.5 <= DNU < 0.5. A special form of the power series is
!         implemented on 0 < X <= X1 while the Miller algorithm for the
!         K Bessel function in terms of the confluent hypergeometric
!         function U(FNU+0.5,2*FNU+1,X) is implemented on X1 < X <= X2.
!         For X > X2, the asymptotic expansion for large X is used.
!         When FNU is a half odd integer, a special formula for
!         DNU=-0.5 and DNU+1.0=0.5 is used to start the recursion.
!
!         BESKNU assumes that a significant digit SINH(X) function is
!         available.
!
!     Description of Arguments
!
!         Input
!           X      - X > 0.0E0
!           FNU    - Order of initial K function, FNU >= 0.0E0
!           N      - Number of members of the sequence, N >= 1
!           KODE   - A parameter to indicate the scaling option
!                    KODE= 1  returns
!                             Y(I)=       K/SUB(FNU+I-1)/(X)
!                                  I=1,...,N
!                        = 2  returns
!                             Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X)
!                                  I=1,...,N
!
!         Output
!           Y      - A vector whose first N components contain values
!                    for the sequence
!                    Y(I)=       K/SUB(FNU+I-1)/(X), I=1,...,N or
!                    Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X), I=1,...,N
!                    depending on KODE
!           NZ     - Number of components set to zero due to
!                    underflow,
!                    NZ= 0   , Normal return
!                    NZ /= 0 , First NZ components of Y set to zero
!                              due to underflow, Y(I)=0.0E0,I=1,...,NZ
!
!     Error Conditions
!         Improper input arguments - a fatal error
!         Overflow - a fatal error
!         Underflow with KODE=1 - a non-fatal error (NZ /= 0)
!
!***SEE ALSO  BESK
!***REFERENCES  N. M. Temme, On the numerical evaluation of the modified
!                 Bessel function of the third kind, Journal of
!                 Computational Physics 19, (1975), pp. 324-337.
!***ROUTINES CALLED  GAMMA, I1MACH, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790201  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!   900727  Added EXTERNAL statement.  (WRB)
!   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BESKNU
!
  INTEGER I, IFLAG, INU, J, K, KK, KODE, KODED, N, NN, NZ
  INTEGER I1MACH
  REAL A, AK, A1, A2, B, BK, CC, CK, COEF, CX, DK, DNU, DNU2, ELIM, &
   ETEST, EX, F, FC, FHS, FK, FKS, FLRX, FMU, FNU, G1, G2, P, PI, &
   PT, P1, P2, Q, RTHPI, RX, S, SMU, SQK, ST, S1, S2, TM, TOL, T1, &
   T2, X, X1, X2, Y
  REAL GAMMA, R1MACH
  DIMENSION A(160), B(160), Y(*), CC(8)
  EXTERNAL GAMMA
  SAVE X1, X2, PI, RTHPI, CC
  DATA X1, X2 / 2.0E0, 17.0E0 /
  DATA PI,RTHPI        / 3.14159265358979E+00, 1.25331413731550E+00/
  DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8) &
                       / 5.77215664901533E-01,-4.20026350340952E-02, &
  -4.21977345555443E-02, 7.21894324666300E-03,-2.15241674114900E-04, &
  -2.01348547807000E-05, 1.13302723200000E-06, 6.11609500000000E-09/
!***FIRST EXECUTABLE STATEMENT  BESKNU
  KK = -I1MACH(12)
  ELIM = 2.303E0*(KK*R1MACH(5)-3.0E0)
  AK = R1MACH(3)
  TOL = MAX(AK,1.0E-15)
  if (X <= 0.0E0) go to 350
  if (FNU < 0.0E0) go to 360
  if (KODE < 1 .OR. KODE > 2) go to 370
  if (N < 1) go to 380
  NZ = 0
  IFLAG = 0
  KODED = KODE
  RX = 2.0E0/X
  INU = INT(FNU+0.5E0)
  DNU = FNU - INU
  if (ABS(DNU) == 0.5E0) go to 120
  DNU2 = 0.0E0
  if (ABS(DNU) < TOL) go to 10
  DNU2 = DNU*DNU
   10 CONTINUE
  if (X > X1) go to 120
!
!     SERIES FOR X <= X1
!
  A1 = 1.0E0 - DNU
  A2 = 1.0E0 + DNU
  T1 = 1.0E0/GAMMA(A1)
  T2 = 1.0E0/GAMMA(A2)
  if (ABS(DNU) > 0.1E0) go to 40
!     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
  S = CC(1)
  AK = 1.0E0
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
  G2 = (T1+T2)*0.5E0
  SMU = 1.0E0
  FC = 1.0E0
  FLRX = LOG(RX)
  FMU = DNU*FLRX
  if (DNU == 0.0E0) go to 60
  FC = DNU*PI
  FC = FC/SIN(FC)
  if (FMU /= 0.0E0) SMU = SINH(FMU)/FMU
   60 CONTINUE
  F = FC*(G1*COSH(FMU)+G2*FLRX*SMU)
  FC = EXP(FMU)
  P = 0.5E0*FC/T2
  Q = 0.5E0/(FC*T1)
  AK = 1.0E0
  CK = 1.0E0
  BK = 1.0E0
  S1 = F
  S2 = P
  if (INU > 0 .OR. N > 1) go to 90
  if (X < TOL) go to 80
  CX = X*X*0.25E0
   70 CONTINUE
  F = (AK*F+P+Q)/(BK-DNU2)
  P = P/(AK-DNU)
  Q = Q/(AK+DNU)
  CK = CK*CX/AK
  T1 = CK*F
  S1 = S1 + T1
  BK = BK + AK + AK + 1.0E0
  AK = AK + 1.0E0
  S = ABS(T1)/(1.0E0+ABS(S1))
  if (S > TOL) go to 70
   80 CONTINUE
  Y(1) = S1
  if (KODED == 1) RETURN
  Y(1) = S1*EXP(X)
  return
   90 CONTINUE
  if (X < TOL) go to 110
  CX = X*X*0.25E0
  100 CONTINUE
  F = (AK*F+P+Q)/(BK-DNU2)
  P = P/(AK-DNU)
  Q = Q/(AK+DNU)
  CK = CK*CX/AK
  T1 = CK*F
  S1 = S1 + T1
  T2 = CK*(P-AK*F)
  S2 = S2 + T2
  BK = BK + AK + AK + 1.0E0
  AK = AK + 1.0E0
  S = ABS(T1)/(1.0E0+ABS(S1)) + ABS(T2)/(1.0E0+ABS(S2))
  if (S > TOL) go to 100
  110 CONTINUE
  S2 = S2*RX
  if (KODED == 1) go to 170
  F = EXP(X)
  S1 = S1*F
  S2 = S2*F
  go to 170
  120 CONTINUE
  COEF = RTHPI/SQRT(X)
  if (KODED == 2) go to 130
  if (X > ELIM) go to 330
  COEF = COEF*EXP(-X)
  130 CONTINUE
  if (ABS(DNU) == 0.5E0) go to 340
  if (X > X2) go to 280
!
!     MILLER ALGORITHM FOR X1 < X <= X2
!
  ETEST = COS(PI*DNU)/(PI*X*TOL)
  FKS = 1.0E0
  FHS = 0.25E0
  FK = 0.0E0
  CK = X + X + 2.0E0
  P1 = 0.0E0
  P2 = 1.0E0
  K = 0
  140 CONTINUE
  K = K + 1
  FK = FK + 1.0E0
  AK = (FHS-DNU2)/(FKS+FK)
  BK = CK/(FK+1.0E0)
  PT = P2
  P2 = BK*P2 - AK*P1
  P1 = PT
  A(K) = AK
  B(K) = BK
  CK = CK + 2.0E0
  FKS = FKS + FK + FK + 1.0E0
  FHS = FHS + FK + FK
  if (ETEST > FK*P1) go to 140
  KK = K
  S = 1.0E0
  P1 = 0.0E0
  P2 = 1.0E0
  DO 150 I=1,K
    PT = P2
    P2 = (B(KK)*P2-P1)/A(KK)
    P1 = PT
    S = S + P2
    KK = KK - 1
  150 CONTINUE
  S1 = COEF*(P2/S)
  if (INU > 0 .OR. N > 1) go to 160
  go to 200
  160 CONTINUE
  S2 = S1*(X+DNU+0.5E0-P1/P2)/X
!
!     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION
!
  170 CONTINUE
  CK = (DNU+DNU+2.0E0)/X
  if (N == 1) INU = INU - 1
  if (INU > 0) go to 180
  if (N > 1) go to 200
  S1 = S2
  go to 200
  180 CONTINUE
  DO 190 I=1,INU
    ST = S2
    S2 = CK*S2 + S1
    S1 = ST
    CK = CK + RX
  190 CONTINUE
  if (N == 1) S1 = S2
  200 CONTINUE
  if (IFLAG == 1) go to 220
  Y(1) = S1
  if (N == 1) RETURN
  Y(2) = S2
  if (N == 2) RETURN
  DO 210 I=3,N
    Y(I) = CK*Y(I-1) + Y(I-2)
    CK = CK + RX
  210 CONTINUE
  return
!     IFLAG=1 CASES
  220 CONTINUE
  S = -X + LOG(S1)
  Y(1) = 0.0E0
  NZ = 1
  if (S < -ELIM) go to 230
  Y(1) = EXP(S)
  NZ = 0
  230 CONTINUE
  if (N == 1) RETURN
  S = -X + LOG(S2)
  Y(2) = 0.0E0
  NZ = NZ + 1
  if (S < -ELIM) go to 240
  NZ = NZ - 1
  Y(2) = EXP(S)
  240 CONTINUE
  if (N == 2) RETURN
  KK = 2
  if (NZ < 2) go to 260
  DO 250 I=3,N
    KK = I
    ST = S2
    S2 = CK*S2 + S1
    S1 = ST
    CK = CK + RX
    S = -X + LOG(S2)
    NZ = NZ + 1
    Y(I) = 0.0E0
    if (S < -ELIM) go to 250
    Y(I) = EXP(S)
    NZ = NZ - 1
    go to 260
  250 CONTINUE
  return
  260 CONTINUE
  if (KK == N) RETURN
  S2 = S2*CK + S1
  CK = CK + RX
  KK = KK + 1
  Y(KK) = EXP(-X+LOG(S2))
  if (KK == N) RETURN
  KK = KK + 1
  DO 270 I=KK,N
    Y(I) = CK*Y(I-1) + Y(I-2)
    CK = CK + RX
  270 CONTINUE
  return
!
!     ASYMPTOTIC EXPANSION FOR LARGE X, X > X2
!
!     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
!     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
!     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
!     RECURSION
  280 CONTINUE
  NN = 2
  if (INU == 0 .AND. N == 1) NN = 1
  DNU2 = DNU + DNU
  FMU = 0.0E0
  if (ABS(DNU2) < TOL) go to 290
  FMU = DNU2*DNU2
  290 CONTINUE
  EX = X*8.0E0
  S2 = 0.0E0
  DO 320 K=1,NN
    S1 = S2
    S = 1.0E0
    AK = 0.0E0
    CK = 1.0E0
    SQK = 1.0E0
    DK = EX
    DO 300 J=1,30
      CK = CK*(FMU-SQK)/DK
      S = S + CK
      DK = DK + EX
      AK = AK + 8.0E0
      SQK = SQK + AK
      if (ABS(CK) < TOL) go to 310
  300   CONTINUE
  310   S2 = S*COEF
    FMU = FMU + 8.0E0*DNU + 4.0E0
  320 CONTINUE
  if (NN > 1) go to 170
  S1 = S2
  go to 200
  330 CONTINUE
  KODED = 2
  IFLAG = 1
  go to 120
!
!     FNU=HALF ODD INTEGER CASE
!
  340 CONTINUE
  S1 = COEF
  S2 = COEF
  go to 170
!
!
  350 call XERMSG ('SLATEC', 'BESKNU', 'X NOT GREATER THAN ZERO', 2, 1)
  return
  360 call XERMSG ('SLATEC', 'BESKNU', 'FNU NOT ZERO OR POSITIVE', 2, &
     1)
  return
  370 call XERMSG ('SLATEC', 'BESKNU', 'KODE NOT 1 OR 2', 2, 1)
  return
  380 call XERMSG ('SLATEC', 'BESKNU', 'N NOT GREATER THAN 0', 2, 1)
  return
end
