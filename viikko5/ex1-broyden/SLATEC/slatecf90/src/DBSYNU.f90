subroutine DBSYNU (X, FNU, N, Y)
!
!! DBSYNU is subsidiary to DBESY.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (BESYNU-S, DBSYNU-D)
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract  **** A DOUBLE PRECISION routine ****
!         DBSYNU computes N member sequences of Y Bessel functions
!         Y/SUB(FNU+I-1)/(X), I=1,N for non-negative orders FNU and
!         positive X. Equations of the references are implemented on
!         small orders DNU for Y/SUB(DNU)/(X) and Y/SUB(DNU+1)/(X).
!         Forward recursion with the three term recursion relation
!         generates higher orders FNU+I-1, I=1,...,N.
!
!         To start the recursion FNU is normalized to the interval
!         -0.5 <= DNU < 0.5. A special form of the power series is
!         implemented on 0 < X <= X1 while the Miller algorithm for the
!         K Bessel function in terms of the confluent hypergeometric
!         function U(FNU+0.5,2*FNU+1,I*X) is implemented on X1 < X <= X
!         Here I is the complex number SQRT(-1.).
!         For X > X2, the asymptotic expansion for large X is used.
!         When FNU is a half odd integer, a special formula for
!         DNU=-0.5 and DNU+1.0=0.5 is used to start the recursion.
!
!         The maximum number of significant digits obtainable
!         is the smaller of 14 and the number of digits carried in
!         DOUBLE PRECISION arithmetic.
!
!         DBSYNU assumes that a significant digit SINH function is
!         available.
!
!     Description of Arguments
!
!         INPUT
!           X      - X > 0.0D0
!           FNU    - Order of initial Y function, FNU >= 0.0D0
!           N      - Number of members of the sequence, N >= 1
!
!         OUTPUT
!           Y      - A vector whose first N components contain values
!                    for the sequence Y(I)=Y/SUB(FNU+I-1), I=1,N.
!
!     Error Conditions
!         Improper input arguments - a fatal error
!         Overflow - a fatal error
!
!***SEE ALSO  DBESY
!***REFERENCES  N. M. Temme, On the numerical evaluation of the ordinary
!                 Bessel function of the second kind, Journal of
!                 Computational Physics 21, (1976), pp. 343-350.
!               N. M. Temme, On the numerical evaluation of the modified
!                 Bessel function of the third kind, Journal of
!                 Computational Physics 19, (1975), pp. 324-337.
!***ROUTINES CALLED  D1MACH, DGAMMA, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800501  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!   900727  Added EXTERNAL statement.  (WRB)
!   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBSYNU
!
  INTEGER I, INU, J, K, KK, N, NN
  DOUBLE PRECISION A,AK,ARG,A1,A2,BK,CB,CBK,CC,CCK,CK,COEF,CPT, &
   CP1, CP2, CS, CS1, CS2, CX, DNU, DNU2, ETEST, ETX, F, FC, FHS, &
   FK, FKS, FLRX, FMU, FN, FNU, FX, G, G1, G2, HPI, P, PI, PT, Q, &
   RB, RBK, RCK, RELB, RPT, RP1, RP2, RS, RS1, RS2, RTHPI, RX, S, &
   SA, SB, SMU, SS, ST, S1, S2, TB, TM, TOL, T1, T2, X, X1, X2, Y
  DIMENSION A(120), RB(120), CB(120), Y(*), CC(8)
  DOUBLE PRECISION DGAMMA, D1MACH
  EXTERNAL DGAMMA
  SAVE X1, X2,PI, RTHPI, HPI, CC
  DATA X1, X2 / 3.0D0, 20.0D0 /
  DATA PI,RTHPI        / 3.14159265358979D+00, 7.97884560802865D-01/
  DATA HPI             / 1.57079632679490D+00/
  DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8) &
                       / 5.77215664901533D-01,-4.20026350340952D-02, &
  -4.21977345555443D-02, 7.21894324666300D-03,-2.15241674114900D-04, &
  -2.01348547807000D-05, 1.13302723200000D-06, 6.11609500000000D-09/
!***FIRST EXECUTABLE STATEMENT  DBSYNU
  AK = D1MACH(3)
  TOL = MAX(AK,1.0D-15)
  if (X <= 0.0D0) go to 270
  if (FNU < 0.0D0) go to 280
  if (N < 1) go to 290
  RX = 2.0D0/X
  INU = INT(FNU+0.5D0)
  DNU = FNU - INU
  if (ABS(DNU) == 0.5D0) go to 260
  DNU2 = 0.0D0
  if (ABS(DNU) < TOL) go to 10
  DNU2 = DNU*DNU
   10 CONTINUE
  if (X > X1) go to 120
!
!     SERIES FOR X <= X1
!
  A1 = 1.0D0 - DNU
  A2 = 1.0D0 + DNU
  T1 = 1.0D0/DGAMMA(A1)
  T2 = 1.0D0/DGAMMA(A2)
  if (ABS(DNU) > 0.1D0) go to 40
!     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
  S = CC(1)
  AK = 1.0D0
  DO 20 K=2,8
    AK = AK*DNU2
    TM = CC(K)*AK
    S = S + TM
    if (ABS(TM) < TOL) go to 30
   20 CONTINUE
   30 G1 = -(S+S)
  go to 50
   40 CONTINUE
  G1 = (T1-T2)/DNU
   50 CONTINUE
  G2 = T1 + T2
  SMU = 1.0D0
  FC = 1.0D0/PI
  FLRX = LOG(RX)
  FMU = DNU*FLRX
  TM = 0.0D0
  if (DNU == 0.0D0) go to 60
  TM = SIN(DNU*HPI)/DNU
  TM = (DNU+DNU)*TM*TM
  FC = DNU/SIN(DNU*PI)
  if (FMU /= 0.0D0) SMU = SINH(FMU)/FMU
   60 CONTINUE
  F = FC*(G1*COSH(FMU)+G2*FLRX*SMU)
  FX = EXP(FMU)
  P = FC*T1*FX
  Q = FC*T2/FX
  G = F + TM*Q
  AK = 1.0D0
  CK = 1.0D0
  BK = 1.0D0
  S1 = G
  S2 = P
  if (INU > 0 .OR. N > 1) go to 90
  if (X < TOL) go to 80
  CX = X*X*0.25D0
   70 CONTINUE
  F = (AK*F+P+Q)/(BK-DNU2)
  P = P/(AK-DNU)
  Q = Q/(AK+DNU)
  G = F + TM*Q
  CK = -CK*CX/AK
  T1 = CK*G
  S1 = S1 + T1
  BK = BK + AK + AK + 1.0D0
  AK = AK + 1.0D0
  S = ABS(T1)/(1.0D0+ABS(S1))
  if (S > TOL) go to 70
   80 CONTINUE
  Y(1) = -S1
  return
   90 CONTINUE
  if (X < TOL) go to 110
  CX = X*X*0.25D0
  100 CONTINUE
  F = (AK*F+P+Q)/(BK-DNU2)
  P = P/(AK-DNU)
  Q = Q/(AK+DNU)
  G = F + TM*Q
  CK = -CK*CX/AK
  T1 = CK*G
  S1 = S1 + T1
  T2 = CK*(P-AK*G)
  S2 = S2 + T2
  BK = BK + AK + AK + 1.0D0
  AK = AK + 1.0D0
  S = ABS(T1)/(1.0D0+ABS(S1)) + ABS(T2)/(1.0D0+ABS(S2))
  if (S > TOL) go to 100
  110 CONTINUE
  S2 = -S2*RX
  S1 = -S1
  go to 160
  120 CONTINUE
  COEF = RTHPI/SQRT(X)
  if (X > X2) go to 210
!
!     MILLER ALGORITHM FOR X1 < X <= X2
!
  ETEST = COS(PI*DNU)/(PI*X*TOL)
  FKS = 1.0D0
  FHS = 0.25D0
  FK = 0.0D0
  RCK = 2.0D0
  CCK = X + X
  RP1 = 0.0D0
  CP1 = 0.0D0
  RP2 = 1.0D0
  CP2 = 0.0D0
  K = 0
  130 CONTINUE
  K = K + 1
  FK = FK + 1.0D0
  AK = (FHS-DNU2)/(FKS+FK)
  PT = FK + 1.0D0
  RBK = RCK/PT
  CBK = CCK/PT
  RPT = RP2
  CPT = CP2
  RP2 = RBK*RPT - CBK*CPT - AK*RP1
  CP2 = CBK*RPT + RBK*CPT - AK*CP1
  RP1 = RPT
  CP1 = CPT
  RB(K) = RBK
  CB(K) = CBK
  A(K) = AK
  RCK = RCK + 2.0D0
  FKS = FKS + FK + FK + 1.0D0
  FHS = FHS + FK + FK
  PT = MAX(ABS(RP1),ABS(CP1))
  FC = (RP1/PT)**2 + (CP1/PT)**2
  PT = PT*SQRT(FC)*FK
  if (ETEST > PT) go to 130
  KK = K
  RS = 1.0D0
  CS = 0.0D0
  RP1 = 0.0D0
  CP1 = 0.0D0
  RP2 = 1.0D0
  CP2 = 0.0D0
  DO 140 I=1,K
    RPT = RP2
    CPT = CP2
    RP2 = (RB(KK)*RPT-CB(KK)*CPT-RP1)/A(KK)
    CP2 = (CB(KK)*RPT+RB(KK)*CPT-CP1)/A(KK)
    RP1 = RPT
    CP1 = CPT
    RS = RS + RP2
    CS = CS + CP2
    KK = KK - 1
  140 CONTINUE
  PT = MAX(ABS(RS),ABS(CS))
  FC = (RS/PT)**2 + (CS/PT)**2
  PT = PT*SQRT(FC)
  RS1 = (RP2*(RS/PT)+CP2*(CS/PT))/PT
  CS1 = (CP2*(RS/PT)-RP2*(CS/PT))/PT
  FC = HPI*(DNU-0.5D0) - X
  P = COS(FC)
  Q = SIN(FC)
  S1 = (CS1*Q-RS1*P)*COEF
  if (INU > 0 .OR. N > 1) go to 150
  Y(1) = S1
  return
  150 CONTINUE
  PT = MAX(ABS(RP2),ABS(CP2))
  FC = (RP2/PT)**2 + (CP2/PT)**2
  PT = PT*SQRT(FC)
  RPT = DNU + 0.5D0 - (RP1*(RP2/PT)+CP1*(CP2/PT))/PT
  CPT = X - (CP1*(RP2/PT)-RP1*(CP2/PT))/PT
  CS2 = CS1*CPT - RS1*RPT
  RS2 = RPT*CS1 + RS1*CPT
  S2 = (RS2*Q+CS2*P)*COEF/X
!
!     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION
!
  160 CONTINUE
  CK = (DNU+DNU+2.0D0)/X
  if (N == 1) INU = INU - 1
  if (INU > 0) go to 170
  if (N > 1) go to 190
  S1 = S2
  go to 190
  170 CONTINUE
  DO 180 I=1,INU
    ST = S2
    S2 = CK*S2 - S1
    S1 = ST
    CK = CK + RX
  180 CONTINUE
  if (N == 1) S1 = S2
  190 CONTINUE
  Y(1) = S1
  if (N == 1) RETURN
  Y(2) = S2
  if (N == 2) RETURN
  DO 200 I=3,N
    Y(I) = CK*Y(I-1) - Y(I-2)
    CK = CK + RX
  200 CONTINUE
  return
!
!     ASYMPTOTIC EXPANSION FOR LARGE X, X > X2
!
  210 CONTINUE
  NN = 2
  if (INU == 0 .AND. N == 1) NN = 1
  DNU2 = DNU + DNU
  FMU = 0.0D0
  if (ABS(DNU2) < TOL) go to 220
  FMU = DNU2*DNU2
  220 CONTINUE
  ARG = X - HPI*(DNU+0.5D0)
  SA = SIN(ARG)
  SB = COS(ARG)
  ETX = 8.0D0*X
  DO 250 K=1,NN
    S1 = S2
    T2 = (FMU-1.0D0)/ETX
    SS = T2
    RELB = TOL*ABS(T2)
    T1 = ETX
    S = 1.0D0
    FN = 1.0D0
    AK = 0.0D0
    DO 230 J=1,13
      T1 = T1 + ETX
      AK = AK + 8.0D0
      FN = FN + AK
      T2 = -T2*(FMU-FN)/T1
      S = S + T2
      T1 = T1 + ETX
      AK = AK + 8.0D0
      FN = FN + AK
      T2 = T2*(FMU-FN)/T1
      SS = SS + T2
      if (ABS(T2) <= RELB) go to 240
  230   CONTINUE
  240   S2 = COEF*(S*SA+SS*SB)
    FMU = FMU + 8.0D0*DNU + 4.0D0
    TB = SA
    SA = -SB
    SB = TB
  250 CONTINUE
  if (NN > 1) go to 160
  S1 = S2
  go to 190
!
!     FNU=HALF ODD INTEGER CASE
!
  260 CONTINUE
  COEF = RTHPI/SQRT(X)
  S1 = COEF*SIN(X)
  S2 = -COEF*COS(X)
  go to 160
!
!
  270 call XERMSG ('SLATEC', 'DBSYNU', 'X NOT GREATER THAN ZERO', 2, 1)
  return
  280 call XERMSG ('SLATEC', 'DBSYNU', 'FNU NOT ZERO OR POSITIVE', 2, &
     1)
  return
  290 call XERMSG ('SLATEC', 'DBSYNU', 'N NOT GREATER THAN 0', 2, 1)
  return
end
