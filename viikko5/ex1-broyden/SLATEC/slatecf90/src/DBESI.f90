subroutine DBESI (X, ALPHA, KODE, N, Y, NZ)
!
!! DBESI computes an N member sequence of I Bessel functions ...
!            I/SUB(ALPHA+K-1)/(X), K=1,...,N or scaled Bessel functions
!            EXP(-X)*I/SUB(ALPHA+K-1)/(X), K=1,...,N for nonnegative
!            ALPHA and X.
!
!***LIBRARY   SLATEC
!***CATEGORY  C10B3
!***TYPE      DOUBLE PRECISION (BESI-S, DBESI-D)
!***KEYWORDS  I BESSEL FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Amos, D. E., (SNLA)
!           Daniel, S. L., (SNLA)
!***DESCRIPTION
!
!     Abstract  **** a double precision routine ****
!         DBESI computes an N member sequence of I Bessel functions
!         I/sub(ALPHA+K-1)/(X), K=1,...,N or scaled Bessel functions
!         EXP(-X)*I/sub(ALPHA+K-1)/(X), K=1,...,N for nonnegative ALPHA
!         and X.  A combination of the power series, the asymptotic
!         expansion for X to infinity, and the uniform asymptotic
!         expansion for NU to infinity are applied over subdivisions of
!         the (NU,X) plane.  For values not covered by one of these
!         formulae, the order is incremented by an integer so that one
!         of these formulae apply.  Backward recursion is used to reduce
!         orders by integer values.  The asymptotic expansion for X to
!         infinity is used only when the entire sequence (specifically
!         the last member) lies within the region covered by the
!         expansion.  Leading terms of these expansions are used to test
!         for over or underflow where appropriate.  If a sequence is
!         requested and the last member would underflow, the result is
!         set to zero and the next lower order tried, etc., until a
!         member comes on scale or all are set to zero.  An overflow
!         cannot occur with scaling.
!
!         The maximum number of significant digits obtainable
!         is the smaller of 14 and the number of digits carried in
!         double precision arithmetic.
!
!     Description of Arguments
!
!         Input      X,ALPHA are double precision
!           X      - X  >=  0.0D0
!           ALPHA  - order of first member of the sequence,
!                    ALPHA  >=  0.0D0
!           KODE   - a parameter to indicate the scaling option
!                    KODE=1 returns
!                           Y(K)=        I/sub(ALPHA+K-1)/(X),
!                                K=1,...,N
!                    KODE=2 returns
!                           Y(K)=EXP(-X)*I/sub(ALPHA+K-1)/(X),
!                                K=1,...,N
!           N      - number of members in the sequence, N  >=  1
!
!         Output     Y is double precision
!           Y      - a vector whose first N components contain
!                    values for I/sub(ALPHA+K-1)/(X) or scaled
!                    values for EXP(-X)*I/sub(ALPHA+K-1)/(X),
!                    K=1,...,N depending on KODE
!           NZ     - number of components of Y set to zero due to
!                    underflow,
!                    NZ=0   , normal return, computation completed
!                    NZ  /=  0, last NZ components of Y set to zero,
!                             Y(K)=0.0D0, K=N-NZ+1,...,N.
!
!     Error Conditions
!         Improper input arguments - a fatal error
!         Overflow with KODE=1 - a fatal error
!         Underflow - a non-fatal error(NZ  /=  0)
!
!***REFERENCES  D. E. Amos, S. L. Daniel and M. K. Weston, CDC 6600
!                 subroutines IBESS and JBESS for Bessel functions
!                 I(NU,X) and J(NU,X), X  >=  0, NU  >=  0, ACM
!                 Transactions on Mathematical Software 3, (1977),
!                 pp. 76-92.
!               F. W. J. Olver, Tables of Bessel Functions of Moderate
!                 or Large Orders, NPL Mathematical Tables 6, Her
!                 Majesty's Stationery Office, London, 1962.
!***ROUTINES CALLED  D1MACH, DASYIK, DLNGAM, I1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBESI
!
  INTEGER I, IALP, IN, INLIM, IS, I1, K, KK, KM, KODE, KT, &
   N, NN, NS, NZ
  INTEGER I1MACH
  DOUBLE PRECISION AIN,AK,AKM,ALPHA,ANS,AP,ARG,ATOL,TOLLN,DFN, &
   DTM, DX, EARG, ELIM, ETX, FLGIK,FN, FNF, FNI,FNP1,FNU,GLN,RA, &
   RTTPI, S, SX, SXO2, S1, S2, T, TA, TB, TEMP, TFN, TM, TOL, &
   TRX, T2, X, XO2, XO2L, Y, Z
  DOUBLE PRECISION D1MACH, DLNGAM
  DIMENSION Y(*), TEMP(3)
  SAVE RTTPI, INLIM
  DATA RTTPI           / 3.98942280401433D-01/
  DATA INLIM           /          80         /
!***FIRST EXECUTABLE STATEMENT  DBESI
  NZ = 0
  KT = 1
!     I1MACH(15) REPLACES I1MACH(12) IN A DOUBLE PRECISION CODE
!     I1MACH(14) REPLACES I1MACH(11) IN A DOUBLE PRECISION CODE
  RA = D1MACH(3)
  TOL = MAX(RA,1.0D-15)
  I1 = -I1MACH(15)
  GLN = D1MACH(5)
  ELIM = 2.303D0*(I1*GLN-3.0D0)
!     TOLLN = -LN(TOL)
  I1 = I1MACH(14)+1
  TOLLN = 2.303D0*GLN*I1
  TOLLN = MIN(TOLLN,34.5388D0)
  if (N-1) 590, 10, 20
   10 KT = 2
   20 NN = N
  if (KODE < 1 .OR. KODE > 2) go to 570
  if (X) 600, 30, 80
   30 if (ALPHA) 580, 40, 50
   40 Y(1) = 1.0D0
  if (N == 1) RETURN
  I1 = 2
  go to 60
   50 I1 = 1
   60 DO 70 I=I1,N
    Y(I) = 0.0D0
   70 CONTINUE
  return
   80 CONTINUE
  if (ALPHA < 0.0D0) go to 580
!
  IALP = INT(ALPHA)
  FNI = IALP + N - 1
  FNF = ALPHA - IALP
  DFN = FNI + FNF
  FNU = DFN
  IN = 0
  XO2 = X*0.5D0
  SXO2 = XO2*XO2
  ETX = KODE - 1
  SX = ETX*X
!
!     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
!     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
!     APPLIED.
!
  if (SXO2 <= (FNU+1.0D0)) go to 90
  if (X <= 12.0D0) go to 110
  FN = 0.55D0*FNU*FNU
  FN = MAX(17.0D0,FN)
  if (X >= FN) go to 430
  ANS = MAX(36.0D0-FNU,0.0D0)
  NS = INT(ANS)
  FNI = FNI + NS
  DFN = FNI + FNF
  FN = DFN
  IS = KT
  KM = N - 1 + NS
  if (KM > 0) IS = 3
  go to 120
   90 FN = FNU
  FNP1 = FN + 1.0D0
  XO2L = LOG(XO2)
  IS = KT
  if (X <= 0.5D0) go to 230
  NS = 0
  100 FNI = FNI + NS
  DFN = FNI + FNF
  FN = DFN
  FNP1 = FN + 1.0D0
  IS = KT
  if (N-1+NS > 0) IS = 3
  go to 230
  110 XO2L = LOG(XO2)
  NS = INT(SXO2-FNU)
  go to 100
  120 CONTINUE
!
!     OVERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION
!
  if (KODE == 2) go to 130
  if (ALPHA < 1.0D0) go to 150
  Z = X/ALPHA
  RA = SQRT(1.0D0+Z*Z)
  GLN = LOG((1.0D0+RA)/Z)
  T = RA*(1.0D0-ETX) + ETX/(Z+RA)
  ARG = ALPHA*(T-GLN)
  if (ARG > ELIM) go to 610
  if (KM == 0) go to 140
  130 CONTINUE
!
!     UNDERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION
!
  Z = X/FN
  RA = SQRT(1.0D0+Z*Z)
  GLN = LOG((1.0D0+RA)/Z)
  T = RA*(1.0D0-ETX) + ETX/(Z+RA)
  ARG = FN*(T-GLN)
  140 if (ARG < (-ELIM)) go to 280
  go to 190
  150 if (X > ELIM) go to 610
  go to 130
!
!     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
!
  160 if (KM /= 0) go to 170
  Y(1) = TEMP(3)
  return
  170 TEMP(1) = TEMP(3)
  IN = NS
  KT = 1
  I1 = 0
  180 CONTINUE
  IS = 2
  FNI = FNI - 1.0D0
  DFN = FNI + FNF
  FN = DFN
  if ( I1 == 2) go to 350
  Z = X/FN
  RA = SQRT(1.0D0+Z*Z)
  GLN = LOG((1.0D0+RA)/Z)
  T = RA*(1.0D0-ETX) + ETX/(Z+RA)
  ARG = FN*(T-GLN)
  190 CONTINUE
  I1 = ABS(3-IS)
  I1 = MAX(I1,1)
  FLGIK = 1.0D0
  call DASYIK(X,FN,KODE,FLGIK,RA,ARG,I1,TEMP(IS))
  go to (180, 350, 510), IS
!
!     SERIES FOR (X/2)**2 <= NU+1
!
  230 CONTINUE
  GLN = DLNGAM(FNP1)
  ARG = FN*XO2L - GLN - SX
  if (ARG < (-ELIM)) go to 300
  EARG = EXP(ARG)
  240 CONTINUE
  S = 1.0D0
  if (X < TOL) go to 260
  AK = 3.0D0
  T2 = 1.0D0
  T = 1.0D0
  S1 = FN
  DO 250 K=1,17
    S2 = T2 + S1
    T = T*SXO2/S2
    S = S + T
    if (ABS(T) < TOL) go to 260
    T2 = T2 + AK
    AK = AK + 2.0D0
    S1 = S1 + FN
  250 CONTINUE
  260 CONTINUE
  TEMP(IS) = S*EARG
  go to (270, 350, 500), IS
  270 EARG = EARG*FN/XO2
  FNI = FNI - 1.0D0
  DFN = FNI + FNF
  FN = DFN
  IS = 2
  go to 240
!
!     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
!
  280 Y(NN) = 0.0D0
  NN = NN - 1
  FNI = FNI - 1.0D0
  DFN = FNI + FNF
  FN = DFN
  if (NN-1) 340, 290, 130
  290 KT = 2
  IS = 2
  go to 130
  300 Y(NN) = 0.0D0
  NN = NN - 1
  FNP1 = FN
  FNI = FNI - 1.0D0
  DFN = FNI + FNF
  FN = DFN
  if (NN-1) 340, 310, 320
  310 KT = 2
  IS = 2
  320 if (SXO2 <= FNP1) go to 330
  go to 130
  330 ARG = ARG - XO2L + LOG(FNP1)
  if (ARG < (-ELIM)) go to 300
  go to 230
  340 NZ = N - NN
  return
!
!     BACKWARD RECURSION SECTION
!
  350 CONTINUE
  NZ = N - NN
  360 CONTINUE
  if ( KT == 2) go to 420
  S1 = TEMP(1)
  S2 = TEMP(2)
  TRX = 2.0D0/X
  DTM = FNI
  TM = (DTM+FNF)*TRX
  if (IN == 0) go to 390
!     BACKWARD RECUR TO INDEX ALPHA+NN-1
  DO 380 I=1,IN
    S = S2
    S2 = TM*S2 + S1
    S1 = S
    DTM = DTM - 1.0D0
    TM = (DTM+FNF)*TRX
  380 CONTINUE
  Y(NN) = S1
  if (NN == 1) RETURN
  Y(NN-1) = S2
  if (NN == 2) RETURN
  go to 400
  390 CONTINUE
!     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
  Y(NN) = S1
  Y(NN-1) = S2
  if (NN == 2) RETURN
  400 K = NN + 1
  DO 410 I=3,NN
    K = K - 1
    Y(K-2) = TM*Y(K-1) + Y(K)
    DTM = DTM - 1.0D0
    TM = (DTM+FNF)*TRX
  410 CONTINUE
  return
  420 Y(1) = TEMP(2)
  return
!
!     ASYMPTOTIC EXPANSION FOR X TO INFINITY
!
  430 CONTINUE
  EARG = RTTPI/SQRT(X)
  if (KODE == 2) go to 440
  if (X > ELIM) go to 610
  EARG = EARG*EXP(X)
  440 ETX = 8.0D0*X
  IS = KT
  IN = 0
  FN = FNU
  450 DX = FNI + FNI
  TM = 0.0D0
  if (FNI == 0.0D0 .AND. ABS(FNF) < TOL) go to 460
  TM = 4.0D0*FNF*(FNI+FNI+FNF)
  460 CONTINUE
  DTM = DX*DX
  S1 = ETX
  TRX = DTM - 1.0D0
  DX = -(TRX+TM)/ETX
  T = DX
  S = 1.0D0 + DX
  ATOL = TOL*ABS(S)
  S2 = 1.0D0
  AK = 8.0D0
  DO 470 K=1,25
    S1 = S1 + ETX
    S2 = S2 + AK
    DX = DTM - S2
    AP = DX + TM
    T = -T*AP/S1
    S = S + T
    if (ABS(T) <= ATOL) go to 480
    AK = AK + 8.0D0
  470 CONTINUE
  480 TEMP(IS) = S*EARG
  if ( IS == 2) go to 360
  IS = 2
  FNI = FNI - 1.0D0
  DFN = FNI + FNF
  FN = DFN
  go to 450
!
!     BACKWARD RECURSION WITH NORMALIZATION BY
!     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.
!
  500 CONTINUE
!     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
  AKM = MAX(3.0D0-FN,0.0D0)
  KM = INT(AKM)
  TFN = FN + KM
  TA = (GLN+TFN-0.9189385332D0-0.0833333333D0/TFN)/(TFN+0.5D0)
  TA = XO2L - TA
  TB = -(1.0D0-1.0D0/TFN)/TFN
  AIN = TOLLN/(-TA+SQRT(TA*TA-TOLLN*TB)) + 1.5D0
  IN = INT(AIN)
  IN = IN + KM
  go to 520
  510 CONTINUE
!     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
  T = 1.0D0/(FN*RA)
  AIN = TOLLN/(GLN+SQRT(GLN*GLN+T*TOLLN)) + 1.5D0
  IN = INT(AIN)
  if (IN > INLIM) go to 160
  520 CONTINUE
  TRX = 2.0D0/X
  DTM = FNI + IN
  TM = (DTM+FNF)*TRX
  TA = 0.0D0
  TB = TOL
  KK = 1
  530 CONTINUE
!
!     BACKWARD RECUR UNINDEXED
!
  DO 540 I=1,IN
    S = TB
    TB = TM*TB + TA
    TA = S
    DTM = DTM - 1.0D0
    TM = (DTM+FNF)*TRX
  540 CONTINUE
!     NORMALIZATION
  if (KK /= 1) go to 550
  TA = (TA/TB)*TEMP(3)
  TB = TEMP(3)
  KK = 2
  IN = NS
  if (NS /= 0) go to 530
  550 Y(NN) = TB
  NZ = N - NN
  if (NN == 1) RETURN
  TB = TM*TB + TA
  K = NN - 1
  Y(K) = TB
  if (NN == 2) RETURN
  DTM = DTM - 1.0D0
  TM = (DTM+FNF)*TRX
  KM = K - 1
!
!     BACKWARD RECUR INDEXED
!
  DO 560 I=1,KM
    Y(K-1) = TM*Y(K) + Y(K+1)
    DTM = DTM - 1.0D0
    TM = (DTM+FNF)*TRX
    K = K - 1
  560 CONTINUE
  return
!
!
!
  570 CONTINUE
  call XERMSG ('SLATEC', 'DBESI', &
     'SCALING OPTION, KODE, NOT 1 OR 2.', 2, 1)
  return
  580 CONTINUE
  call XERMSG ('SLATEC', 'DBESI', 'ORDER, ALPHA, LESS THAN ZERO.', &
     2, 1)
  return
  590 CONTINUE
  call XERMSG ('SLATEC', 'DBESI', 'N LESS THAN ONE.', 2, 1)
  return
  600 CONTINUE
  call XERMSG ('SLATEC', 'DBESI', 'X LESS THAN ZERO.', 2, 1)
  return
  610 CONTINUE
  call XERMSG ('SLATEC', 'DBESI', &
     'OVERFLOW, X TOO LARGE FOR KODE = 1.', 6, 1)
  return
end
