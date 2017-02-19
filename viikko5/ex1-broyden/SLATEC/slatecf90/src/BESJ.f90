subroutine BESJ (X, ALPHA, N, Y, NZ)
!
!! BESJ computes an N member sequence of J Bessel functions ...
!  J/SUB(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA and X.
!
!***LIBRARY   SLATEC
!***CATEGORY  C10A3
!***TYPE      SINGLE PRECISION (BESJ-S, DBESJ-D)
!***KEYWORDS  J BESSEL FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Amos, D. E., (SNLA)
!           Daniel, S. L., (SNLA)
!           Weston, M. K., (SNLA)
!***DESCRIPTION
!
!     Abstract
!         BESJ computes an N member sequence of J Bessel functions
!         J/sub(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA and X.
!         A combination of the power series, the asymptotic expansion
!         for X to infinity and the uniform asymptotic expansion for
!         NU to infinity are applied over subdivisions of the (NU,X)
!         plane.  For values of (NU,X) not covered by one of these
!         formulae, the order is incremented or decremented by integer
!         values into a region where one of the formulae apply. Backward
!         recursion is applied to reduce orders by integer values except
!         where the entire sequence lies in the oscillatory region.  In
!         this case forward recursion is stable and values from the
!         asymptotic expansion for X to infinity start the recursion
!         when it is efficient to do so.  Leading terms of the series
!         and uniform expansion are tested for underflow.  If a sequence
!         is requested and the last member would underflow, the result
!         is set to zero and the next lower order tried, etc., until a
!         member comes on scale or all members are set to zero.
!         Overflow cannot occur.
!
!     Description of Arguments
!
!         Input
!           X      - X  >=  0.0E0
!           ALPHA  - order of first member of the sequence,
!                    ALPHA  >=  0.0E0
!           N      - number of members in the sequence, N  >=  1
!
!         Output
!           Y      - a vector whose first  N components contain
!                    values for J/sub(ALPHA+K-1)/(X), K=1,...,N
!           NZ     - number of components of Y set to zero due to
!                    underflow,
!                    NZ=0   , normal return, computation completed
!                    NZ  /=  0, last NZ components of Y set to zero,
!                             Y(K)=0.0E0, K=N-NZ+1,...,N.
!
!     Error Conditions
!         Improper input arguments - a fatal error
!         Underflow  - a non-fatal error (NZ  /=  0)
!
!***REFERENCES  D. E. Amos, S. L. Daniel and M. K. Weston, CDC 6600
!                 subroutines IBESS and JBESS for Bessel functions
!                 I(NU,X) and J(NU,X), X  >=  0, NU  >=  0, ACM
!                 Transactions on Mathematical Software 3, (1977),
!                 pp. 76-92.
!               F. W. J. Olver, Tables of Bessel Functions of Moderate
!                 or Large Orders, NPL Mathematical Tables 6, Her
!                 Majesty's Stationery Office, London, 1962.
!***ROUTINES CALLED  ALNGAM, ASYJY, I1MACH, JAIRY, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BESJ
  EXTERNAL JAIRY
  INTEGER I,IALP,IDALP,IFLW,IN,INLIM,IS,I1,I2,K,KK,KM,KT,N,NN, &
          NS,NZ
  INTEGER I1MACH
  REAL       AK,AKM,ALPHA,ANS,AP,ARG,COEF,DALPHA,DFN,DTM,EARG, &
             ELIM1,ETX,FIDAL,FLGJY,FN,FNF,FNI,FNP1,FNU,FNULIM, &
             GLN,PDF,PIDT,PP,RDEN,RELB,RTTP,RTWO,RTX,RZDEN, &
             S,SA,SB,SXO2,S1,S2,T,TA,TAU,TB,TEMP,TFN,TM,TOL, &
             TOLLN,TRX,TX,T1,T2,WK,X,XO2,XO2L,Y,RTOL,SLIM
  SAVE RTWO, PDF, RTTP, PIDT, PP, INLIM, FNULIM
  REAL R1MACH, ALNGAM
  DIMENSION Y(*), TEMP(3), FNULIM(2), PP(4), WK(7)
  DATA RTWO,PDF,RTTP,PIDT                    / 1.34839972492648E+00, &
   7.85398163397448E-01, 7.97884560802865E-01, 1.57079632679490E+00/
  DATA  PP(1),  PP(2),  PP(3),  PP(4)        / 8.72909153935547E+00, &
   2.65693932265030E-01, 1.24578576865586E-01, 7.70133747430388E-04/
  DATA INLIM           /      150            /
  DATA FNULIM(1), FNULIM(2) /      100.0E0,     60.0E0     /
!***FIRST EXECUTABLE STATEMENT  BESJ
  NZ = 0
  KT = 1
  NS=0
!     I1MACH(14) REPLACES I1MACH(11) IN A DOUBLE PRECISION CODE
!     I1MACH(15) REPLACES I1MACH(12) IN A DOUBLE PRECISION CODE
  TA = R1MACH(3)
  TOL = MAX(TA,1.0E-15)
  I1 = I1MACH(11) + 1
  I2 = I1MACH(12)
  TB = R1MACH(5)
  ELIM1 = -2.303E0*(I2*TB+3.0E0)
  RTOL=1.0E0/TOL
  SLIM=R1MACH(1)*1.0E+3*RTOL
!     TOLLN = -LN(TOL)
  TOLLN = 2.303E0*TB*I1
  TOLLN = MIN(TOLLN,34.5388E0)
  if (N-1) 720, 10, 20
   10 KT = 2
   20 NN = N
  if (X) 730, 30, 80
   30 if (ALPHA) 710, 40, 50
   40 Y(1) = 1.0E0
  if (N == 1) RETURN
  I1 = 2
  go to 60
   50 I1 = 1
   60 DO 70 I=I1,N
    Y(I) = 0.0E0
   70 CONTINUE
  return
   80 CONTINUE
  if (ALPHA < 0.0E0) go to 710
!
  IALP = INT(ALPHA)
  FNI = IALP + N - 1
  FNF = ALPHA - IALP
  DFN = FNI + FNF
  FNU = DFN
  XO2 = X*0.5E0
  SXO2 = XO2*XO2
!
!     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
!     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
!     APPLIED.
!
  if (SXO2 <= (FNU+1.0E0)) go to 90
  TA = MAX(20.0E0,FNU)
  if (X > TA) go to 120
  if (X > 12.0E0) go to 110
  XO2L = LOG(XO2)
  NS = INT(SXO2-FNU) + 1
  go to 100
   90 FN = FNU
  FNP1 = FN + 1.0E0
  XO2L = LOG(XO2)
  IS = KT
  if (X <= 0.50E0) go to 330
  NS = 0
  100 FNI = FNI + NS
  DFN = FNI + FNF
  FN = DFN
  FNP1 = FN + 1.0E0
  IS = KT
  if (N-1+NS > 0) IS = 3
  go to 330
  110 ANS = MAX(36.0E0-FNU,0.0E0)
  NS = INT(ANS)
  FNI = FNI + NS
  DFN = FNI + FNF
  FN = DFN
  IS = KT
  if (N-1+NS > 0) IS = 3
  go to 130
  120 CONTINUE
  RTX = SQRT(X)
  TAU = RTWO*RTX
  TA = TAU + FNULIM(KT)
  if (FNU <= TA) go to 480
  FN = FNU
  IS = KT
!
!     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
!
  130 CONTINUE
  I1 = ABS(3-IS)
  I1 = MAX(I1,1)
  FLGJY = 1.0E0
  call ASYJY(JAIRY,X,FN,FLGJY,I1,TEMP(IS),WK,IFLW)
  if ( IFLW /= 0) go to 380
  go to (320, 450, 620), IS
  310 TEMP(1) = TEMP(3)
  KT = 1
  320 IS = 2
  FNI = FNI - 1.0E0
  DFN = FNI + FNF
  FN = DFN
  if ( I1 == 2) go to 450
  go to 130
!
!     SERIES FOR (X/2)**2 <= NU+1
!
  330 CONTINUE
  GLN = ALNGAM(FNP1)
  ARG = FN*XO2L - GLN
  if (ARG < (-ELIM1)) go to 400
  EARG = EXP(ARG)
  340 CONTINUE
  S = 1.0E0
  if (X < TOL) go to 360
  AK = 3.0E0
  T2 = 1.0E0
  T = 1.0E0
  S1 = FN
  DO 350 K=1,17
    S2 = T2 + S1
    T = -T*SXO2/S2
    S = S + T
    if (ABS(T) < TOL) go to 360
    T2 = T2 + AK
    AK = AK + 2.0E0
    S1 = S1 + FN
  350 CONTINUE
  360 CONTINUE
  TEMP(IS) = S*EARG
  go to (370, 450, 610), IS
  370 EARG = EARG*FN/XO2
  FNI = FNI - 1.0E0
  DFN = FNI + FNF
  FN = DFN
  IS = 2
  go to 340
!
!     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
!     UNDERFLOW CAN ONLY OCCUR FOR NS=0 SINCE THE ORDER MUST BE
!     LARGER THAN 36. THEREFORE, NS NEED NOT BE CONSIDERED.
!
  380 Y(NN) = 0.0E0
  NN = NN - 1
  FNI = FNI - 1.0E0
  DFN = FNI + FNF
  FN = DFN
  if (NN-1) 440, 390, 130
  390 KT = 2
  IS = 2
  go to 130
  400 Y(NN) = 0.0E0
  NN = NN - 1
  FNP1 = FN
  FNI = FNI - 1.0E0
  DFN = FNI + FNF
  FN = DFN
  if (NN-1) 440, 410, 420
  410 KT = 2
  IS = 2
  420 if (SXO2 <= FNP1) go to 430
  go to 130
  430 ARG = ARG - XO2L + LOG(FNP1)
  if (ARG < (-ELIM1)) go to 400
  go to 330
  440 NZ = N - NN
  return
!
!     BACKWARD RECURSION SECTION
!
  450 CONTINUE
  if ( NS /= 0) go to 451
  NZ = N - NN
  if (KT == 2) go to 470
!     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
  Y(NN) = TEMP(1)
  Y(NN-1) = TEMP(2)
  if (NN == 2) RETURN
  451 CONTINUE
  TRX = 2.0E0/X
  DTM = FNI
  TM = (DTM+FNF)*TRX
  AK=1.0E0
  TA=TEMP(1)
  TB=TEMP(2)
  if ( ABS(TA) > SLIM) go to 455
  TA=TA*RTOL
  TB=TB*RTOL
  AK=TOL
  455 CONTINUE
  KK=2
  IN=NS-1
  if ( IN == 0) go to 690
  if ( NS /= 0) go to 670
  K=NN-2
  DO 460 I=3,NN
    S=TB
    TB=TM*TB-TA
    TA=S
    Y(K)=TB*AK
    K=K-1
    DTM = DTM - 1.0E0
    TM = (DTM+FNF)*TRX
  460 CONTINUE
  return
  470 Y(1) = TEMP(2)
  return
!
!     ASYMPTOTIC EXPANSION FOR X TO INFINITY WITH FORWARD RECURSION IN
!     OSCILLATORY REGION X > MAX(20, NU), PROVIDED THE LAST MEMBER
!     OF THE SEQUENCE IS ALSO IN THE REGION.
!
  480 CONTINUE
  IN = INT(ALPHA-TAU+2.0E0)
  if (IN <= 0) go to 490
  IDALP = IALP - IN - 1
  KT = 1
  go to 500
  490 CONTINUE
  IDALP = IALP
  IN = 0
  500 IS = KT
  FIDAL = IDALP
  DALPHA = FIDAL + FNF
  ARG = X - PIDT*DALPHA - PDF
  SA = SIN(ARG)
  SB = COS(ARG)
  COEF = RTTP/RTX
  ETX = 8.0E0*X
  510 CONTINUE
  DTM = FIDAL + FIDAL
  DTM = DTM*DTM
  TM = 0.0E0
  if (FIDAL == 0.0E0 .AND. ABS(FNF) < TOL) go to 520
  TM = 4.0E0*FNF*(FIDAL+FIDAL+FNF)
  520 CONTINUE
  TRX = DTM - 1.0E0
  T2 = (TRX+TM)/ETX
  S2 = T2
  RELB = TOL*ABS(T2)
  T1 = ETX
  S1 = 1.0E0
  FN = 1.0E0
  AK = 8.0E0
  DO 530 K=1,13
    T1 = T1 + ETX
    FN = FN + AK
    TRX = DTM - FN
    AP = TRX + TM
    T2 = -T2*AP/T1
    S1 = S1 + T2
    T1 = T1 + ETX
    AK = AK + 8.0E0
    FN = FN + AK
    TRX = DTM - FN
    AP = TRX + TM
    T2 = T2*AP/T1
    S2 = S2 + T2
    if (ABS(T2) <= RELB) go to 540
    AK = AK + 8.0E0
  530 CONTINUE
  540 TEMP(IS) = COEF*(S1*SB-S2*SA)
  if ( IS == 2) go to 560
  FIDAL = FIDAL + 1.0E0
  DALPHA = FIDAL + FNF
  IS = 2
  TB = SA
  SA = -SB
  SB = TB
  go to 510
!
!     FORWARD RECURSION SECTION
!
  560 if (KT == 2) go to 470
  S1 = TEMP(1)
  S2 = TEMP(2)
  TX = 2.0E0/X
  TM = DALPHA*TX
  if (IN == 0) go to 580
!
!     FORWARD RECUR TO INDEX ALPHA
!
  DO 570 I=1,IN
    S = S2
    S2 = TM*S2 - S1
    TM = TM + TX
    S1 = S
  570 CONTINUE
  if (NN == 1) go to 600
  S = S2
  S2 = TM*S2 - S1
  TM = TM + TX
  S1 = S
  580 CONTINUE
!
!     FORWARD RECUR FROM INDEX ALPHA TO ALPHA+N-1
!
  Y(1) = S1
  Y(2) = S2
  if (NN == 2) RETURN
  DO 590 I=3,NN
    Y(I) = TM*Y(I-1) - Y(I-2)
    TM = TM + TX
  590 CONTINUE
  return
  600 Y(1) = S2
  return
!
!     BACKWARD RECURSION WITH NORMALIZATION BY
!     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.
!
  610 CONTINUE
!     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
  AKM = MAX(3.0E0-FN,0.0E0)
  KM = INT(AKM)
  TFN = FN + KM
  TA = (GLN+TFN-0.9189385332E0-0.0833333333E0/TFN)/(TFN+0.5E0)
  TA = XO2L - TA
  TB = -(1.0E0-1.5E0/TFN)/TFN
  AKM = TOLLN/(-TA+SQRT(TA*TA-TOLLN*TB)) + 1.5E0
  IN = KM + INT(AKM)
  go to 660
  620 CONTINUE
!     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
  GLN = WK(3) + WK(2)
  if (WK(6) > 30.0E0) go to 640
  RDEN = (PP(4)*WK(6)+PP(3))*WK(6) + 1.0E0
  RZDEN = PP(1) + PP(2)*WK(6)
  TA = RZDEN/RDEN
  if (WK(1) < 0.10E0) go to 630
  TB = GLN/WK(5)
  go to 650
  630 TB=(1.259921049E0+(0.1679894730E0+0.0887944358E0*WK(1))*WK(1)) &
   /WK(7)
  go to 650
  640 CONTINUE
  TA = 0.5E0*TOLLN/WK(4)
  TA=((0.0493827160E0*TA-0.1111111111E0)*TA+0.6666666667E0)*TA*WK(6)
  if (WK(1) < 0.10E0) go to 630
  TB = GLN/WK(5)
  650 IN = INT(TA/TB+1.5E0)
  if (IN > INLIM) go to 310
  660 CONTINUE
  DTM = FNI + IN
  TRX = 2.0E0/X
  TM = (DTM+FNF)*TRX
  TA = 0.0E0
  TB = TOL
  KK = 1
  AK=1.0E0
  670 CONTINUE
!
!     BACKWARD RECUR UNINDEXED AND SCALE WHEN MAGNITUDES ARE CLOSE TO
!     UNDERFLOW LIMITS (LESS THAN SLIM=R1MACH(1)*1.0E+3/TOL)
!
  DO 680 I=1,IN
    S = TB
    TB = TM*TB - TA
    TA = S
    DTM = DTM - 1.0E0
    TM = (DTM+FNF)*TRX
  680 CONTINUE
!     NORMALIZATION
  if (KK /= 1) go to 690
  S=TEMP(3)
  SA=TA/TB
  TA=S
  TB=S
  if ( ABS(S) > SLIM) go to 685
  TA=TA*RTOL
  TB=TB*RTOL
  AK=TOL
  685 CONTINUE
  TA=TA*SA
  KK = 2
  IN = NS
  if (NS /= 0) go to 670
  690 Y(NN) = TB*AK
  NZ = N - NN
  if (NN == 1) RETURN
  K = NN - 1
  S=TB
  TB = TM*TB - TA
  TA=S
  Y(K)=TB*AK
  if (NN == 2) RETURN
  DTM = DTM - 1.0E0
  TM = (DTM+FNF)*TRX
  K=NN-2
!
!     BACKWARD RECUR INDEXED
!
  DO 700 I=3,NN
    S=TB
    TB = TM*TB - TA
    TA=S
    Y(K)=TB*AK
    DTM = DTM - 1.0E0
    TM = (DTM+FNF)*TRX
    K = K - 1
  700 CONTINUE
  return
!
!
!
  710 CONTINUE
  call XERMSG ('SLATEC', 'BESJ', 'ORDER, ALPHA, LESS THAN ZERO.', &
     2, 1)
  return
  720 CONTINUE
  call XERMSG ('SLATEC', 'BESJ', 'N LESS THAN ONE.', 2, 1)
  return
  730 CONTINUE
  call XERMSG ('SLATEC', 'BESJ', 'X LESS THAN ZERO.', 2, 1)
  return
end
