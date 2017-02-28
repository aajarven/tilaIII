subroutine BKIAS (X, N, KTRMS, T, ANS, IND, MS, GMRN, H, IERR)
!
!! BKIAS is subsidiary to BSKIN.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (BKIAS-S, DBKIAS-D)
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     BKIAS computes repeated integrals of the K0 Bessel function
!     by the asymptotic expansion
!
!***SEE ALSO  BSKIN
!***ROUTINES CALLED  BDIFF, GAMRN, HKSEQ, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   820601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  BKIAS
  INTEGER I, II, IND, J, JMI, JN, K, KK, KM, KTRMS, MM, MP, MS, N, &
   IERR
  REAL ANS, B, BND, DEN1, DEN2, DEN3, ER, ERR, FJ, FK, FLN, FM1, &
   GMRN, G1, GS, H, HN, HRTPI, RAT, RG1, RXP, RZ, RZX, S, SS, SUMI, &
   SUMJ, T, TOL, V, W, X, XP, Z
  REAL GAMRN, R1MACH
  DIMENSION B(120), XP(16), S(31), H(*), V(52), W(52), T(50), &
   BND(15)
  SAVE B, BND, HRTPI
!-----------------------------------------------------------------------
!             COEFFICIENTS OF POLYNOMIAL P(J-1,X), J=1,15
!-----------------------------------------------------------------------
  DATA B(1), B(2), B(3), B(4), B(5), B(6), B(7), B(8), B(9), B(10), &
   B(11), B(12), B(13), B(14), B(15), B(16), B(17), B(18), B(19), &
   B(20), B(21), B(22), B(23), B(24) /1.00000000000000000E+00, &
   1.00000000000000000E+00,-2.00000000000000000E+00, &
   1.00000000000000000E+00,-8.00000000000000000E+00, &
   6.00000000000000000E+00,1.00000000000000000E+00, &
   -2.20000000000000000E+01,5.80000000000000000E+01, &
   -2.40000000000000000E+01,1.00000000000000000E+00, &
   -5.20000000000000000E+01,3.28000000000000000E+02, &
   -4.44000000000000000E+02,1.20000000000000000E+02, &
   1.00000000000000000E+00,-1.14000000000000000E+02, &
   1.45200000000000000E+03,-4.40000000000000000E+03, &
   3.70800000000000000E+03,-7.20000000000000000E+02, &
   1.00000000000000000E+00,-2.40000000000000000E+02, &
   5.61000000000000000E+03/
  DATA B(25), B(26), B(27), B(28), B(29), B(30), B(31), B(32), &
   B(33), B(34), B(35), B(36), B(37), B(38), B(39), B(40), B(41), &
   B(42), B(43), B(44), B(45), B(46), B(47), B(48) &
   /-3.21200000000000000E+04,5.81400000000000000E+04, &
   -3.39840000000000000E+04,5.04000000000000000E+03, &
   1.00000000000000000E+00,-4.94000000000000000E+02, &
   1.99500000000000000E+04,-1.95800000000000000E+05, &
   6.44020000000000000E+05,-7.85304000000000000E+05, &
   3.41136000000000000E+05,-4.03200000000000000E+04, &
   1.00000000000000000E+00,-1.00400000000000000E+03, &
   6.72600000000000000E+04,-1.06250000000000000E+06, &
   5.76550000000000000E+06,-1.24400640000000000E+07, &
   1.10262960000000000E+07,-3.73392000000000000E+06, &
   3.62880000000000000E+05,1.00000000000000000E+00, &
   -2.02600000000000000E+03,2.18848000000000000E+05/
  DATA B(49), B(50), B(51), B(52), B(53), B(54), B(55), B(56), &
   B(57), B(58), B(59), B(60), B(61), B(62), B(63), B(64), B(65), &
   B(66), B(67), B(68), B(69), B(70), B(71), B(72) &
   /-5.32616000000000000E+06,4.47650000000000000E+07, &
   -1.55357384000000000E+08,2.38904904000000000E+08, &
   -1.62186912000000000E+08,4.43390400000000000E+07, &
   -3.62880000000000000E+06,1.00000000000000000E+00, &
   -4.07200000000000000E+03,6.95038000000000000E+05, &
   -2.52439040000000000E+07,3.14369720000000000E+08, &
   -1.64838430400000000E+09,4.00269508800000000E+09, &
   -4.64216395200000000E+09,2.50748121600000000E+09, &
   -5.68356480000000000E+08,3.99168000000000000E+07, &
   1.00000000000000000E+00,-8.16600000000000000E+03, &
   2.17062600000000000E+06,-1.14876376000000000E+08, &
   2.05148277600000000E+09,-1.55489607840000000E+10/
  DATA B(73), B(74), B(75), B(76), B(77), B(78), B(79), B(80), &
   B(81), B(82), B(83), B(84), B(85), B(86), B(87), B(88), B(89), &
   B(90), B(91), B(92), B(93), B(94), B(95), B(96) &
   /5.60413987840000000E+10,-1.01180433024000000E+11, &
   9.21997902240000000E+10,-4.07883018240000000E+10, &
   7.82771904000000000E+09,-4.79001600000000000E+08, &
   1.00000000000000000E+00,-1.63560000000000000E+04, &
   6.69969600000000000E+06,-5.07259276000000000E+08, &
   1.26698177760000000E+10,-1.34323420224000000E+11, &
   6.87720046384000000E+11,-1.81818864230400000E+12, &
   2.54986547342400000E+12,-1.88307966182400000E+12, &
   6.97929436800000000E+11,-1.15336085760000000E+11, &
   6.22702080000000000E+09,1.00000000000000000E+00, &
   -3.27380000000000000E+04,2.05079880000000000E+07, &
   -2.18982980800000000E+09,7.50160522280000000E+10/
  DATA B(97), B(98), B(99), B(100), B(101), B(102), B(103), B(104), &
   B(105), B(106), B(107), B(108), B(109), B(110), B(111), B(112), &
   B(113), B(114), B(115), B(116), B(117), B(118) &
   /-1.08467651241600000E+12,7.63483214939200000E+12, &
   -2.82999100661120000E+13,5.74943734645920000E+13, &
   -6.47283751398720000E+13,3.96895780558080000E+13, &
   -1.25509040179200000E+13,1.81099255680000000E+12, &
   -8.71782912000000000E+10,1.00000000000000000E+00, &
   -6.55040000000000000E+04,6.24078900000000000E+07, &
   -9.29252692000000000E+09,4.29826006340000000E+11, &
   -8.30844432796800000E+12,7.83913848313120000E+13, &
   -3.94365587815520000E+14,1.11174747256968000E+15, &
   -1.79717122069056000E+15,1.66642448627145600E+15, &
   -8.65023253219584000E+14,2.36908271543040000E+14/
  DATA B(119), B(120) /-3.01963769856000000E+13, &
   1.30767436800000000E+12/
!-----------------------------------------------------------------------
!             BOUNDS B(M,K) , K=M-3
!-----------------------------------------------------------------------
  DATA BND(1), BND(2), BND(3), BND(4), BND(5), BND(6), BND(7), &
   BND(8), BND(9), BND(10), BND(11), BND(12), BND(13), BND(14), &
   BND(15) /1.0E0,1.0E0,1.0E0,1.0E0,3.10E0,5.18E0,11.7E0,29.8E0, &
   90.4E0,297.0E0,1070.0E0,4290.0E0,18100.0E0,84700.0E0,408000.0E0/
  DATA HRTPI /8.86226925452758014E-01/
!
!***FIRST EXECUTABLE STATEMENT  BKIAS
  IERR=0
  TOL = MAX(R1MACH(4),1.0E-18)
  FLN = N
  RZ = 1.0E0/(X+FLN)
  RZX = X*RZ
  Z = 0.5E0*(X+FLN)
  if (IND > 1) go to 10
  GMRN = GAMRN(Z)
   10 CONTINUE
  GS = HRTPI*GMRN
  G1 = GS + GS
  RG1 = 1.0E0/G1
  GMRN = (RZ+RZ)/GMRN
  if (IND > 1) go to 70
!-----------------------------------------------------------------------
!     EVALUATE ERROR FOR M=MS
!-----------------------------------------------------------------------
  HN = 0.5E0*FLN
  DEN2 = KTRMS + KTRMS + N
  DEN3 = DEN2 - 2.0E0
  DEN1 = X + DEN2
  ERR = RG1*(X+X)/(DEN1-1.0E0)
  if (N == 0) go to 20
  RAT = 1.0E0/(FLN*FLN)
   20 CONTINUE
  if (KTRMS == 0) go to 30
  FJ = KTRMS
  RAT = 0.25E0/(HRTPI*DEN3*SQRT(FJ))
   30 CONTINUE
  ERR = ERR*RAT
  FJ = -3.0E0
  DO 50 J=1,15
    if (J <= 5) ERR = ERR/DEN1
    FM1 = MAX(1.0E0,FJ)
    FJ = FJ + 1.0E0
    ER = BND(J)*ERR
    if (KTRMS == 0) go to 40
    ER = ER/FM1
    if (ER < TOL) go to 60
    if (J >= 5) ERR = ERR/DEN3
    go to 50
   40   CONTINUE
    ER = ER*(1.0E0+HN/FM1)
    if (ER < TOL) go to 60
    if (J >= 5) ERR = ERR/FLN
   50 CONTINUE
  go to 200
   60 CONTINUE
  MS = J
   70 CONTINUE
  MM = MS + MS
  MP = MM + 1
!-----------------------------------------------------------------------
!     H(K)=(-Z)**(K)*(PSI(K-1,Z)-PSI(K-1,Z+0.5))/GAMMA(K) , K=1,2,...,MM
!-----------------------------------------------------------------------
  if (IND > 1) go to 80
  call HKSEQ(Z, MM, H, IERR)
  go to 100
   80 CONTINUE
  RAT = Z/(Z-0.5E0)
  RXP = RAT
  DO 90 I=1,MM
    H(I) = RXP*(1.0E0-H(I))
    RXP = RXP*RAT
   90 CONTINUE
  100 CONTINUE
!-----------------------------------------------------------------------
!     SCALED S SEQUENCE
!-----------------------------------------------------------------------
  S(1) = 1.0E0
  FK = 1.0E0
  DO 120 K=2,MP
    SS = 0.0E0
    KM = K - 1
    I = KM
    DO 110 J=1,KM
      SS = SS + S(J)*H(I)
      I = I - 1
  110   CONTINUE
    S(K) = SS/FK
    FK = FK + 1.0E0
  120 CONTINUE
!-----------------------------------------------------------------------
!     SCALED S-TILDA SEQUENCE
!-----------------------------------------------------------------------
  if (KTRMS == 0) go to 160
  FK = 0.0E0
  SS = 0.0E0
  RG1 = RG1/Z
  DO 130 K=1,KTRMS
    V(K) = Z/(Z+FK)
    W(K) = T(K)*V(K)
    SS = SS + W(K)
    FK = FK + 1.0E0
  130 CONTINUE
  S(1) = S(1) - SS*RG1
  DO 150 I=2,MP
    SS = 0.0E0
    DO 140 K=1,KTRMS
      W(K) = W(K)*V(K)
      SS = SS + W(K)
  140   CONTINUE
    S(I) = S(I) - SS*RG1
  150 CONTINUE
  160 CONTINUE
!-----------------------------------------------------------------------
!     SUM ON J
!-----------------------------------------------------------------------
  SUMJ = 0.0E0
  JN = 1
  RXP = 1.0E0
  XP(1) = 1.0E0
  DO 190 J=1,MS
    JN = JN + J - 1
    XP(J+1) = XP(J)*RZX
    RXP = RXP*RZ
!-----------------------------------------------------------------------
!     SUM ON I
!-----------------------------------------------------------------------
    SUMI = 0.0E0
    II = JN
    DO 180 I=1,J
      JMI = J - I + 1
      KK = J + I + 1
      DO 170 K=1,JMI
        V(K) = S(KK)*XP(K)
        KK = KK + 1
  170     CONTINUE
      call BDIFF(JMI, V)
      SUMI = SUMI + B(II)*V(JMI)*XP(I+1)
      II = II + 1
  180   CONTINUE
    SUMJ = SUMJ + SUMI*RXP
  190 CONTINUE
  ANS = GS*(S(1)-SUMJ)
  return
  200 CONTINUE
  IERR=2
  return
end