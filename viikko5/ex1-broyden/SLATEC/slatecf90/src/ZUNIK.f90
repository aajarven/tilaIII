subroutine ZUNIK (ZRR, ZRI, FNU, IKFLG, IPMTR, TOL, INIT, PHIR, &
     PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
!
!! ZUNIK is subsidiary to ZBESI and ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNIK-A, ZUNIK-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!        ZUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
!        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
!        RESPECTIVELY BY
!
!        W(FNU,ZR) = PHI*EXP(ZETA)*SUM
!
!        WHERE       ZETA=-ZETA1 + ZETA2       OR
!                          ZETA1 - ZETA2
!
!        THE FIRST call MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
!        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG=
!        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
!        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
!        ZETA1,ZETA2.
!
!***SEE ALSO  ZBESI, ZBESK
!***ROUTINES CALLED  D1MACH, ZDIV, ZLOG, ZSQRT
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   930122  Added EXTERNAL statement with ZLOG and ZSQRT.  (RWC)
!***END PROLOGUE  ZUNIK
!     COMPLEX CFN,CON,CONE,CRFN,CWRK,CZERO,PHI,S,SR,SUM,T,T2,ZETA1,
!    *ZETA2,ZN,ZR
  DOUBLE PRECISION AC, C, CON, CONEI, CONER, CRFNI, CRFNR, CWRKI, &
   CWRKR, FNU, PHII, PHIR, RFN, SI, SR, SRI, SRR, STI, STR, SUMI, &
   SUMR, TEST, TI, TOL, TR, T2I, T2R, ZEROI, ZEROR, ZETA1I, ZETA1R, &
   ZETA2I, ZETA2R, ZNI, ZNR, ZRI, ZRR, D1MACH
  INTEGER I, IDUM, IKFLG, INIT, IPMTR, J, K, L
  DIMENSION C(120), CWRKR(16), CWRKI(16), CON(2)
  EXTERNAL ZLOG, ZSQRT
  DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
  DATA CON(1), CON(2)  / &
   3.98942280401432678D-01,  1.25331413731550025D+00 /
  DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10), &
       C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18), &
       C(19), C(20), C(21), C(22), C(23), C(24)/ &
       1.00000000000000000D+00,    -2.08333333333333333D-01, &
       1.25000000000000000D-01,     3.34201388888888889D-01, &
      -4.01041666666666667D-01,     7.03125000000000000D-02, &
      -1.02581259645061728D+00,     1.84646267361111111D+00, &
      -8.91210937500000000D-01,     7.32421875000000000D-02, &
       4.66958442342624743D+00,    -1.12070026162229938D+01, &
       8.78912353515625000D+00,    -2.36408691406250000D+00, &
       1.12152099609375000D-01,    -2.82120725582002449D+01, &
       8.46362176746007346D+01,    -9.18182415432400174D+01, &
       4.25349987453884549D+01,    -7.36879435947963170D+00, &
       2.27108001708984375D-01,     2.12570130039217123D+02, &
      -7.65252468141181642D+02,     1.05999045252799988D+03/
  DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32), &
       C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40), &
       C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/ &
      -6.99579627376132541D+02,     2.18190511744211590D+02, &
      -2.64914304869515555D+01,     5.72501420974731445D-01, &
      -1.91945766231840700D+03,     8.06172218173730938D+03, &
      -1.35865500064341374D+04,     1.16553933368645332D+04, &
      -5.30564697861340311D+03,     1.20090291321635246D+03, &
      -1.08090919788394656D+02,     1.72772750258445740D+00, &
       2.02042913309661486D+04,    -9.69805983886375135D+04, &
       1.92547001232531532D+05,    -2.03400177280415534D+05, &
       1.22200464983017460D+05,    -4.11926549688975513D+04, &
       7.10951430248936372D+03,    -4.93915304773088012D+02, &
       6.07404200127348304D+00,    -2.42919187900551333D+05, &
       1.31176361466297720D+06,    -2.99801591853810675D+06/
  DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56), &
       C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64), &
       C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/ &
       3.76327129765640400D+06,    -2.81356322658653411D+06, &
       1.26836527332162478D+06,    -3.31645172484563578D+05, &
       4.52187689813627263D+04,    -2.49983048181120962D+03, &
       2.43805296995560639D+01,     3.28446985307203782D+06, &
      -1.97068191184322269D+07,     5.09526024926646422D+07, &
      -7.41051482115326577D+07,     6.63445122747290267D+07, &
      -3.75671766607633513D+07,     1.32887671664218183D+07, &
      -2.78561812808645469D+06,     3.08186404612662398D+05, &
      -1.38860897537170405D+04,     1.10017140269246738D+02, &
      -4.93292536645099620D+07,     3.25573074185765749D+08, &
      -9.39462359681578403D+08,     1.55359689957058006D+09, &
      -1.62108055210833708D+09,     1.10684281682301447D+09/
  DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80), &
       C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88), &
       C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/ &
      -4.95889784275030309D+08,     1.42062907797533095D+08, &
      -2.44740627257387285D+07,     2.24376817792244943D+06, &
      -8.40054336030240853D+04,     5.51335896122020586D+02, &
       8.14789096118312115D+08,    -5.86648149205184723D+09, &
       1.86882075092958249D+10,    -3.46320433881587779D+10, &
       4.12801855797539740D+10,    -3.30265997498007231D+10, &
       1.79542137311556001D+10,    -6.56329379261928433D+09, &
       1.55927986487925751D+09,    -2.25105661889415278D+08, &
       1.73951075539781645D+07,    -5.49842327572288687D+05, &
       3.03809051092238427D+03,    -1.46792612476956167D+10, &
       1.14498237732025810D+11,    -3.99096175224466498D+11, &
       8.19218669548577329D+11,    -1.09837515608122331D+12/
  DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104), &
       C(105), C(106), C(107), C(108), C(109), C(110), C(111), &
       C(112), C(113), C(114), C(115), C(116), C(117), C(118)/ &
       1.00815810686538209D+12,    -6.45364869245376503D+11, &
       2.87900649906150589D+11,    -8.78670721780232657D+10, &
       1.76347306068349694D+10,    -2.16716498322379509D+09, &
       1.43157876718888981D+08,    -3.87183344257261262D+06, &
       1.82577554742931747D+04,     2.86464035717679043D+11, &
      -2.40629790002850396D+12,     9.10934118523989896D+12, &
      -2.05168994109344374D+13,     3.05651255199353206D+13, &
      -3.16670885847851584D+13,     2.33483640445818409D+13, &
      -1.23204913055982872D+13,     4.61272578084913197D+12, &
      -1.19655288019618160D+12,     2.05914503232410016D+11, &
      -2.18229277575292237D+10,     1.24700929351271032D+09/
  DATA C(119), C(120)/ &
      -2.91883881222208134D+07,     1.18838426256783253D+05/
!***FIRST EXECUTABLE STATEMENT  ZUNIK
  if (INIT /= 0) go to 40
!-----------------------------------------------------------------------
!     INITIALIZE ALL VARIABLES
!-----------------------------------------------------------------------
  RFN = 1.0D0/FNU
!-----------------------------------------------------------------------
!     OVERFLOW TEST (ZR/FNU TOO SMALL)
!-----------------------------------------------------------------------
  TEST = D1MACH(1)*1.0D+3
  AC = FNU*TEST
  if (ABS(ZRR) > AC .OR. ABS(ZRI) > AC) go to 15
  ZETA1R = 2.0D0*ABS(LOG(TEST))+FNU
  ZETA1I = 0.0D0
  ZETA2R = FNU
  ZETA2I = 0.0D0
  PHIR = 1.0D0
  PHII = 0.0D0
  return
   15 CONTINUE
  TR = ZRR*RFN
  TI = ZRI*RFN
  SR = CONER + (TR*TR-TI*TI)
  SI = CONEI + (TR*TI+TI*TR)
  call ZSQRT(SR, SI, SRR, SRI)
  STR = CONER + SRR
  STI = CONEI + SRI
  call ZDIV(STR, STI, TR, TI, ZNR, ZNI)
  call ZLOG(ZNR, ZNI, STR, STI, IDUM)
  ZETA1R = FNU*STR
  ZETA1I = FNU*STI
  ZETA2R = FNU*SRR
  ZETA2I = FNU*SRI
  call ZDIV(CONER, CONEI, SRR, SRI, TR, TI)
  SRR = TR*RFN
  SRI = TI*RFN
  call ZSQRT(SRR, SRI, CWRKR(16), CWRKI(16))
  PHIR = CWRKR(16)*CON(IKFLG)
  PHII = CWRKI(16)*CON(IKFLG)
  if (IPMTR /= 0) RETURN
  call ZDIV(CONER, CONEI, SR, SI, T2R, T2I)
  CWRKR(1) = CONER
  CWRKI(1) = CONEI
  CRFNR = CONER
  CRFNI = CONEI
  AC = 1.0D0
  L = 1
  DO 20 K=2,15
    SR = ZEROR
    SI = ZEROI
    DO 10 J=1,K
      L = L + 1
      STR = SR*T2R - SI*T2I + C(L)
      SI = SR*T2I + SI*T2R
      SR = STR
   10   CONTINUE
    STR = CRFNR*SRR - CRFNI*SRI
    CRFNI = CRFNR*SRI + CRFNI*SRR
    CRFNR = STR
    CWRKR(K) = CRFNR*SR - CRFNI*SI
    CWRKI(K) = CRFNR*SI + CRFNI*SR
    AC = AC*RFN
    TEST = ABS(CWRKR(K)) + ABS(CWRKI(K))
    if (AC < TOL .AND. TEST < TOL) go to 30
   20 CONTINUE
  K = 15
   30 CONTINUE
  INIT = K
   40 CONTINUE
  if (IKFLG == 2) go to 60
!-----------------------------------------------------------------------
!     COMPUTE SUM FOR THE I FUNCTION
!-----------------------------------------------------------------------
  SR = ZEROR
  SI = ZEROI
  DO 50 I=1,INIT
    SR = SR + CWRKR(I)
    SI = SI + CWRKI(I)
   50 CONTINUE
  SUMR = SR
  SUMI = SI
  PHIR = CWRKR(16)*CON(1)
  PHII = CWRKI(16)*CON(1)
  return
   60 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE SUM FOR THE K FUNCTION
!-----------------------------------------------------------------------
  SR = ZEROR
  SI = ZEROI
  TR = CONER
  DO 70 I=1,INIT
    SR = SR + TR*CWRKR(I)
    SI = SI + TR*CWRKI(I)
    TR = -TR
   70 CONTINUE
  SUMR = SR
  SUMI = SI
  PHIR = CWRKR(16)*CON(2)
  PHII = CWRKI(16)*CON(2)
  return
end