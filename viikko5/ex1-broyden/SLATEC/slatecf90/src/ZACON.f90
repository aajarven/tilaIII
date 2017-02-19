subroutine ZACON (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, FNUL, &
     TOL, ELIM, ALIM)
!
!! ZACON is subsidiary to ZBESH and ZBESK
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CACON-A, ZACON-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZACON APPLIES THE ANALYTIC CONTINUATION FORMULA
!
!         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
!                 MP=PI*MR*CMPLX(0.0,1.0)
!
!     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
!     HALF Z PLANE
!
!***SEE ALSO  ZBESH, ZBESK
!***ROUTINES CALLED  D1MACH, ZABS, ZBINU, ZBKNU, ZMLT, ZS1S2
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ZACON
!     COMPLEX CK,CONE,CSCL,CSCR,CSGN,CSPN,CY,CZERO,C1,C2,RZ,SC1,SC2,ST,
!    *S1,S2,Y,Z,ZN
  DOUBLE PRECISION ALIM, ARG, ASCLE, AS2, AZN, BRY, BSCLE, CKI, &
   CKR, CONER, CPN, CSCL, CSCR, CSGNI, CSGNR, CSPNI, CSPNR, &
   CSR, CSRR, CSSR, CYI, CYR, C1I, C1M, C1R, C2I, C2R, ELIM, FMR, &
   FN, FNU, FNUL, PI, PTI, PTR, RAZN, RL, RZI, RZR, SC1I, SC1R, &
   SC2I, SC2R, SGN, SPN, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, &
   YY, ZEROR, ZI, ZNI, ZNR, ZR, D1MACH, ZABS
  INTEGER I, INU, IUF, KFLAG, KODE, MR, N, NN, NW, NZ
  DIMENSION YR(N), YI(N), CYR(2), CYI(2), CSSR(3), CSRR(3), BRY(3)
  EXTERNAL ZABS
  DATA PI / 3.14159265358979324D0 /
  DATA ZEROR,CONER / 0.0D0,1.0D0 /
!***FIRST EXECUTABLE STATEMENT  ZACON
  NZ = 0
  ZNR = -ZR
  ZNI = -ZI
  NN = N
  call ZBINU(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, FNUL, TOL, &
   ELIM, ALIM)
  if (NW < 0) go to 90
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
!-----------------------------------------------------------------------
  NN = MIN(2,N)
  call ZBKNU(ZNR, ZNI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
  if (NW /= 0) go to 90
  S1R = CYR(1)
  S1I = CYI(1)
  FMR = MR
  SGN = -DSIGN(PI,FMR)
  CSGNR = ZEROR
  CSGNI = SGN
  if (KODE == 1) go to 10
  YY = -ZNI
  CPN = COS(YY)
  SPN = SIN(YY)
  call ZMLT(CSGNR, CSGNI, CPN, SPN, CSGNR, CSGNI)
   10 CONTINUE
!-----------------------------------------------------------------------
!     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
  INU = FNU
  ARG = (FNU-INU)*SGN
  CPN = COS(ARG)
  SPN = SIN(ARG)
  CSPNR = CPN
  CSPNI = SPN
  if (MOD(INU,2) == 0) go to 20
  CSPNR = -CSPNR
  CSPNI = -CSPNI
   20 CONTINUE
  IUF = 0
  C1R = S1R
  C1I = S1I
  C2R = YR(1)
  C2I = YI(1)
  ASCLE = 1.0D+3*D1MACH(1)/TOL
  if (KODE == 1) go to 30
  call ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
  NZ = NZ + NW
  SC1R = C1R
  SC1I = C1I
   30 CONTINUE
  call ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI)
  call ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI)
  YR(1) = STR + PTR
  YI(1) = STI + PTI
  if (N == 1) RETURN
  CSPNR = -CSPNR
  CSPNI = -CSPNI
  S2R = CYR(2)
  S2I = CYI(2)
  C1R = S2R
  C1I = S2I
  C2R = YR(2)
  C2I = YI(2)
  if (KODE == 1) go to 40
  call ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
  NZ = NZ + NW
  SC2R = C1R
  SC2I = C1I
   40 CONTINUE
  call ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI)
  call ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI)
  YR(2) = STR + PTR
  YI(2) = STI + PTI
  if (N == 2) RETURN
  CSPNR = -CSPNR
  CSPNI = -CSPNI
  AZN = ZABS(ZNR,ZNI)
  RAZN = 1.0D0/AZN
  STR = ZNR*RAZN
  STI = -ZNI*RAZN
  RZR = (STR+STR)*RAZN
  RZI = (STI+STI)*RAZN
  FN = FNU + 1.0D0
  CKR = FN*RZR
  CKI = FN*RZI
!-----------------------------------------------------------------------
!     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
!-----------------------------------------------------------------------
  CSCL = 1.0D0/TOL
  CSCR = TOL
  CSSR(1) = CSCL
  CSSR(2) = CONER
  CSSR(3) = CSCR
  CSRR(1) = CSCR
  CSRR(2) = CONER
  CSRR(3) = CSCL
  BRY(1) = ASCLE
  BRY(2) = 1.0D0/ASCLE
  BRY(3) = D1MACH(2)
  AS2 = ZABS(S2R,S2I)
  KFLAG = 2
  if (AS2 > BRY(1)) go to 50
  KFLAG = 1
  go to 60
   50 CONTINUE
  if (AS2 < BRY(2)) go to 60
  KFLAG = 3
   60 CONTINUE
  BSCLE = BRY(KFLAG)
  S1R = S1R*CSSR(KFLAG)
  S1I = S1I*CSSR(KFLAG)
  S2R = S2R*CSSR(KFLAG)
  S2I = S2I*CSSR(KFLAG)
  CSR = CSRR(KFLAG)
  DO 80 I=3,N
    STR = S2R
    STI = S2I
    S2R = CKR*STR - CKI*STI + S1R
    S2I = CKR*STI + CKI*STR + S1I
    S1R = STR
    S1I = STI
    C1R = S2R*CSR
    C1I = S2I*CSR
    STR = C1R
    STI = C1I
    C2R = YR(I)
    C2I = YI(I)
    if (KODE == 1) go to 70
    if (IUF < 0) go to 70
    call ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
    NZ = NZ + NW
    SC1R = SC2R
    SC1I = SC2I
    SC2R = C1R
    SC2I = C1I
    if (IUF /= 3) go to 70
    IUF = -4
    S1R = SC1R*CSSR(KFLAG)
    S1I = SC1I*CSSR(KFLAG)
    S2R = SC2R*CSSR(KFLAG)
    S2I = SC2I*CSSR(KFLAG)
    STR = SC2R
    STI = SC2I
   70   CONTINUE
    PTR = CSPNR*C1R - CSPNI*C1I
    PTI = CSPNR*C1I + CSPNI*C1R
    YR(I) = PTR + CSGNR*C2R - CSGNI*C2I
    YI(I) = PTI + CSGNR*C2I + CSGNI*C2R
    CKR = CKR + RZR
    CKI = CKI + RZI
    CSPNR = -CSPNR
    CSPNI = -CSPNI
    if (KFLAG >= 3) go to 80
    PTR = ABS(C1R)
    PTI = ABS(C1I)
    C1M = MAX(PTR,PTI)
    if (C1M <= BSCLE) go to 80
    KFLAG = KFLAG + 1
    BSCLE = BRY(KFLAG)
    S1R = S1R*CSR
    S1I = S1I*CSR
    S2R = STR
    S2I = STI
    S1R = S1R*CSSR(KFLAG)
    S1I = S1I*CSSR(KFLAG)
    S2R = S2R*CSSR(KFLAG)
    S2I = S2I*CSSR(KFLAG)
    CSR = CSRR(KFLAG)
   80 CONTINUE
  return
   90 CONTINUE
  NZ = -1
  if ( NW == (-2)) NZ=-2
  return
end
