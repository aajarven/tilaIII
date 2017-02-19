subroutine CACON (Z, FNU, KODE, MR, N, Y, NZ, RL, FNUL, TOL, ELIM, &
     ALIM)
!
!! CACON is subsidiary to CBESH and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CACON-A, ZACON-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CACON APPLIES THE ANALYTIC CONTINUATION FORMULA
!
!         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
!                 MP=PI*MR*CMPLX(0.0,1.0)
!
!     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
!     HALF Z PLANE
!
!***SEE ALSO  CBESH, CBESK
!***ROUTINES CALLED  CBINU, CBKNU, CS1S2, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CACON
  COMPLEX CK, CONE, CS, CSCL, CSCR, CSGN, CSPN, CSS, CSR, C1, C2, &
   RZ, SC1, SC2, ST, S1, S2, Y, Z, ZN, CY
  REAL ALIM, ARG, ASCLE, AS2, BSCLE, BRY, CPN, C1I, C1M, C1R, ELIM, &
   FMR, FNU, FNUL, PI, RL, SGN, SPN, TOL, YY, R1MACH
  INTEGER I, INU, IUF, KFLAG, KODE, MR, N, NN, NW, NZ
  DIMENSION Y(N), CY(2), CSS(3), CSR(3), BRY(3)
  DATA PI / 3.14159265358979324E0 /
  DATA CONE / (1.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CACON
  NZ = 0
  ZN = -Z
  NN = N
  call CBINU(ZN, FNU, KODE, NN, Y, NW, RL, FNUL, TOL, ELIM, ALIM)
  if (NW < 0) go to 80
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
!-----------------------------------------------------------------------
  NN = MIN(2,N)
  call CBKNU(ZN, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
  if (NW /= 0) go to 80
  S1 = CY(1)
  FMR = MR
  SGN = -SIGN(PI,FMR)
  CSGN = CMPLX(0.0E0,SGN)
  if (KODE == 1) go to 10
  YY = -AIMAG(ZN)
  CPN = COS(YY)
  SPN = SIN(YY)
  CSGN = CSGN*CMPLX(CPN,SPN)
   10 CONTINUE
!-----------------------------------------------------------------------
!     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
  INU = FNU
  ARG = (FNU-INU)*SGN
  CPN = COS(ARG)
  SPN = SIN(ARG)
  CSPN = CMPLX(CPN,SPN)
  if (MOD(INU,2) == 1) CSPN = -CSPN
  IUF = 0
  C1 = S1
  C2 = Y(1)
  ASCLE = 1.0E+3*R1MACH(1)/TOL
  if (KODE == 1) go to 20
  call CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
  NZ = NZ + NW
  SC1 = C1
   20 CONTINUE
  Y(1) = CSPN*C1 + CSGN*C2
  if (N == 1) RETURN
  CSPN = -CSPN
  S2 = CY(2)
  C1 = S2
  C2 = Y(2)
  if (KODE == 1) go to 30
  call CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
  NZ = NZ + NW
  SC2 = C1
   30 CONTINUE
  Y(2) = CSPN*C1 + CSGN*C2
  if (N == 2) RETURN
  CSPN = -CSPN
  RZ = CMPLX(2.0E0,0.0E0)/ZN
  CK = CMPLX(FNU+1.0E0,0.0E0)*RZ
!-----------------------------------------------------------------------
!     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
!-----------------------------------------------------------------------
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CSCR = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CSCR
  CSR(1) = CSCR
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = ASCLE
  BRY(2) = 1.0E0/ASCLE
  BRY(3) = R1MACH(2)
  AS2 = ABS(S2)
  KFLAG = 2
  if (AS2 > BRY(1)) go to 40
  KFLAG = 1
  go to 50
   40 CONTINUE
  if (AS2 < BRY(2)) go to 50
  KFLAG = 3
   50 CONTINUE
  BSCLE = BRY(KFLAG)
  S1 = S1*CSS(KFLAG)
  S2 = S2*CSS(KFLAG)
  CS = CSR(KFLAG)
  DO 70 I=3,N
    ST = S2
    S2 = CK*S2 + S1
    S1 = ST
    C1 = S2*CS
    ST = C1
    C2 = Y(I)
    if (KODE == 1) go to 60
    if (IUF < 0) go to 60
    call CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
    NZ = NZ + NW
    SC1 = SC2
    SC2 = C1
    if (IUF /= 3) go to 60
    IUF = -4
    S1 = SC1*CSS(KFLAG)
    S2 = SC2*CSS(KFLAG)
    ST = SC2
   60   CONTINUE
    Y(I) = CSPN*C1 + CSGN*C2
    CK = CK + RZ
    CSPN = -CSPN
    if (KFLAG >= 3) go to 70
    C1R = REAL(C1)
    C1I = AIMAG(C1)
    C1R = ABS(C1R)
    C1I = ABS(C1I)
    C1M = MAX(C1R,C1I)
    if (C1M <= BSCLE) go to 70
    KFLAG = KFLAG + 1
    BSCLE = BRY(KFLAG)
    S1 = S1*CS
    S2 = ST
    S1 = S1*CSS(KFLAG)
    S2 = S2*CSS(KFLAG)
    CS = CSR(KFLAG)
   70 CONTINUE
  return
   80 CONTINUE
  NZ = -1
  if ( NW == (-2)) NZ=-2
  return
end
