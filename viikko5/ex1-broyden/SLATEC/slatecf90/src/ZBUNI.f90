subroutine ZBUNI (ZR, ZI, FNU, KODE, N, YR, YI, NZ, NUI, NLAST, &
     FNUL, TOL, ELIM, ALIM)
!
!! ZBUNI is subsidiary to ZBESI and ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CBUNI-A, ZBUNI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE ABS(Z) >
!     FNUL AND FNU+N-1 < FNUL. THE ORDER IS INCREASED FROM
!     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
!     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
!     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
!
!***SEE ALSO  ZBESI, ZBESK
!***ROUTINES CALLED  D1MACH, ZABS, ZUNI1, ZUNI2
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ZBUNI
!     COMPLEX CSCL,CSCR,CY,RZ,ST,S1,S2,Y,Z
  DOUBLE PRECISION ALIM, AX, AY, CSCLR, CSCRR, CYI, CYR, DFNU, &
   ELIM, FNU, FNUI, FNUL, GNU, RAZ, RZI, RZR, STI, STR, S1I, S1R, &
   S2I, S2R, TOL, YI, YR, ZI, ZR, ZABS, ASCLE, BRY, C1R, C1I, C1M, &
   D1MACH
  INTEGER I, IFLAG, IFORM, K, KODE, N, NL, NLAST, NUI, NW, NZ
  DIMENSION YR(N), YI(N), CYR(2), CYI(2), BRY(3)
  EXTERNAL ZABS
!***FIRST EXECUTABLE STATEMENT  ZBUNI
  NZ = 0
  AX = ABS(ZR)*1.7321D0
  AY = ABS(ZI)
  IFORM = 1
  if (AY > AX) IFORM = 2
  if (NUI == 0) go to 60
  FNUI = NUI
  DFNU = FNU + (N-1)
  GNU = DFNU + FNUI
  if (IFORM == 2) go to 10
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3 <= ARG(Z) <= PI/3
!-----------------------------------------------------------------------
  call ZUNI1(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL, &
   ELIM, ALIM)
  go to 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
!     AND HPI=PI/2
!-----------------------------------------------------------------------
  call ZUNI2(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL, &
   ELIM, ALIM)
   20 CONTINUE
  if (NW < 0) go to 50
  if (NW /= 0) go to 90
  STR = ZABS(CYR(1),CYI(1))
!----------------------------------------------------------------------
!     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
!----------------------------------------------------------------------
  BRY(1)=1.0D+3*D1MACH(1)/TOL
  BRY(2) = 1.0D0/BRY(1)
  BRY(3) = BRY(2)
  IFLAG = 2
  ASCLE = BRY(2)
  CSCLR = 1.0D0
  if (STR > BRY(1)) go to 21
  IFLAG = 1
  ASCLE = BRY(1)
  CSCLR = 1.0D0/TOL
  go to 25
   21 CONTINUE
  if (STR < BRY(2)) go to 25
  IFLAG = 3
  ASCLE=BRY(3)
  CSCLR = TOL
   25 CONTINUE
  CSCRR = 1.0D0/CSCLR
  S1R = CYR(2)*CSCLR
  S1I = CYI(2)*CSCLR
  S2R = CYR(1)*CSCLR
  S2I = CYI(1)*CSCLR
  RAZ = 1.0D0/ZABS(ZR,ZI)
  STR = ZR*RAZ
  STI = -ZI*RAZ
  RZR = (STR+STR)*RAZ
  RZI = (STI+STI)*RAZ
  DO 30 I=1,NUI
    STR = S2R
    STI = S2I
    S2R = (DFNU+FNUI)*(RZR*STR-RZI*STI) + S1R
    S2I = (DFNU+FNUI)*(RZR*STI+RZI*STR) + S1I
    S1R = STR
    S1I = STI
    FNUI = FNUI - 1.0D0
    if (IFLAG >= 3) go to 30
    STR = S2R*CSCRR
    STI = S2I*CSCRR
    C1R = ABS(STR)
    C1I = ABS(STI)
    C1M = MAX(C1R,C1I)
    if (C1M <= ASCLE) go to 30
    IFLAG = IFLAG+1
    ASCLE = BRY(IFLAG)
    S1R = S1R*CSCRR
    S1I = S1I*CSCRR
    S2R = STR
    S2I = STI
    CSCLR = CSCLR*TOL
    CSCRR = 1.0D0/CSCLR
    S1R = S1R*CSCLR
    S1I = S1I*CSCLR
    S2R = S2R*CSCLR
    S2I = S2I*CSCLR
   30 CONTINUE
  YR(N) = S2R*CSCRR
  YI(N) = S2I*CSCRR
  if (N == 1) RETURN
  NL = N - 1
  FNUI = NL
  K = NL
  DO 40 I=1,NL
    STR = S2R
    STI = S2I
    S2R = (FNU+FNUI)*(RZR*STR-RZI*STI) + S1R
    S2I = (FNU+FNUI)*(RZR*STI+RZI*STR) + S1I
    S1R = STR
    S1I = STI
    STR = S2R*CSCRR
    STI = S2I*CSCRR
    YR(K) = STR
    YI(K) = STI
    FNUI = FNUI - 1.0D0
    K = K - 1
    if (IFLAG >= 3) go to 40
    C1R = ABS(STR)
    C1I = ABS(STI)
    C1M = MAX(C1R,C1I)
    if (C1M <= ASCLE) go to 40
    IFLAG = IFLAG+1
    ASCLE = BRY(IFLAG)
    S1R = S1R*CSCRR
    S1I = S1I*CSCRR
    S2R = STR
    S2I = STI
    CSCLR = CSCLR*TOL
    CSCRR = 1.0D0/CSCLR
    S1R = S1R*CSCLR
    S1I = S1I*CSCLR
    S2R = S2R*CSCLR
    S2I = S2I*CSCLR
   40 CONTINUE
  return
   50 CONTINUE
  NZ = -1
  if ( NW == (-2)) NZ=-2
  return
   60 CONTINUE
  if (IFORM == 2) go to 70
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3 <= ARG(Z) <= PI/3
!-----------------------------------------------------------------------
  call ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL, &
   ELIM, ALIM)
  go to 80
   70 CONTINUE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
!     AND HPI=PI/2
!-----------------------------------------------------------------------
  call ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL, &
   ELIM, ALIM)
   80 CONTINUE
  if (NW < 0) go to 50
  NZ = NW
  return
   90 CONTINUE
  NLAST = N
  return
end
