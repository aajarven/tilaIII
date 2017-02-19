subroutine CBUNI (Z, FNU, KODE, N, Y, NZ, NUI, NLAST, FNUL, TOL, &
     ELIM, ALIM)
!
!! CBUNI is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CBUNI-A, ZBUNI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE ABS(Z) >
!     FNUL AND FNU+N-1 < FNUL. THE ORDER IS INCREASED FROM
!     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
!     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
!     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  CUNI1, CUNI2, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CBUNI
  COMPLEX CSCL, CSCR, CY, RZ, ST, S1, S2, Y, Z
  REAL ALIM, AX, AY, DFNU, ELIM, FNU, FNUI, FNUL, GNU, TOL, XX, YY, &
   ASCLE, BRY, STR, STI, STM, R1MACH
  INTEGER I, IFLAG, IFORM, K, KODE, N, NL, NLAST, NUI, NW, NZ
  DIMENSION Y(N), CY(2), BRY(3)
!***FIRST EXECUTABLE STATEMENT  CBUNI
  NZ = 0
  XX = REAL(Z)
  YY = AIMAG(Z)
  AX = ABS(XX)*1.7321E0
  AY = ABS(YY)
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
  call CUNI1(Z, GNU, KODE, 2, CY, NW, NLAST, FNUL, TOL, ELIM, ALIM)
  go to 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
!     AND HPI=PI/2
!-----------------------------------------------------------------------
  call CUNI2(Z, GNU, KODE, 2, CY, NW, NLAST, FNUL, TOL, ELIM, ALIM)
   20 CONTINUE
  if (NW < 0) go to 50
  if (NW /= 0) go to 90
  AY = ABS(CY(1))
!----------------------------------------------------------------------
!     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
!----------------------------------------------------------------------
  BRY(1) = 1.0E+3*R1MACH(1)/TOL
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = BRY(2)
  IFLAG = 2
  ASCLE = BRY(2)
  AX = 1.0E0
  CSCL = CMPLX(AX,0.0E0)
  if (AY > BRY(1)) go to 21
  IFLAG = 1
  ASCLE = BRY(1)
  AX = 1.0E0/TOL
  CSCL = CMPLX(AX,0.0E0)
  go to 25
   21 CONTINUE
  if (AY < BRY(2)) go to 25
  IFLAG = 3
  ASCLE = BRY(3)
  AX = TOL
  CSCL = CMPLX(AX,0.0E0)
   25 CONTINUE
  AY = 1.0E0/AX
  CSCR = CMPLX(AY,0.0E0)
  S1 = CY(2)*CSCL
  S2 = CY(1)*CSCL
  RZ = CMPLX(2.0E0,0.0E0)/Z
  DO 30 I=1,NUI
    ST = S2
    S2 = CMPLX(DFNU+FNUI,0.0E0)*RZ*S2 + S1
    S1 = ST
    FNUI = FNUI - 1.0E0
    if (IFLAG >= 3) go to 30
    ST = S2*CSCR
    STR = REAL(ST)
    STI = AIMAG(ST)
    STR = ABS(STR)
    STI = ABS(STI)
    STM = MAX(STR,STI)
    if (STM <= ASCLE) go to 30
    IFLAG = IFLAG+1
    ASCLE = BRY(IFLAG)
    S1 = S1*CSCR
    S2 = ST
    AX = AX*TOL
    AY = 1.0E0/AX
    CSCL = CMPLX(AX,0.0E0)
    CSCR = CMPLX(AY,0.0E0)
    S1 = S1*CSCL
    S2 = S2*CSCL
   30 CONTINUE
  Y(N) = S2*CSCR
  if (N == 1) RETURN
  NL = N - 1
  FNUI = NL
  K = NL
  DO 40 I=1,NL
    ST = S2
    S2 = CMPLX(FNU+FNUI,0.0E0)*RZ*S2 + S1
    S1 = ST
    ST = S2*CSCR
    Y(K) = ST
    FNUI = FNUI - 1.0E0
    K = K - 1
    if (IFLAG >= 3) go to 40
    STR = REAL(ST)
    STI = AIMAG(ST)
    STR = ABS(STR)
    STI = ABS(STI)
    STM = MAX(STR,STI)
    if (STM <= ASCLE) go to 40
    IFLAG = IFLAG+1
    ASCLE = BRY(IFLAG)
    S1 = S1*CSCR
    S2 = ST
    AX = AX*TOL
    AY = 1.0E0/AX
    CSCL = CMPLX(AX,0.0E0)
    CSCR = CMPLX(AY,0.0E0)
    S1 = S1*CSCL
    S2 = S2*CSCL
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
  call CUNI1(Z, FNU, KODE, N, Y, NW, NLAST, FNUL, TOL, ELIM, ALIM)
  go to 80
   70 CONTINUE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
!     AND HPI=PI/2
!-----------------------------------------------------------------------
  call CUNI2(Z, FNU, KODE, N, Y, NW, NLAST, FNUL, TOL, ELIM, ALIM)
   80 CONTINUE
  if (NW < 0) go to 50
  NZ = NW
  return
   90 CONTINUE
  NLAST = N
  return
end
