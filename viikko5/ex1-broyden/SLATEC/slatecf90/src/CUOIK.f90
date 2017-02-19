subroutine CUOIK (Z, FNU, KODE, IKFLG, N, Y, NUF, TOL, ELIM, ALIM)
!
!! CUOIK is subsidiary to CBESH, CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUOIK-A, ZUOIK-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
!     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
!     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
!     WHERE ALIM < ELIM. if THE MAGNITUDE, BASED ON THE LEADING
!     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
!     THE RESULT IS ON SCALE. if NOT, THEN A REFINED TEST USING OTHER
!     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
!     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
!     EXP(-ELIM)/TOL
!
!     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
!          =2 MEANS THE K SEQUENCE IS TESTED
!     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
!         =-1 MEANS AN OVERFLOW WOULD OCCUR
!     IKFLG=1 AND NUF > 0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
!             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
!     IKFLG=2 AND NUF == N MEANS ALL Y VALUES WERE SET TO ZERO
!     IKFLG=2 AND 0 < NUF < N NOT CONSIDERED. Y MUST BE SET BY
!             ANOTHER ROUTINE
!
!***SEE ALSO  CBESH, CBESI, CBESK
!***ROUTINES CALLED  CUCHK, CUNHJ, CUNIK, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUOIK
  COMPLEX ARG, ASUM, BSUM, CWRK, CZ, CZERO, PHI, SUM, Y, Z, ZB, &
   ZETA1, ZETA2, ZN, ZR
  REAL AARG, AIC, ALIM, APHI, ASCLE, AX, AY, ELIM, FNN, FNU, GNN, &
   GNU, RCZ, TOL, X, YY, R1MACH
  INTEGER I, IFORM, IKFLG, INIT, KODE, N, NN, NUF, NW
  DIMENSION Y(N), CWRK(16)
  DATA CZERO / (0.0E0,0.0E0) /
  DATA AIC / 1.265512123484645396E+00 /
!***FIRST EXECUTABLE STATEMENT  CUOIK
  NUF = 0
  NN = N
  X = REAL(Z)
  ZR = Z
  if (X < 0.0E0) ZR = -Z
  ZB = ZR
  YY = AIMAG(ZR)
  AX = ABS(X)*1.7321E0
  AY = ABS(YY)
  IFORM = 1
  if (AY > AX) IFORM = 2
  GNU = MAX(FNU,1.0E0)
  if (IKFLG == 1) go to 10
  FNN = NN
  GNN = FNU + FNN - 1.0E0
  GNU = MAX(GNN,FNN)
   10 CONTINUE
!-----------------------------------------------------------------------
!     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
!     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
!     THE SIGN OF THE IMAGINARY PART CORRECT.
!-----------------------------------------------------------------------
  if (IFORM == 2) go to 20
  INIT = 0
  call CUNIK(ZR, GNU, IKFLG, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM, &
   CWRK)
  CZ = -ZETA1 + ZETA2
  go to 40
   20 CONTINUE
  ZN = -ZR*CMPLX(0.0E0,1.0E0)
  if (YY > 0.0E0) go to 30
  ZN = CONJG(-ZN)
   30 CONTINUE
  call CUNHJ(ZN, GNU, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
  CZ = -ZETA1 + ZETA2
  AARG = ABS(ARG)
   40 CONTINUE
  if (KODE == 2) CZ = CZ - ZB
  if (IKFLG == 2) CZ = -CZ
  APHI = ABS(PHI)
  RCZ = REAL(CZ)
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
  if (RCZ > ELIM) go to 170
  if (RCZ < ALIM) go to 50
  RCZ = RCZ + ALOG(APHI)
  if (IFORM == 2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
  if (RCZ > ELIM) go to 170
  go to 100
   50 CONTINUE
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
  if (RCZ < (-ELIM)) go to 60
  if (RCZ > (-ALIM)) go to 100
  RCZ = RCZ + ALOG(APHI)
  if (IFORM == 2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
  if (RCZ > (-ELIM)) go to 80
   60 CONTINUE
  DO 70 I=1,NN
    Y(I) = CZERO
   70 CONTINUE
  NUF = NN
  return
   80 CONTINUE
  ASCLE = 1.0E+3*R1MACH(1)/TOL
  CZ = CZ + CLOG(PHI)
  if (IFORM == 1) go to 90
  CZ = CZ - CMPLX(0.25E0,0.0E0)*CLOG(ARG) - CMPLX(AIC,0.0E0)
   90 CONTINUE
  AX = EXP(RCZ)/TOL
  AY = AIMAG(CZ)
  CZ = CMPLX(AX,0.0E0)*CMPLX(COS(AY),SIN(AY))
  call CUCHK(CZ, NW, ASCLE, TOL)
  if (NW == 1) go to 60
  100 CONTINUE
  if (IKFLG == 2) RETURN
  if (N == 1) RETURN
!-----------------------------------------------------------------------
!     SET UNDERFLOWS ON I SEQUENCE
!-----------------------------------------------------------------------
  110 CONTINUE
  GNU = FNU + (NN-1)
  if (IFORM == 2) go to 120
  INIT = 0
  call CUNIK(ZR, GNU, IKFLG, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM, &
   CWRK)
  CZ = -ZETA1 + ZETA2
  go to 130
  120 CONTINUE
  call CUNHJ(ZN, GNU, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
  CZ = -ZETA1 + ZETA2
  AARG = ABS(ARG)
  130 CONTINUE
  if (KODE == 2) CZ = CZ - ZB
  APHI = ABS(PHI)
  RCZ = REAL(CZ)
  if (RCZ < (-ELIM)) go to 140
  if (RCZ > (-ALIM)) RETURN
  RCZ = RCZ + ALOG(APHI)
  if (IFORM == 2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
  if (RCZ > (-ELIM)) go to 150
  140 CONTINUE
  Y(NN) = CZERO
  NN = NN - 1
  NUF = NUF + 1
  if (NN == 0) RETURN
  go to 110
  150 CONTINUE
  ASCLE = 1.0E+3*R1MACH(1)/TOL
  CZ = CZ + CLOG(PHI)
  if (IFORM == 1) go to 160
  CZ = CZ - CMPLX(0.25E0,0.0E0)*CLOG(ARG) - CMPLX(AIC,0.0E0)
  160 CONTINUE
  AX = EXP(RCZ)/TOL
  AY = AIMAG(CZ)
  CZ = CMPLX(AX,0.0E0)*CMPLX(COS(AY),SIN(AY))
  call CUCHK(CZ, NW, ASCLE, TOL)
  if (NW == 1) go to 140
  return
  170 CONTINUE
  NUF = -1
  return
end
