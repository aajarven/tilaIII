subroutine CACAI (Z, FNU, KODE, MR, N, Y, NZ, RL, TOL, ELIM, ALIM)
!
!! CACAI is subsidiary to CAIRY.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CACAI-A, ZACAI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
!
!         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
!                 MP=PI*MR*CMPLX(0.0,1.0)
!
!     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
!     HALF Z PLANE FOR USE WITH CAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
!     CACAI IS THE SAME AS CACON WITH THE PARTS FOR LARGER ORDERS AND
!     RECURRENCE REMOVED. A RECURSIVE call TO CACON CAN RESULT if CACON
!     IS CALLED FROM CAIRY.
!
!***SEE ALSO  CAIRY
!***ROUTINES CALLED  CASYI, CBKNU, CMLRI, CS1S2, CSERI, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CACAI
  COMPLEX CSGN, CSPN, C1, C2, Y, Z, ZN, CY
  REAL ALIM, ARG, ASCLE, AZ, CPN, DFNU, ELIM, FMR, FNU, PI, RL, &
   SGN, SPN, TOL, YY, R1MACH
  INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
  DIMENSION Y(N), CY(2)
  DATA PI / 3.14159265358979324E0 /
!***FIRST EXECUTABLE STATEMENT  CACAI
  NZ = 0
  ZN = -Z
  AZ = ABS(Z)
  NN = N
  DFNU = FNU + (N-1)
  if (AZ <= 2.0E0) go to 10
  if (AZ*AZ*0.25E0 > DFNU+1.0E0) go to 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     POWER SERIES FOR THE I FUNCTION
!-----------------------------------------------------------------------
  call CSERI(ZN, FNU, KODE, NN, Y, NW, TOL, ELIM, ALIM)
  go to 40
   20 CONTINUE
  if (AZ < RL) go to 30
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
!-----------------------------------------------------------------------
  call CASYI(ZN, FNU, KODE, NN, Y, NW, RL, TOL, ELIM, ALIM)
  if (NW < 0) go to 70
  go to 40
   30 CONTINUE
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
!-----------------------------------------------------------------------
  call CMLRI(ZN, FNU, KODE, NN, Y, NW, TOL)
  if ( NW < 0) go to 70
   40 CONTINUE
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
!-----------------------------------------------------------------------
  call CBKNU(ZN, FNU, KODE, 1, CY, NW, TOL, ELIM, ALIM)
  if (NW /= 0) go to 70
  FMR = MR
  SGN = -SIGN(PI,FMR)
  CSGN = CMPLX(0.0E0,SGN)
  if (KODE == 1) go to 50
  YY = -AIMAG(ZN)
  CPN = COS(YY)
  SPN = SIN(YY)
  CSGN = CSGN*CMPLX(CPN,SPN)
   50 CONTINUE
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
  C1 = CY(1)
  C2 = Y(1)
  if (KODE == 1) go to 60
  IUF = 0
  ASCLE = 1.0E+3*R1MACH(1)/TOL
  call CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
  NZ = NZ + NW
   60 CONTINUE
  Y(1) = CSPN*C1 + CSGN*C2
  return
   70 CONTINUE
  NZ = -1
  if ( NW == (-2)) NZ=-2
  return
end
