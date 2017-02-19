subroutine ZACAI (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, TOL, &
     ELIM, ALIM)
!
!! ZACAI is subsidiary to ZAIRY
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CACAI-A, ZACAI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
!
!         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
!                 MP=PI*MR*CMPLX(0.0,1.0)
!
!     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
!     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
!     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
!     RECURRENCE REMOVED. A RECURSIVE call TO ZACON CAN RESULT if ZACON
!     IS CALLED FROM ZAIRY.
!
!***SEE ALSO  ZAIRY
!***ROUTINES CALLED  D1MACH, ZABS, ZASYI, ZBKNU, ZMLRI, ZS1S2, ZSERI
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ZACAI
!     COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY
  DOUBLE PRECISION ALIM, ARG, ASCLE, AZ, CSGNR, CSGNI, CSPNR, &
   CSPNI, C1R, C1I, C2R, C2I, CYR, CYI, DFNU, ELIM, FMR, FNU, PI, &
   RL, SGN, TOL, YY, YR, YI, ZR, ZI, ZNR, ZNI, D1MACH, ZABS
  INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
  DIMENSION YR(N), YI(N), CYR(2), CYI(2)
  EXTERNAL ZABS
  DATA PI / 3.14159265358979324D0 /
!***FIRST EXECUTABLE STATEMENT  ZACAI
  NZ = 0
  ZNR = -ZR
  ZNI = -ZI
  AZ = ZABS(ZR,ZI)
  NN = N
  DFNU = FNU + (N-1)
  if (AZ <= 2.0D0) go to 10
  if (AZ*AZ*0.25D0 > DFNU+1.0D0) go to 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     POWER SERIES FOR THE I FUNCTION
!-----------------------------------------------------------------------
  call ZSERI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL, ELIM, ALIM)
  go to 40
   20 CONTINUE
  if (AZ < RL) go to 30
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
!-----------------------------------------------------------------------
  call ZASYI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, TOL, ELIM, &
   ALIM)
  if (NW < 0) go to 80
  go to 40
   30 CONTINUE
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
!-----------------------------------------------------------------------
  call ZMLRI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL)
  if ( NW < 0) go to 80
   40 CONTINUE
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
!-----------------------------------------------------------------------
  call ZBKNU(ZNR, ZNI, FNU, KODE, 1, CYR, CYI, NW, TOL, ELIM, ALIM)
  if (NW /= 0) go to 80
  FMR = MR
  SGN = -DSIGN(PI,FMR)
  CSGNR = 0.0D0
  CSGNI = SGN
  if (KODE == 1) go to 50
  YY = -ZNI
  CSGNR = -CSGNI*SIN(YY)
  CSGNI = CSGNI*COS(YY)
   50 CONTINUE
!-----------------------------------------------------------------------
!     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
  INU = FNU
  ARG = (FNU-INU)*SGN
  CSPNR = COS(ARG)
  CSPNI = SIN(ARG)
  if (MOD(INU,2) == 0) go to 60
  CSPNR = -CSPNR
  CSPNI = -CSPNI
   60 CONTINUE
  C1R = CYR(1)
  C1I = CYI(1)
  C2R = YR(1)
  C2I = YI(1)
  if (KODE == 1) go to 70
  IUF = 0
  ASCLE = 1.0D+3*D1MACH(1)/TOL
  call ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
  NZ = NZ + NW
   70 CONTINUE
  YR(1) = CSPNR*C1R - CSPNI*C1I + CSGNR*C2R - CSGNI*C2I
  YI(1) = CSPNR*C1I + CSPNI*C1R + CSGNR*C2I + CSGNI*C2R
  return
   80 CONTINUE
  NZ = -1
  if ( NW == (-2)) NZ=-2
  return
end
