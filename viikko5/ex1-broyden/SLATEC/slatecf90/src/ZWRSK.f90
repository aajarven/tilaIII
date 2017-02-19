subroutine ZWRSK (ZRR, ZRI, FNU, KODE, N, YR, YI, NZ, CWR, CWI, &
     TOL, ELIM, ALIM)
!
!! ZWRSK is subsidiary to ZBESI and ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CWRSK-A, ZWRSK-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY
!     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
!
!***SEE ALSO  ZBESI, ZBESK
!***ROUTINES CALLED  D1MACH, ZABS, ZBKNU, ZRATI
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ZWRSK
!     COMPLEX CINU,CSCL,CT,CW,C1,C2,RCT,ST,Y,ZR
  DOUBLE PRECISION ACT, ACW, ALIM, ASCLE, CINUI, CINUR, CSCLR, CTI, &
   CTR, CWI, CWR, C1I, C1R, C2I, C2R, ELIM, FNU, PTI, PTR, RACT, &
   STI, STR, TOL, YI, YR, ZRI, ZRR, ZABS, D1MACH
  INTEGER I, KODE, N, NW, NZ
  DIMENSION YR(N), YI(N), CWR(2), CWI(2)
  EXTERNAL ZABS
!***FIRST EXECUTABLE STATEMENT  ZWRSK
!-----------------------------------------------------------------------
!     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
!     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
!     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
!-----------------------------------------------------------------------
!
  NZ = 0
  call ZBKNU(ZRR, ZRI, FNU, KODE, 2, CWR, CWI, NW, TOL, ELIM, ALIM)
  if (NW /= 0) go to 50
  call ZRATI(ZRR, ZRI, FNU, N, YR, YI, TOL)
!-----------------------------------------------------------------------
!     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
!     R(FNU+J-1,Z)=Y(J),  J=1,...,N
!-----------------------------------------------------------------------
  CINUR = 1.0D0
  CINUI = 0.0D0
  if (KODE == 1) go to 10
  CINUR = COS(ZRI)
  CINUI = SIN(ZRI)
   10 CONTINUE
!-----------------------------------------------------------------------
!     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
!     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
!     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
!     THE RESULT IS ON SCALE.
!-----------------------------------------------------------------------
  ACW = ZABS(CWR(2),CWI(2))
  ASCLE = 1.0D+3*D1MACH(1)/TOL
  CSCLR = 1.0D0
  if (ACW > ASCLE) go to 20
  CSCLR = 1.0D0/TOL
  go to 30
   20 CONTINUE
  ASCLE = 1.0D0/ASCLE
  if (ACW < ASCLE) go to 30
  CSCLR = TOL
   30 CONTINUE
  C1R = CWR(1)*CSCLR
  C1I = CWI(1)*CSCLR
  C2R = CWR(2)*CSCLR
  C2I = CWI(2)*CSCLR
  STR = YR(1)
  STI = YI(1)
!-----------------------------------------------------------------------
!     CINU=CINU*(CONJG(CT)/ABS(CT))*(1.0D0/ABS(CT) PREVENTS
!     UNDER- OR OVERFLOW PREMATURELY BY SQUARING ABS(CT)
!-----------------------------------------------------------------------
  PTR = STR*C1R - STI*C1I
  PTI = STR*C1I + STI*C1R
  PTR = PTR + C2R
  PTI = PTI + C2I
  CTR = ZRR*PTR - ZRI*PTI
  CTI = ZRR*PTI + ZRI*PTR
  ACT = ZABS(CTR,CTI)
  RACT = 1.0D0/ACT
  CTR = CTR*RACT
  CTI = -CTI*RACT
  PTR = CINUR*RACT
  PTI = CINUI*RACT
  CINUR = PTR*CTR - PTI*CTI
  CINUI = PTR*CTI + PTI*CTR
  YR(1) = CINUR*CSCLR
  YI(1) = CINUI*CSCLR
  if (N == 1) RETURN
  DO 40 I=2,N
    PTR = STR*CINUR - STI*CINUI
    CINUI = STR*CINUI + STI*CINUR
    CINUR = PTR
    STR = YR(I)
    STI = YI(I)
    YR(I) = CINUR*CSCLR
    YI(I) = CINUI*CSCLR
   40 CONTINUE
  return
   50 CONTINUE
  NZ = -1
  if ( NW == (-2)) NZ=-2
  return
end
