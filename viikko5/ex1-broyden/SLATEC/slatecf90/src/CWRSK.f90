subroutine CWRSK (ZR, FNU, KODE, N, Y, NZ, CW, TOL, ELIM, ALIM)
!
!! CWRSK is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CWRSK-A, ZWRSK-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY
!     NORMALIZING THE I FUNCTION RATIOS FROM CRATI BY THE WRONSKIAN
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  CBKNU, CRATI, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CWRSK
  COMPLEX CINU, CSCL, CT, CW, C1, C2, RCT, ST, Y, ZR
  REAL ACT, ACW, ALIM, ASCLE, ELIM, FNU, S1, S2, TOL, YY, R1MACH
  INTEGER I, KODE, N, NW, NZ
  DIMENSION Y(N), CW(2)
!***FIRST EXECUTABLE STATEMENT  CWRSK
!
!     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
!     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
!     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
!
  NZ = 0
  call CBKNU(ZR, FNU, KODE, 2, CW, NW, TOL, ELIM, ALIM)
  if (NW /= 0) go to 50
  call CRATI(ZR, FNU, N, Y, TOL)
!-----------------------------------------------------------------------
!     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
!     R(FNU+J-1,Z)=Y(J),  J=1,...,N
!-----------------------------------------------------------------------
  CINU = CMPLX(1.0E0,0.0E0)
  if (KODE == 1) go to 10
  YY = AIMAG(ZR)
  S1 = COS(YY)
  S2 = SIN(YY)
  CINU = CMPLX(S1,S2)
   10 CONTINUE
!-----------------------------------------------------------------------
!     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
!     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
!     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
!     THE RESULT IS ON SCALE.
!-----------------------------------------------------------------------
  ACW = ABS(CW(2))
  ASCLE = 1.0E+3*R1MACH(1)/TOL
  CSCL = CMPLX(1.0E0,0.0E0)
  if (ACW > ASCLE) go to 20
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  go to 30
   20 CONTINUE
  ASCLE = 1.0E0/ASCLE
  if (ACW < ASCLE) go to 30
  CSCL = CMPLX(TOL,0.0E0)
   30 CONTINUE
  C1 = CW(1)*CSCL
  C2 = CW(2)*CSCL
  ST = Y(1)
!-----------------------------------------------------------------------
!     CINU=CINU*(CONJG(CT)/ABS(CT))*(1.0E0/ABS(CT) PREVENTS
!     UNDER- OR OVERFLOW PREMATURELY BY SQUARING ABS(CT)
!-----------------------------------------------------------------------
  CT = ZR*(C2+ST*C1)
  ACT = ABS(CT)
  RCT = CMPLX(1.0E0/ACT,0.0E0)
  CT = CONJG(CT)*RCT
  CINU = CINU*RCT*CT
  Y(1) = CINU*CSCL
  if (N == 1) RETURN
  DO 40 I=2,N
    CINU = ST*CINU
    ST = Y(I)
    Y(I) = CINU*CSCL
   40 CONTINUE
  return
   50 CONTINUE
  NZ = -1
  if ( NW == (-2)) NZ=-2
  return
end
