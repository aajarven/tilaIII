subroutine CAIRY (Z, ID, KODE, AI, NZ, IERR)
!
!! CAIRY computes the Airy function Ai(z) or its derivative dAi/dz ...
!            for complex argument z.  A scaling option is available ...
!            to help avoid underflow and overflow.
!
!***LIBRARY   SLATEC
!***CATEGORY  C10D
!***TYPE      COMPLEX (CAIRY-C, ZAIRY-C)
!***KEYWORDS  AIRY FUNCTION, BESSEL FUNCTION OF ORDER ONE THIRD,
!             BESSEL FUNCTION OF ORDER TWO THIRDS
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!         On KODE=1, CAIRY computes the complex Airy function Ai(z)
!         or its derivative dAi/dz on ID=0 or ID=1 respectively. On
!         KODE=2, a scaling option exp(zeta)*Ai(z) or exp(zeta)*dAi/dz
!         is provided to remove the exponential decay in -pi/3<arg(z)
!         <pi/3 and the exponential growth in pi/3<abs(arg(z))<pi where
!         zeta=(2/3)*z**(3/2).
!
!         While the Airy functions Ai(z) and dAi/dz are analytic in
!         the whole z-plane, the corresponding scaled functions defined
!         for KODE=2 have a cut along the negative real axis.
!
!         Input
!           Z      - Argument of type COMPLEX
!           ID     - Order of derivative, ID=0 or ID=1
!           KODE   - A parameter to indicate the scaling option
!                    KODE=1  returns
!                            AI=Ai(z)  on ID=0
!                            AI=dAi/dz on ID=1
!                            at z=Z
!                        =2  returns
!                            AI=exp(zeta)*Ai(z)  on ID=0
!                            AI=exp(zeta)*dAi/dz on ID=1
!                            at z=Z where zeta=(2/3)*z**(3/2)
!
!         Output
!           AI     - Result of type COMPLEX
!           NZ     - Underflow indicator
!                    NZ=0    Normal return
!                    NZ=1    AI=0 due to underflow in
!                            -pi/3<arg(Z)<pi/3 on KODE=1
!           IERR   - Error flag
!                    IERR=0  Normal return     - COMPUTATION COMPLETED
!                    IERR=1  Input error       - NO COMPUTATION
!                    IERR=2  Overflow          - NO COMPUTATION
!                            (Re(Z) too large with KODE=1)
!                    IERR=3  Precision warning - COMPUTATION COMPLETED
!                            (Result has less than half precision)
!                    IERR=4  Precision error   - NO COMPUTATION
!                            (Result has no precision)
!                    IERR=5  Algorithmic error - NO COMPUTATION
!                            (Termination condition not met)
!
! *Long Description:
!
!         Ai(z) and dAi/dz are computed from K Bessel functions by
!
!                Ai(z) =  c*sqrt(z)*K(1/3,zeta)
!               dAi/dz = -c*   z   *K(2/3,zeta)
!                    c =  1/(pi*sqrt(3))
!                 zeta =  (2/3)*z**(3/2)
!
!         when abs(z)>1 and from power series when abs(z)<=1.
!
!         In most complex variable computation, one must evaluate ele-
!         mentary functions.  When the magnitude of Z is large, losses
!         of significance by argument reduction occur.  Consequently, if
!         the magnitude of ZETA=(2/3)*Z**(3/2) exceeds U1=SQRT(0.5/UR),
!         then losses exceeding half precision are likely and an error
!         flag IERR=3 is triggered where UR=R1MACH(4)=UNIT ROUNDOFF.
!         Also, if the magnitude of ZETA is larger than U2=0.5/UR, then
!         all significance is lost and IERR=4.  In order to use the INT
!         function, ZETA must be further restricted not to exceed
!         U3=I1MACH(9)=LARGEST INTEGER.  Thus, the magnitude of ZETA
!         must be restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2,
!         and U3 are approximately 2.0E+3, 4.2E+6, 2.1E+9 in single
!         precision and 4.7E+7, 2.3E+15, 2.1E+9 in double precision.
!         This makes U2 limiting is single precision and U3 limiting
!         in double precision.  This means that the magnitude of Z
!         cannot exceed approximately 3.4E+4 in single precision and
!         2.1E+6 in double precision.  This also means that one can
!         expect to retain, in the worst cases on 32-bit machines,
!         no digits in single precision and only 6 digits in double
!         precision.
!
!         The approximate relative error in the magnitude of a complex
!         Bessel function can be expressed as P*10**S where P=MAX(UNIT
!         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre-
!         sents the increase in error due to argument reduction in the
!         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))),
!         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF
!         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may
!         have only absolute accuracy.  This is most likely to occur
!         when one component (in magnitude) is larger than the other by
!         several orders of magnitude.  If one component is 10**K larger
!         than the other, then one can expect only MAX(ABS(LOG10(P))-K,
!         0) significant digits; or, stated another way, when K exceeds
!         the exponent of P, no significant digits remain in the smaller
!         component.  However, the phase angle retains absolute accuracy
!         because, in complex arithmetic with precision P, the smaller
!         component will not (as a rule) decrease below P times the
!         magnitude of the larger component. In these extreme cases,
!         the principal phase angle is on the order of +P, -P, PI/2-P,
!         or -PI/2+P.
!
!***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
!                 matical Functions, National Bureau of Standards
!                 Applied Mathematics Series 55, U. S. Department
!                 of Commerce, Tenth Printing (1972) or later.
!               2. D. E. Amos, Computation of Bessel Functions of
!                 Complex Argument and Large Order, Report SAND83-0643,
!                 Sandia National Laboratories, Albuquerque, NM, May
!                 1983.
!               3. D. E. Amos, A Subroutine Package for Bessel Functions
!                 of a Complex Argument and Nonnegative Order, Report
!                 SAND85-1018, Sandia National Laboratory, Albuquerque,
!                 NM, May 1985.
!               4. D. E. Amos, A portable package for Bessel functions
!                 of a complex argument and nonnegative order, ACM
!                 Transactions on Mathematical Software, 12 (September
!                 1986), pp. 265-273.
!
!***ROUTINES CALLED  CACAI, CBKNU, I1MACH, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   890801  REVISION DATE from Version 3.2
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   920128  Category corrected.  (WRB)
!   920811  Prologue revised.  (DWL)
!***END PROLOGUE  CAIRY
  COMPLEX AI, CONE, CSQ, CY, S1, S2, TRM1, TRM2, Z, ZTA, Z3
  REAL AA, AD, AK, ALIM, ATRM, AZ, AZ3, BK, CK, COEF, C1, C2, DIG, &
   DK, D1, D2, ELIM, FID, FNU, RL, R1M5, SFAC, TOL, TTH, ZI, ZR, &
   Z3I, Z3R, R1MACH, BB, ALAZ
  INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, I1MACH
  DIMENSION CY(1)
  DATA TTH, C1, C2, COEF /6.66666666666666667E-01, &
   3.55028053887817240E-01,2.58819403792806799E-01, &
   1.83776298473930683E-01/
  DATA  CONE / (1.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CAIRY
  IERR = 0
  NZ=0
  if (ID < 0 .OR. ID > 1) IERR=1
  if (KODE < 1 .OR. KODE > 2) IERR=1
  if (IERR /= 0) RETURN
  AZ = ABS(Z)
  TOL = MAX(R1MACH(4),1.0E-18)
  FID = ID
  if (AZ > 1.0E0) go to 60
!-----------------------------------------------------------------------
!     POWER SERIES FOR ABS(Z) <= 1.
!-----------------------------------------------------------------------
  S1 = CONE
  S2 = CONE
  if (AZ < TOL) go to 160
  AA = AZ*AZ
  if (AA < TOL/AZ) go to 40
  TRM1 = CONE
  TRM2 = CONE
  ATRM = 1.0E0
  Z3 = Z*Z*Z
  AZ3 = AZ*AA
  AK = 2.0E0 + FID
  BK = 3.0E0 - FID - FID
  CK = 4.0E0 - FID
  DK = 3.0E0 + FID + FID
  D1 = AK*DK
  D2 = BK*CK
  AD = MIN(D1,D2)
  AK = 24.0E0 + 9.0E0*FID
  BK = 30.0E0 - 9.0E0*FID
  Z3R = REAL(Z3)
  Z3I = AIMAG(Z3)
  DO 30 K=1,25
    TRM1 = TRM1*CMPLX(Z3R/D1,Z3I/D1)
    S1 = S1 + TRM1
    TRM2 = TRM2*CMPLX(Z3R/D2,Z3I/D2)
    S2 = S2 + TRM2
    ATRM = ATRM*AZ3/AD
    D1 = D1 + AK
    D2 = D2 + BK
    AD = MIN(D1,D2)
    if (ATRM < TOL*AD) go to 40
    AK = AK + 18.0E0
    BK = BK + 18.0E0
   30 CONTINUE
   40 CONTINUE
  if (ID == 1) go to 50
  AI = S1*CMPLX(C1,0.0E0) - Z*S2*CMPLX(C2,0.0E0)
  if (KODE == 1) RETURN
  ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
  AI = AI*CEXP(ZTA)
  return
   50 CONTINUE
  AI = -S2*CMPLX(C2,0.0E0)
  if (AZ > TOL) AI = AI + Z*Z*S1*CMPLX(C1/(1.0E0+FID),0.0E0)
  if (KODE == 1) RETURN
  ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
  AI = AI*CEXP(ZTA)
  return
!-----------------------------------------------------------------------
!     CASE FOR ABS(Z) > 1.0
!-----------------------------------------------------------------------
   60 CONTINUE
  FNU = (1.0E0+FID)/3.0E0
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!-----------------------------------------------------------------------
  K1 = I1MACH(12)
  K2 = I1MACH(13)
  R1M5 = R1MACH(5)
  K = MIN(ABS(K1),ABS(K2))
  ELIM = 2.303E0*(K*R1M5-3.0E0)
  K1 = I1MACH(11) - 1
  AA = R1M5*K1
  DIG = MIN(AA,18.0E0)
  AA = AA*2.303E0
  ALIM = ELIM + MAX(-AA,-41.45E0)
  RL = 1.2E0*DIG + 3.0E0
  ALAZ=ALOG(AZ)
!-----------------------------------------------------------------------
!     TEST FOR RANGE
!-----------------------------------------------------------------------
  AA=0.5E0/TOL
  BB=I1MACH(9)*0.5E0
  AA=MIN(AA,BB)
  AA=AA**TTH
  if (AZ > AA) go to 260
  AA=SQRT(AA)
  if (AZ > AA) IERR=3
  CSQ=CSQRT(Z)
  ZTA=Z*CSQ*CMPLX(TTH,0.0E0)
!-----------------------------------------------------------------------
!     RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
!-----------------------------------------------------------------------
  IFLAG = 0
  SFAC = 1.0E0
  ZI = AIMAG(Z)
  ZR = REAL(Z)
  AK = AIMAG(ZTA)
  if (ZR >= 0.0E0) go to 70
  BK = REAL(ZTA)
  CK = -ABS(BK)
  ZTA = CMPLX(CK,AK)
   70 CONTINUE
  if (ZI /= 0.0E0) go to 80
  if (ZR > 0.0E0) go to 80
  ZTA = CMPLX(0.0E0,AK)
   80 CONTINUE
  AA = REAL(ZTA)
  if (AA >= 0.0E0 .AND. ZR > 0.0E0) go to 100
  if (KODE == 2) go to 90
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
  if (AA > (-ALIM)) go to 90
  AA = -AA + 0.25E0*ALAZ
  IFLAG = 1
  SFAC = TOL
  if (AA > ELIM) go to 240
   90 CONTINUE
!-----------------------------------------------------------------------
!     CBKNU AND CACAI RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
!-----------------------------------------------------------------------
  MR = 1
  if (ZI < 0.0E0) MR = -1
  call CACAI(ZTA, FNU, KODE, MR, 1, CY, NN, RL, TOL, ELIM, ALIM)
  if (NN < 0) go to 250
  NZ = NZ + NN
  go to 120
  100 CONTINUE
  if (KODE == 2) go to 110
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
  if (AA < ALIM) go to 110
  AA = -AA - 0.25E0*ALAZ
  IFLAG = 2
  SFAC = 1.0E0/TOL
  if (AA < (-ELIM)) go to 180
  110 CONTINUE
  call CBKNU(ZTA, FNU, KODE, 1, CY, NZ, TOL, ELIM, ALIM)
  120 CONTINUE
  S1 = CY(1)*CMPLX(COEF,0.0E0)
  if (IFLAG /= 0) go to 140
  if (ID == 1) go to 130
  AI = CSQ*S1
  return
  130 AI = -Z*S1
  return
  140 CONTINUE
  S1 = S1*CMPLX(SFAC,0.0E0)
  if (ID == 1) go to 150
  S1 = S1*CSQ
  AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
  return
  150 CONTINUE
  S1 = -S1*Z
  AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
  return
  160 CONTINUE
  AA = 1.0E+3*R1MACH(1)
  S1 = CMPLX(0.0E0,0.0E0)
  if (ID == 1) go to 170
  if (AZ > AA) S1 = CMPLX(C2,0.0E0)*Z
  AI = CMPLX(C1,0.0E0) - S1
  return
  170 CONTINUE
  AI = -CMPLX(C2,0.0E0)
  AA = SQRT(AA)
  if (AZ > AA) S1 = Z*Z*CMPLX(0.5E0,0.0E0)
  AI = AI + S1*CMPLX(C1,0.0E0)
  return
  180 CONTINUE
  NZ = 1
  AI = CMPLX(0.0E0,0.0E0)
  return
  240 CONTINUE
  NZ = 0
  IERR=2
  return
  250 CONTINUE
  if ( NN == (-1)) go to 240
  NZ=0
  IERR=5
  return
  260 CONTINUE
  IERR=4
  NZ=0
  return
end
