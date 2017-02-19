subroutine CBIRY (Z, ID, KODE, BI, IERR)
!
!! CBIRY computes the Airy function Bi(z) or its derivative dBi/dz ...
!            for complex argument z.  A scaling option is available ...
!            to help avoid overflow.
!
!***LIBRARY   SLATEC
!***CATEGORY  C10D
!***TYPE      COMPLEX (CBIRY-C, ZBIRY-C)
!***KEYWORDS  AIRY FUNCTION, BESSEL FUNCTION OF ORDER ONE THIRD,
!             BESSEL FUNCTION OF ORDER TWO THIRDS
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!         On KODE=1, CBIRY computes the complex Airy function Bi(z)
!         or its derivative dBi/dz on ID=0 or ID=1 respectively.
!         On KODE=2, a scaling option exp(abs(Re(zeta)))*Bi(z) or
!         exp(abs(Re(zeta)))*dBi/dz is provided to remove the
!         exponential behavior in both the left and right half planes
!         where zeta=(2/3)*z**(3/2).
!
!         The Airy functions Bi(z) and dBi/dz are analytic in the
!         whole z-plane, and the scaling option does not destroy this
!         property.
!
!         Input
!           Z      - Argument of type COMPLEX
!           ID     - Order of derivative, ID=0 or ID=1
!           KODE   - A parameter to indicate the scaling option
!                    KODE=1  returns
!                            BI=Bi(z)  on ID=0
!                            BI=dBi/dz on ID=1
!                            at z=Z
!                        =2  returns
!                            BI=exp(abs(Re(zeta)))*Bi(z)  on ID=0
!                            BI=exp(abs(Re(zeta)))*dBi/dz on ID=1
!                            at z=Z where zeta=(2/3)*z**(3/2)
!
!         Output
!           BI     - Result of type COMPLEX
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
!         Bi(z) and dBi/dz are computed from I Bessel functions by
!
!                Bi(z) =  c*sqrt(z)*( I(-1/3,zeta) + I(1/3,zeta) )
!               dBi/dz =  c*   z   *( I(-2/3,zeta) + I(2/3,zeta) )
!                    c =  1/sqrt(3)
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
!***ROUTINES CALLED  CBINU, I1MACH, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   890801  REVISION DATE from Version 3.2
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   920128  Category corrected.  (WRB)
!   920811  Prologue revised.  (DWL)
!***END PROLOGUE  CBIRY
  COMPLEX BI, CONE, CSQ, CY, S1, S2, TRM1, TRM2, Z, ZTA, Z3
  REAL AA, AD, AK, ALIM, ATRM, AZ, AZ3, BB, BK, CK, COEF, C1, C2, &
   DIG, DK, D1, D2, ELIM, FID, FMR, FNU, FNUL, PI, RL, R1M5, SFAC, &
   TOL, TTH, ZI, ZR, Z3I, Z3R, R1MACH
  INTEGER ID, IERR, K, KODE, K1, K2, NZ, I1MACH
  DIMENSION CY(2)
  DATA TTH, C1, C2, COEF, PI /6.66666666666666667E-01, &
   6.14926627446000736E-01,4.48288357353826359E-01, &
   5.77350269189625765E-01,3.14159265358979324E+00/
  DATA  CONE / (1.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CBIRY
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
  if (AZ < TOL) go to 110
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
  BI = S1*CMPLX(C1,0.0E0) + Z*S2*CMPLX(C2,0.0E0)
  if (KODE == 1) RETURN
  ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
  AA = REAL(ZTA)
  AA = -ABS(AA)
  BI = BI*CMPLX(EXP(AA),0.0E0)
  return
   50 CONTINUE
  BI = S2*CMPLX(C2,0.0E0)
  if (AZ > TOL) BI = BI + Z*Z*S1*CMPLX(C1/(1.0E0+FID),0.0E0)
  if (KODE == 1) RETURN
  ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
  AA = REAL(ZTA)
  AA = -ABS(AA)
  BI = BI*CMPLX(EXP(AA),0.0E0)
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
!     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
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
  FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
!-----------------------------------------------------------------------
!     TEST FOR RANGE
!-----------------------------------------------------------------------
  AA=0.5E0/TOL
  BB=I1MACH(9)*0.5E0
  AA=MIN(AA,BB)
  AA=AA**TTH
  if (AZ > AA) go to 190
  AA=SQRT(AA)
  if (AZ > AA) IERR=3
  CSQ=CSQRT(Z)
  ZTA=Z*CSQ*CMPLX(TTH,0.0E0)
!-----------------------------------------------------------------------
!     RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
!-----------------------------------------------------------------------
  SFAC = 1.0E0
  ZI = AIMAG(Z)
  ZR = REAL(Z)
  AK = AIMAG(ZTA)
  if (ZR >= 0.0E0) go to 70
  BK = REAL(ZTA)
  CK = -ABS(BK)
  ZTA = CMPLX(CK,AK)
   70 CONTINUE
  if (ZI == 0.0E0 .AND. ZR <= 0.0E0) ZTA = CMPLX(0.0E0,AK)
  AA = REAL(ZTA)
  if (KODE == 2) go to 80
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
  BB = ABS(AA)
  if (BB < ALIM) go to 80
  BB = BB + 0.25E0*ALOG(AZ)
  SFAC = TOL
  if (BB > ELIM) go to 170
   80 CONTINUE
  FMR = 0.0E0
  if (AA >= 0.0E0 .AND. ZR > 0.0E0) go to 90
  FMR = PI
  if (ZI < 0.0E0) FMR = -PI
  ZTA = -ZTA
   90 CONTINUE
!-----------------------------------------------------------------------
!     AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA)
!     KODE=2 RETURNS EXP(-ABS(XZTA))*I(FNU,ZTA) FROM CBINU
!-----------------------------------------------------------------------
  call CBINU(ZTA, FNU, KODE, 1, CY, NZ, RL, FNUL, TOL, ELIM, ALIM)
  if (NZ < 0) go to 180
  AA = FMR*FNU
  Z3 = CMPLX(SFAC,0.0E0)
  S1 = CY(1)*CMPLX(COS(AA),SIN(AA))*Z3
  FNU = (2.0E0-FID)/3.0E0
  call CBINU(ZTA, FNU, KODE, 2, CY, NZ, RL, FNUL, TOL, ELIM, ALIM)
  CY(1) = CY(1)*Z3
  CY(2) = CY(2)*Z3
!-----------------------------------------------------------------------
!     BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3
!-----------------------------------------------------------------------
  S2 = CY(1)*CMPLX(FNU+FNU,0.0E0)/ZTA + CY(2)
  AA = FMR*(FNU-1.0E0)
  S1 = (S1+S2*CMPLX(COS(AA),SIN(AA)))*CMPLX(COEF,0.0E0)
  if (ID == 1) go to 100
  S1 = CSQ*S1
  BI = S1*CMPLX(1.0E0/SFAC,0.0E0)
  return
  100 CONTINUE
  S1 = Z*S1
  BI = S1*CMPLX(1.0E0/SFAC,0.0E0)
  return
  110 CONTINUE
  AA = C1*(1.0E0-FID) + FID*C2
  BI = CMPLX(AA,0.0E0)
  return
  170 CONTINUE
  NZ=0
  IERR=2
  return
  180 CONTINUE
  if ( NZ == (-1)) go to 170
  NZ=0
  IERR=5
  return
  190 CONTINUE
  IERR=4
  NZ=0
  return
end
