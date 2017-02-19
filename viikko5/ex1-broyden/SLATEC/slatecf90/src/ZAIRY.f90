subroutine ZAIRY (ZR, ZI, ID, KODE, AIR, AII, NZ, IERR)
!
!! ZAIRY computes the Airy function Ai(z) or its derivative dAi/dz ...
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
!                      ***A DOUBLE PRECISION ROUTINE***
!         On KODE=1, ZAIRY computes the complex Airy function Ai(z)
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
!           ZR     - DOUBLE PRECISION real part of argument Z
!           ZI     - DOUBLE PRECISION imag part of argument Z
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
!           AIR    - DOUBLE PRECISION real part of result
!           AII    - DOUBLE PRECISION imag part of result
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
!         flag IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is
!         double precision unit roundoff limited to 18 digits precision.
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
!***ROUTINES CALLED  D1MACH, I1MACH, ZABS, ZACAI, ZBKNU, ZEXP, ZSQRT
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   890801  REVISION DATE from Version 3.2
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   920128  Category corrected.  (WRB)
!   920811  Prologue revised.  (DWL)
!   930122  Added ZEXP and ZSQRT to EXTERNAL statement.  (RWC)
!***END PROLOGUE  ZAIRY
!     COMPLEX AI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3
  DOUBLE PRECISION AA, AD, AII, AIR, AK, ALIM, ATRM, AZ, AZ3, BK, &
   CC, CK, COEF, CONEI, CONER, CSQI, CSQR, CYI, CYR, C1, C2, DIG, &
   DK, D1, D2, ELIM, FID, FNU, PTR, RL, R1M5, SFAC, STI, STR, &
   S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I, TRM2R, TTH, ZEROI, &
   ZEROR, ZI, ZR, ZTAI, ZTAR, Z3I, Z3R, D1MACH, ZABS, ALAZ, BB
  INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, I1MACH
  DIMENSION CYR(1), CYI(1)
  EXTERNAL ZABS, ZEXP, ZSQRT
  DATA TTH, C1, C2, COEF /6.66666666666666667D-01, &
   3.55028053887817240D-01,2.58819403792806799D-01, &
   1.83776298473930683D-01/
  DATA ZEROR, ZEROI, CONER, CONEI /0.0D0,0.0D0,1.0D0,0.0D0/
!***FIRST EXECUTABLE STATEMENT  ZAIRY
  IERR = 0
  NZ=0
  if (ID < 0 .OR. ID > 1) IERR=1
  if (KODE < 1 .OR. KODE > 2) IERR=1
  if (IERR /= 0) RETURN
  AZ = ZABS(ZR,ZI)
  TOL = MAX(D1MACH(4),1.0D-18)
  FID = ID
  if (AZ > 1.0D0) go to 70
!-----------------------------------------------------------------------
!     POWER SERIES FOR ABS(Z) <= 1.
!-----------------------------------------------------------------------
  S1R = CONER
  S1I = CONEI
  S2R = CONER
  S2I = CONEI
  if (AZ < TOL) go to 170
  AA = AZ*AZ
  if (AA < TOL/AZ) go to 40
  TRM1R = CONER
  TRM1I = CONEI
  TRM2R = CONER
  TRM2I = CONEI
  ATRM = 1.0D0
  STR = ZR*ZR - ZI*ZI
  STI = ZR*ZI + ZI*ZR
  Z3R = STR*ZR - STI*ZI
  Z3I = STR*ZI + STI*ZR
  AZ3 = AZ*AA
  AK = 2.0D0 + FID
  BK = 3.0D0 - FID - FID
  CK = 4.0D0 - FID
  DK = 3.0D0 + FID + FID
  D1 = AK*DK
  D2 = BK*CK
  AD = MIN(D1,D2)
  AK = 24.0D0 + 9.0D0*FID
  BK = 30.0D0 - 9.0D0*FID
  DO 30 K=1,25
    STR = (TRM1R*Z3R-TRM1I*Z3I)/D1
    TRM1I = (TRM1R*Z3I+TRM1I*Z3R)/D1
    TRM1R = STR
    S1R = S1R + TRM1R
    S1I = S1I + TRM1I
    STR = (TRM2R*Z3R-TRM2I*Z3I)/D2
    TRM2I = (TRM2R*Z3I+TRM2I*Z3R)/D2
    TRM2R = STR
    S2R = S2R + TRM2R
    S2I = S2I + TRM2I
    ATRM = ATRM*AZ3/AD
    D1 = D1 + AK
    D2 = D2 + BK
    AD = MIN(D1,D2)
    if (ATRM < TOL*AD) go to 40
    AK = AK + 18.0D0
    BK = BK + 18.0D0
   30 CONTINUE
   40 CONTINUE
  if (ID == 1) go to 50
  AIR = S1R*C1 - C2*(ZR*S2R-ZI*S2I)
  AII = S1I*C1 - C2*(ZR*S2I+ZI*S2R)
  if (KODE == 1) RETURN
  call ZSQRT(ZR, ZI, STR, STI)
  ZTAR = TTH*(ZR*STR-ZI*STI)
  ZTAI = TTH*(ZR*STI+ZI*STR)
  call ZEXP(ZTAR, ZTAI, STR, STI)
  PTR = AIR*STR - AII*STI
  AII = AIR*STI + AII*STR
  AIR = PTR
  return
   50 CONTINUE
  AIR = -S2R*C2
  AII = -S2I*C2
  if (AZ <= TOL) go to 60
  STR = ZR*S1R - ZI*S1I
  STI = ZR*S1I + ZI*S1R
  CC = C1/(1.0D0+FID)
  AIR = AIR + CC*(STR*ZR-STI*ZI)
  AII = AII + CC*(STR*ZI+STI*ZR)
   60 CONTINUE
  if (KODE == 1) RETURN
  call ZSQRT(ZR, ZI, STR, STI)
  ZTAR = TTH*(ZR*STR-ZI*STI)
  ZTAI = TTH*(ZR*STI+ZI*STR)
  call ZEXP(ZTAR, ZTAI, STR, STI)
  PTR = STR*AIR - STI*AII
  AII = STR*AII + STI*AIR
  AIR = PTR
  return
!-----------------------------------------------------------------------
!     CASE FOR ABS(Z) > 1.0
!-----------------------------------------------------------------------
   70 CONTINUE
  FNU = (1.0D0+FID)/3.0D0
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0D-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!-----------------------------------------------------------------------
  K1 = I1MACH(15)
  K2 = I1MACH(16)
  R1M5 = D1MACH(5)
  K = MIN(ABS(K1),ABS(K2))
  ELIM = 2.303D0*(K*R1M5-3.0D0)
  K1 = I1MACH(14) - 1
  AA = R1M5*K1
  DIG = MIN(AA,18.0D0)
  AA = AA*2.303D0
  ALIM = ELIM + MAX(-AA,-41.45D0)
  RL = 1.2D0*DIG + 3.0D0
  ALAZ = LOG(AZ)
!-----------------------------------------------------------------------
!     TEST FOR PROPER RANGE
!-----------------------------------------------------------------------
  AA=0.5D0/TOL
  BB=I1MACH(9)*0.5D0
  AA=MIN(AA,BB)
  AA=AA**TTH
  if (AZ > AA) go to 260
  AA=SQRT(AA)
  if (AZ > AA) IERR=3
  call ZSQRT(ZR, ZI, CSQR, CSQI)
  ZTAR = TTH*(ZR*CSQR-ZI*CSQI)
  ZTAI = TTH*(ZR*CSQI+ZI*CSQR)
!-----------------------------------------------------------------------
!     RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
!-----------------------------------------------------------------------
  IFLAG = 0
  SFAC = 1.0D0
  AK = ZTAI
  if (ZR >= 0.0D0) go to 80
  BK = ZTAR
  CK = -ABS(BK)
  ZTAR = CK
  ZTAI = AK
   80 CONTINUE
  if (ZI /= 0.0D0) go to 90
  if (ZR > 0.0D0) go to 90
  ZTAR = 0.0D0
  ZTAI = AK
   90 CONTINUE
  AA = ZTAR
  if (AA >= 0.0D0 .AND. ZR > 0.0D0) go to 110
  if (KODE == 2) go to 100
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
  if (AA > (-ALIM)) go to 100
  AA = -AA + 0.25D0*ALAZ
  IFLAG = 1
  SFAC = TOL
  if (AA > ELIM) go to 270
  100 CONTINUE
!-----------------------------------------------------------------------
!     CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
!-----------------------------------------------------------------------
  MR = 1
  if (ZI < 0.0D0) MR = -1
  call ZACAI(ZTAR, ZTAI, FNU, KODE, MR, 1, CYR, CYI, NN, RL, TOL, &
   ELIM, ALIM)
  if (NN < 0) go to 280
  NZ = NZ + NN
  go to 130
  110 CONTINUE
  if (KODE == 2) go to 120
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
  if (AA < ALIM) go to 120
  AA = -AA - 0.25D0*ALAZ
  IFLAG = 2
  SFAC = 1.0D0/TOL
  if (AA < (-ELIM)) go to 210
  120 CONTINUE
  call ZBKNU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, TOL, ELIM, &
   ALIM)
  130 CONTINUE
  S1R = CYR(1)*COEF
  S1I = CYI(1)*COEF
  if (IFLAG /= 0) go to 150
  if (ID == 1) go to 140
  AIR = CSQR*S1R - CSQI*S1I
  AII = CSQR*S1I + CSQI*S1R
  return
  140 CONTINUE
  AIR = -(ZR*S1R-ZI*S1I)
  AII = -(ZR*S1I+ZI*S1R)
  return
  150 CONTINUE
  S1R = S1R*SFAC
  S1I = S1I*SFAC
  if (ID == 1) go to 160
  STR = S1R*CSQR - S1I*CSQI
  S1I = S1R*CSQI + S1I*CSQR
  S1R = STR
  AIR = S1R/SFAC
  AII = S1I/SFAC
  return
  160 CONTINUE
  STR = -(S1R*ZR-S1I*ZI)
  S1I = -(S1R*ZI+S1I*ZR)
  S1R = STR
  AIR = S1R/SFAC
  AII = S1I/SFAC
  return
  170 CONTINUE
  AA = 1.0D+3*D1MACH(1)
  S1R = ZEROR
  S1I = ZEROI
  if (ID == 1) go to 190
  if (AZ <= AA) go to 180
  S1R = C2*ZR
  S1I = C2*ZI
  180 CONTINUE
  AIR = C1 - S1R
  AII = -S1I
  return
  190 CONTINUE
  AIR = -C2
  AII = 0.0D0
  AA = SQRT(AA)
  if (AZ <= AA) go to 200
  S1R = 0.5D0*(ZR*ZR-ZI*ZI)
  S1I = ZR*ZI
  200 CONTINUE
  AIR = AIR + C1*S1R
  AII = AII + C1*S1I
  return
  210 CONTINUE
  NZ = 1
  AIR = ZEROR
  AII = ZEROI
  return
  270 CONTINUE
  NZ = 0
  IERR=2
  return
  280 CONTINUE
  if ( NN == (-1)) go to 270
  NZ=0
  IERR=5
  return
  260 CONTINUE
  IERR=4
  NZ=0
  return
end
