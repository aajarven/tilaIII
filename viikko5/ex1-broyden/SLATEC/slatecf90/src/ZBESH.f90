subroutine ZBESH (ZR, ZI, FNU, KODE, M, N, CYR, CYI, NZ, IERR)
!
!! ZBESH computes a sequence of the Hankel functions H(m,a,z) ...
!            for superscript m=1 or 2, real nonnegative orders a=b, ...
!            b+1,... where b>0, and nonzero complex argument z.  A ...
!            scaling option is available to help avoid overflow.
!
!***LIBRARY   SLATEC
!***CATEGORY  C10A4
!***TYPE      COMPLEX (CBESH-C, ZBESH-C)
!***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
!             BESSEL FUNCTIONS OF THE THIRD KIND, H BESSEL FUNCTIONS,
!             HANKEL FUNCTIONS
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!                      ***A DOUBLE PRECISION ROUTINE***
!         On KODE=1, ZBESH computes an N member sequence of complex
!         Hankel (Bessel) functions CY(L)=H(M,FNU+L-1,Z) for super-
!         script M=1 or 2, real nonnegative orders FNU+L-1, L=1,...,
!         N, and complex nonzero Z in the cut plane -pi<arg(Z)<=pi
!         where Z=ZR+i*ZI.  On KODE=2, CBESH returns the scaled
!         functions
!
!            CY(L) = H(M,FNU+L-1,Z)*exp(-(3-2*M)*Z*i),  i**2=-1
!
!         which removes the exponential behavior in both the upper
!         and lower half planes.  Definitions and notation are found
!         in the NBS Handbook of Mathematical Functions (Ref. 1).
!
!         Input
!           ZR     - DOUBLE PRECISION real part of nonzero argument Z
!           ZI     - DOUBLE PRECISION imag part of nonzero argument Z
!           FNU    - DOUBLE PRECISION initial order, FNU>=0
!           KODE   - A parameter to indicate the scaling option
!                    KODE=1  returns
!                            CY(L)=H(M,FNU+L-1,Z), L=1,...,N
!                        =2  returns
!                            CY(L)=H(M,FNU+L-1,Z)*exp(-(3-2M)*Z*i),
!                            L=1,...,N
!           M      - Superscript of Hankel function, M=1 or 2
!           N      - Number of terms in the sequence, N>=1
!
!         Output
!           CYR    - DOUBLE PRECISION real part of result vector
!           CYI    - DOUBLE PRECISION imag part of result vector
!           NZ     - Number of underflows set to zero
!                    NZ=0    Normal return
!                    NZ>0    CY(L)=0 for NZ values of L (if M=1 and
!                            Im(Z)>0 or if M=2 and Im(Z)<0, then
!                            CY(L)=0 for L=1,...,NZ; in the com-
!                            plementary half planes, the underflows
!                            may not be in an uninterrupted sequence)
!           IERR   - Error flag
!                    IERR=0  Normal return     - COMPUTATION COMPLETED
!                    IERR=1  Input error       - NO COMPUTATION
!                    IERR=2  Overflow          - NO COMPUTATION
!                            (abs(Z) too small and/or FNU+N-1
!                            too large)
!                    IERR=3  Precision warning - COMPUTATION COMPLETED
!                            (Result has half precision or less
!                            because abs(Z) or FNU+N-1 is large)
!                    IERR=4  Precision error   - NO COMPUTATION
!                            (Result has no precision because
!                            abs(Z) or FNU+N-1 is too large)
!                    IERR=5  Algorithmic error - NO COMPUTATION
!                            (Termination condition not met)
!
! *Long Description:
!
!         The computation is carried out by the formula
!
!            H(m,a,z) = (1/t)*exp(-a*t)*K(a,z*exp(-t))
!                   t = (3-2*m)*i*pi/2
!
!         where the K Bessel function is computed as described in the
!         prologue to CBESK.
!
!         Exponential decay of H(m,a,z) occurs in the upper half z
!         plane for m=1 and the lower half z plane for m=2.  Exponential
!         growth occurs in the complementary half planes.  Scaling
!         by exp(-(3-2*m)*z*i) removes the exponential behavior in the
!         whole z plane as z goes to infinity.
!
!         For negative orders, the formula
!
!            H(m,-a,z) = H(m,a,z)*exp((3-2*m)*a*pi*i)
!
!         can be used.
!
!         In most complex variable computation, one must evaluate ele-
!         mentary functions.  When the magnitude of Z or FNU+N-1 is
!         large, losses of significance by argument reduction occur.
!         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
!         losses exceeding half precision are likely and an error flag
!         IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is double
!         precision unit roundoff limited to 18 digits precision.  Also,
!         if either is larger than U2=0.5/UR, then all significance is
!         lost and IERR=4.  In order to use the INT function, arguments
!         must be further restricted not to exceed the largest machine
!         integer, U3=I1MACH(9).  Thus, the magnitude of Z and FNU+N-1
!         is restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, and
!         U3 approximate 2.0E+3, 4.2E+6, 2.1E+9 in single precision
!         and 4.7E+7, 2.3E+15 and 2.1E+9 in double precision.  This
!         makes U2 limiting in single precision and U3 limiting in
!         double precision.  This means that one can expect to retain,
!         in the worst cases on IEEE machines, no digits in single pre-
!         cision and only 6 digits in double precision.  Similar con-
!         siderations hold for other machines.
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
!         magnitude of the larger component.  In these extreme cases,
!         the principal phase angle is on the order of +P, -P, PI/2-P,
!         or -PI/2+P.
!
!***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
!                 matical Functions, National Bureau of Standards
!                 Applied Mathematics Series 55, U. S. Department
!                 of Commerce, Tenth Printing (1972) or later.
!               2. D. E. Amos, Computation of Bessel Functions of
!                 Complex Argument, Report SAND83-0086, Sandia National
!                 Laboratories, Albuquerque, NM, May 1983.
!               3. D. E. Amos, Computation of Bessel Functions of
!                 Complex Argument and Large Order, Report SAND83-0643,
!                 Sandia National Laboratories, Albuquerque, NM, May
!                 1983.
!               4. D. E. Amos, A Subroutine Package for Bessel Functions
!                 of a Complex Argument and Nonnegative Order, Report
!                 SAND85-1018, Sandia National Laboratory, Albuquerque,
!                 NM, May 1985.
!               5. D. E. Amos, A portable package for Bessel functions
!                 of a complex argument and nonnegative order, ACM
!                 Transactions on Mathematical Software, 12 (September
!                 1986), pp. 265-273.
!
!***ROUTINES CALLED  D1MACH, I1MACH, ZABS, ZACON, ZBKNU, ZBUNK, ZUOIK
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   890801  REVISION DATE from Version 3.2
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   920128  Category corrected.  (WRB)
!   920811  Prologue revised.  (DWL)
!***END PROLOGUE  ZBESH
!
!     COMPLEX CY,Z,ZN,ZT,CSGN
  DOUBLE PRECISION AA, ALIM, ALN, ARG, AZ, CYI, CYR, DIG, ELIM, &
   FMM, FN, FNU, FNUL, HPI, RHPI, RL, R1M5, SGN, STR, TOL, UFL, ZI, &
   ZNI, ZNR, ZR, ZTI, D1MACH, ZABS, BB, ASCLE, RTOL, ATOL, STI, &
   CSGNR, CSGNI
  INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, M, &
   MM, MR, N, NN, NUF, NW, NZ, I1MACH
  DIMENSION CYR(N), CYI(N)
  EXTERNAL ZABS
!
  DATA HPI /1.57079632679489662D0/
!
!***FIRST EXECUTABLE STATEMENT  ZBESH
  IERR = 0
  NZ=0
  if (ZR == 0.0D0 .AND. ZI == 0.0D0) IERR=1
  if (FNU < 0.0D0) IERR=1
  if (M < 1 .OR. M > 2) IERR=1
  if (KODE < 1 .OR. KODE > 2) IERR=1
  if (N < 1) IERR=1
  if (IERR /= 0) RETURN
  NN = N
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
!-----------------------------------------------------------------------
  TOL = MAX(D1MACH(4),1.0D-18)
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
  FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
  RL = 1.2D0*DIG + 3.0D0
  FN = FNU + (NN-1)
  MM = 3 - M - M
  FMM = MM
  ZNR = FMM*ZI
  ZNI = -FMM*ZR
!-----------------------------------------------------------------------
!     TEST FOR PROPER RANGE
!-----------------------------------------------------------------------
  AZ = ZABS(ZR,ZI)
  AA = 0.5D0/TOL
  BB = I1MACH(9)*0.5D0
  AA = MIN(AA,BB)
  if (AZ > AA) go to 260
  if (FN > AA) go to 260
  AA = SQRT(AA)
  if (AZ > AA) IERR=3
  if (FN > AA) IERR=3
!-----------------------------------------------------------------------
!     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
!-----------------------------------------------------------------------
  UFL = D1MACH(1)*1.0D+3
  if (AZ < UFL) go to 230
  if (FNU > FNUL) go to 90
  if (FN <= 1.0D0) go to 70
  if (FN > 2.0D0) go to 60
  if (AZ > TOL) go to 70
  ARG = 0.5D0*AZ
  ALN = -FN*LOG(ARG)
  if (ALN > ELIM) go to 230
  go to 70
   60 CONTINUE
  call ZUOIK(ZNR, ZNI, FNU, KODE, 2, NN, CYR, CYI, NUF, TOL, ELIM, &
   ALIM)
  if (NUF < 0) go to 230
  NZ = NZ + NUF
  NN = NN - NUF
!-----------------------------------------------------------------------
!     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
!     if NUF=NN, THEN CY(I)=CZERO FOR ALL I
!-----------------------------------------------------------------------
  if (NN == 0) go to 140
   70 CONTINUE
  if ((ZNR < 0.0D0) .OR. (ZNR == 0.0D0 .AND. ZNI < 0.0D0 .AND. &
   M == 2)) go to 80
!-----------------------------------------------------------------------
!     RIGHT HALF PLANE COMPUTATION, XN >= 0. .AND. (XN /= 0. .OR.
!     YN >= 0. .OR. M=1)
!-----------------------------------------------------------------------
  call ZBKNU(ZNR, ZNI, FNU, KODE, NN, CYR, CYI, NZ, TOL, ELIM, ALIM)
  go to 110
!-----------------------------------------------------------------------
!     LEFT HALF PLANE COMPUTATION
!-----------------------------------------------------------------------
   80 CONTINUE
  MR = -MM
  call ZACON(ZNR, ZNI, FNU, KODE, MR, NN, CYR, CYI, NW, RL, FNUL, &
   TOL, ELIM, ALIM)
  if (NW < 0) go to 240
  NZ=NW
  go to 110
   90 CONTINUE
!-----------------------------------------------------------------------
!     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU > FNUL
!-----------------------------------------------------------------------
  MR = 0
  if ((ZNR >= 0.0D0) .AND. (ZNR /= 0.0D0 .OR. ZNI >= 0.0D0 .OR. &
   M /= 2)) go to 100
  MR = -MM
  if (ZNR /= 0.0D0 .OR. ZNI >= 0.0D0) go to 100
  ZNR = -ZNR
  ZNI = -ZNI
  100 CONTINUE
  call ZBUNK(ZNR, ZNI, FNU, KODE, MR, NN, CYR, CYI, NW, TOL, ELIM, &
   ALIM)
  if (NW < 0) go to 240
  NZ = NZ + NW
  110 CONTINUE
!-----------------------------------------------------------------------
!     H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)
!
!     ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
!-----------------------------------------------------------------------
  SGN = DSIGN(HPI,-FMM)
!-----------------------------------------------------------------------
!     CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
  INU = FNU
  INUH = INU/2
  IR = INU - 2*INUH
  ARG = (FNU-(INU-IR))*SGN
  RHPI = 1.0D0/SGN
!     ZNI = RHPI*COS(ARG)
!     ZNR = -RHPI*SIN(ARG)
  CSGNI = RHPI*COS(ARG)
  CSGNR = -RHPI*SIN(ARG)
  if (MOD(INUH,2) == 0) go to 120
!     ZNR = -ZNR
!     ZNI = -ZNI
  CSGNR = -CSGNR
  CSGNI = -CSGNI
  120 CONTINUE
  ZTI = -FMM
  RTOL = 1.0D0/TOL
  ASCLE = UFL*RTOL
  DO 130 I=1,NN
!       STR = CYR(I)*ZNR - CYI(I)*ZNI
!       CYI(I) = CYR(I)*ZNI + CYI(I)*ZNR
!       CYR(I) = STR
!       STR = -ZNI*ZTI
!       ZNI = ZNR*ZTI
!       ZNR = STR
    AA = CYR(I)
    BB = CYI(I)
    ATOL = 1.0D0
    if (MAX(ABS(AA),ABS(BB)) > ASCLE) go to 135
      AA = AA*RTOL
      BB = BB*RTOL
      ATOL = TOL
  135 CONTINUE
  STR = AA*CSGNR - BB*CSGNI
  STI = AA*CSGNI + BB*CSGNR
  CYR(I) = STR*ATOL
  CYI(I) = STI*ATOL
  STR = -CSGNI*ZTI
  CSGNI = CSGNR*ZTI
  CSGNR = STR
  130 CONTINUE
  return
  140 CONTINUE
  if (ZNR < 0.0D0) go to 230
  return
  230 CONTINUE
  NZ=0
  IERR=2
  return
  240 CONTINUE
  if ( NW == (-1)) go to 230
  NZ=0
  IERR=5
  return
  260 CONTINUE
  NZ=0
  IERR=4
  return
end
