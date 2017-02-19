subroutine CBESY (Z, FNU, KODE, N, CY, NZ, CWRK, IERR)
!
!! CBESY computes a sequence of the Bessel functions Y(a,z) for ...
!            complex argument z and real nonnegative orders a=b,b+1, ...
!            b+2,... where b>0.  A scaling option is available to ...
!            help avoid overflow.
!
!***LIBRARY   SLATEC
!***CATEGORY  C10A4
!***TYPE      COMPLEX (CBESY-C, ZBESY-C)
!***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
!             BESSEL FUNCTIONS OF SECOND KIND, WEBER'S FUNCTION,
!             Y BESSEL FUNCTIONS
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!         On KODE=1, CBESY computes an N member sequence of complex
!         Bessel functions CY(L)=Y(FNU+L-1,Z) for real nonnegative
!         orders FNU+L-1, L=1,...,N and complex Z in the cut plane
!         -pi<arg(Z)<=pi.  On KODE=2, CBESY returns the scaled
!         functions
!
!            CY(L) = exp(-abs(Y))*Y(FNU+L-1,Z),  L=1,...,N, Y=Im(Z)
!
!         which remove the exponential growth in both the upper and
!         lower half planes as Z goes to infinity.  Definitions and
!         notation are found in the NBS Handbook of Mathematical
!         Functions (Ref. 1).
!
!         Input
!           Z      - Nonzero argument of type COMPLEX
!           FNU    - Initial order of type REAL, FNU>=0
!           KODE   - A parameter to indicate the scaling option
!                    KODE=1  returns
!                            CY(L)=Y(FNU+L-1,Z), L=1,...,N
!                        =2  returns
!                            CY(L)=Y(FNU+L-1,Z)*exp(-abs(Y)), L=1,...,N
!                            where Y=Im(Z)
!           N      - Number of terms in the sequence, N>=1
!           CWRK   - A work vector of type COMPLEX and dimension N
!
!         Output
!           CY     - Result vector of type COMPLEX
!           NZ     - Number of underflows set to zero
!                    NZ=0    Normal return
!                    NZ>0    CY(L)=0 for NZ values of L, usually on
!                            KODE=2 (the underflows may not be in an
!                            uninterrupted sequence)
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
!            Y(a,z) = (H(1,a,z) - H(2,a,z))/(2*i)
!
!         where the Hankel functions are computed as described in CBESH.
!
!         For negative orders, the formula
!
!            Y(-a,z) = Y(a,z)*cos(a*pi) + J(a,z)*sin(a*pi)
!
!         can be used.  However, for large orders close to half odd
!         integers the function changes radically.  When a is a large
!         positive half odd integer, the magnitude of Y(-a,z)=J(a,z)*
!         sin(a*pi) is a large negative power of ten.  But when a is
!         not a half odd integer, Y(a,z) dominates in magnitude with a
!         large positive power of ten and the most that the second term
!         can be reduced is by unit roundoff from the coefficient.
!         Thus,  wide changes can occur within unit roundoff of a large
!         half odd integer.  Here, large means a>abs(z).
!
!         In most complex variable computation, one must evaluate ele-
!         mentary functions.  When the magnitude of Z or FNU+N-1 is
!         large, losses of significance by argument reduction occur.
!         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
!         losses exceeding half precision are likely and an error flag
!         IERR=3 is triggered where UR=R1MACH(4)=UNIT ROUNDOFF.  Also,
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
!***ROUTINES CALLED  CBESH, I1MACH, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   890801  REVISION DATE from Version 3.2
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   920128  Category corrected.  (WRB)
!   920811  Prologue revised.  (DWL)
!***END PROLOGUE  CBESY
!
  COMPLEX CWRK, CY, C1, C2, EX, HCI, Z, ZU, ZV
  REAL ELIM, EY, FNU, R1, R2, TAY, XX, YY, R1MACH, R1M5, ASCLE, &
    RTOL, ATOL, TOL, AA, BB
  INTEGER I, IERR, K, KODE, K1, K2, N, NZ, NZ1, NZ2, I1MACH
  DIMENSION CY(N), CWRK(N)
!***FIRST EXECUTABLE STATEMENT  CBESY
  XX = REAL(Z)
  YY = AIMAG(Z)
  IERR = 0
  NZ=0
  if (XX == 0.0E0 .AND. YY == 0.0E0) IERR=1
  if (FNU < 0.0E0) IERR=1
  if (KODE < 1 .OR. KODE > 2) IERR=1
  if (N < 1) IERR=1
  if (IERR /= 0) RETURN
  HCI = CMPLX(0.0E0,0.5E0)
  call CBESH(Z, FNU, KODE, 1, N, CY, NZ1, IERR)
  if (IERR /= 0.AND.IERR /= 3) go to 170
  call CBESH(Z, FNU, KODE, 2, N, CWRK, NZ2, IERR)
  if (IERR /= 0.AND.IERR /= 3) go to 170
  NZ = MIN(NZ1,NZ2)
  if (KODE == 2) go to 60
  DO 50 I=1,N
    CY(I) = HCI*(CWRK(I)-CY(I))
   50 CONTINUE
  return
   60 CONTINUE
  TOL = MAX(R1MACH(4),1.0E-18)
  K1 = I1MACH(12)
  K2 = I1MACH(13)
  K = MIN(ABS(K1),ABS(K2))
  R1M5 = R1MACH(5)
!-----------------------------------------------------------------------
!     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT
!-----------------------------------------------------------------------
  ELIM = 2.303E0*(K*R1M5-3.0E0)
  R1 = COS(XX)
  R2 = SIN(XX)
  EX = CMPLX(R1,R2)
  EY = 0.0E0
  TAY = ABS(YY+YY)
  if (TAY < ELIM) EY = EXP(-TAY)
  if (YY < 0.0E0) go to 90
  C1 = EX*CMPLX(EY,0.0E0)
  C2 = CONJG(EX)
   70 CONTINUE
  NZ = 0
  RTOL = 1.0E0/TOL
  ASCLE = R1MACH(1)*RTOL*1.0E+3
  DO 80 I=1,N
!       CY(I) = HCI*(C2*CWRK(I)-C1*CY(I))
    ZV = CWRK(I)
    AA=REAL(ZV)
    BB=AIMAG(ZV)
    ATOL=1.0E0
    if (MAX(ABS(AA),ABS(BB)) > ASCLE) go to 75
      ZV = ZV*CMPLX(RTOL,0.0E0)
      ATOL = TOL
   75   CONTINUE
    ZV = ZV*C2*HCI
    ZV = ZV*CMPLX(ATOL,0.0E0)
    ZU=CY(I)
    AA=REAL(ZU)
    BB=AIMAG(ZU)
    ATOL=1.0E0
    if (MAX(ABS(AA),ABS(BB)) > ASCLE) go to 85
      ZU = ZU*CMPLX(RTOL,0.0E0)
      ATOL = TOL
   85   CONTINUE
    ZU = ZU*C1*HCI
    ZU = ZU*CMPLX(ATOL,0.0E0)
    CY(I) = ZV - ZU
    if (CY(I) == CMPLX(0.0E0,0.0E0) .AND. EY == 0.0E0) NZ = NZ + 1
   80 CONTINUE
  return
   90 CONTINUE
  C1 = EX
  C2 = CONJG(EX)*CMPLX(EY,0.0E0)
  go to 70
  170 CONTINUE
  NZ = 0
  return
end
