subroutine CBESI (Z, FNU, KODE, N, CY, NZ, IERR)
!
!! CBESI computes a sequence of the Bessel functions I(a,z) for ...
!            complex argument z and real nonnegative orders a=b,b+1, ...
!            b+2,... where b>0.  A scaling option is available to ...
!            help avoid overflow.
!
!***LIBRARY   SLATEC
!***CATEGORY  C10B4
!***TYPE      COMPLEX (CBESI-C, ZBESI-C)
!***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT, I BESSEL FUNCTIONS,
!             MODIFIED BESSEL FUNCTIONS
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!         On KODE=1, CBESI computes an N-member sequence of complex
!         Bessel functions CY(L)=I(FNU+L-1,Z) for real nonnegative
!         orders FNU+L-1, L=1,...,N and complex Z in the cut plane
!         -pi<arg(Z)<=pi.  On KODE=2, CBESI returns the scaled functions
!
!            CY(L) = exp(-abs(X))*I(FNU+L-1,Z), L=1,...,N and X=Re(Z)
!
!         which removes the exponential growth in both the left and
!         right half-planes as Z goes to infinity.
!
!         Input
!           Z      - Argument of type COMPLEX
!           FNU    - Initial order of type REAL, FNU>=0
!           KODE   - A parameter to indicate the scaling option
!                    KODE=1  returns
!                            CY(L)=I(FNU+L-1,Z), L=1,...,N
!                        =2  returns
!                            CY(L)=exp(-abs(X))*I(FNU+L-1,Z), L=1,...,N
!                            where X=Re(Z)
!           N      - Number of terms in the sequence, N>=1
!
!         Output
!           CY     - Result vector of type COMPLEX
!           NZ     - Number of underflows set to zero
!                    NZ=0    Normal return
!                    NZ>0    CY(L)=0, L=N-NZ+1,...,N
!           IERR   - Error flag
!                    IERR=0  Normal return     - COMPUTATION COMPLETED
!                    IERR=1  Input error       - NO COMPUTATION
!                    IERR=2  Overflow          - NO COMPUTATION
!                            (Re(Z) too large on KODE=1)
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
!         The computation of I(a,z) is carried out by the power series
!         for small abs(z), the asymptotic expansion for large abs(z),
!         the Miller algorithm normalized by the Wronskian and a
!         Neumann series for intermediate magnitudes of z, and the
!         uniform asymptotic expansions for I(a,z) and J(a,z) for
!         large orders a.  Backward recurrence is used to generate
!         sequences or reduce orders when necessary.
!
!         The calculations above are done in the right half plane and
!         continued into the left half plane by the formula
!
!            I(a,z*exp(t)) = exp(t*a)*I(a,z), Re(z)>0
!                        t = i*pi or -i*pi
!
!         For negative orders, the formula
!
!            I(-a,z) = I(a,z) + (2/pi)*sin(pi*a)*K(a,z)
!
!         can be used.  However, for large orders close to integers the
!         the function changes radically.  When a is a large positive
!         integer, the magnitude of I(-a,z)=I(a,z) is a large
!         negative power of ten. But when a is not an integer,
!         K(a,z) dominates in magnitude with a large positive power of
!         ten and the most that the second term can be reduced is by
!         unit roundoff from the coefficient. Thus, wide changes can
!         occur within unit roundoff of a large integer for a. Here,
!         large means a>abs(z).
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
!***ROUTINES CALLED  CBINU, I1MACH, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   890801  REVISION DATE from Version 3.2
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   920128  Category corrected.  (WRB)
!   920811  Prologue revised.  (DWL)
!***END PROLOGUE  CBESI
  COMPLEX CONE, CSGN, CY, Z, ZN
  REAL AA, ALIM, ARG, DIG, ELIM, FNU, FNUL, PI, RL, R1M5, S1, S2, &
   TOL, XX, YY, R1MACH, AZ, FN, BB, ASCLE, RTOL, ATOL
  INTEGER I, IERR, INU, K, KODE, K1, K2, N, NN, NZ, I1MACH
  DIMENSION CY(N)
  DATA PI /3.14159265358979324E0/
  DATA CONE / (1.0E0,0.0E0) /
!
!***FIRST EXECUTABLE STATEMENT  CBESI
  IERR = 0
  NZ=0
  if (FNU < 0.0E0) IERR=1
  if (KODE < 1 .OR. KODE > 2) IERR=1
  if (N < 1) IERR=1
  if (IERR /= 0) RETURN
  XX = REAL(Z)
  YY = AIMAG(Z)
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
  TOL = MAX(R1MACH(4),1.0E-18)
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
  AZ = ABS(Z)
!-----------------------------------------------------------------------
!     TEST FOR RANGE
!-----------------------------------------------------------------------
  AA = 0.5E0/TOL
  BB=I1MACH(9)*0.5E0
  AA=MIN(AA,BB)
  if ( AZ > AA) go to 140
  FN=FNU+(N-1)
  if ( FN > AA) go to 140
  AA=SQRT(AA)
  if ( AZ > AA) IERR=3
  if ( FN > AA) IERR=3
  ZN = Z
  CSGN = CONE
  if (XX >= 0.0E0) go to 40
  ZN = -Z
!-----------------------------------------------------------------------
!     CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
  INU = FNU
  ARG = (FNU-INU)*PI
  if (YY < 0.0E0) ARG = -ARG
  S1 = COS(ARG)
  S2 = SIN(ARG)
  CSGN = CMPLX(S1,S2)
  if (MOD(INU,2) == 1) CSGN = -CSGN
   40 CONTINUE
  call CBINU(ZN, FNU, KODE, N, CY, NZ, RL, FNUL, TOL, ELIM, ALIM)
  if (NZ < 0) go to 120
  if (XX >= 0.0E0) RETURN
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE
!-----------------------------------------------------------------------
  NN = N - NZ
  if (NN == 0) RETURN
  RTOL = 1.0E0/TOL
  ASCLE = R1MACH(1)*RTOL*1.0E+3
  DO 50 I=1,NN
!       CY(I) = CY(I)*CSGN
    ZN=CY(I)
    AA=REAL(ZN)
    BB=AIMAG(ZN)
    ATOL=1.0E0
    if (MAX(ABS(AA),ABS(BB)) > ASCLE) go to 55
      ZN = ZN*CMPLX(RTOL,0.0E0)
      ATOL = TOL
   55   CONTINUE
    ZN = ZN*CSGN
    CY(I) = ZN*CMPLX(ATOL,0.0E0)
    CSGN = -CSGN
   50 CONTINUE
  return
  120 CONTINUE
  if ( NZ == (-2)) go to 130
  NZ = 0
  IERR=2
  return
  130 CONTINUE
  NZ=0
  IERR=5
  return
  140 CONTINUE
  NZ=0
  IERR=4
  return
end
