subroutine HSTSSP (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, &
     BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
!
!! HSTSSP solves the standard five-point finite difference approximation ...
!  on a staggered grid to the Helmholtz equation in spherical coordinates
!  and on the surface of the unit sphere (radius of 1).
!
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B1A1A
!***TYPE      SINGLE PRECISION (HSTSSP-S)
!***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SPHERICAL
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!     HSTSSP solves the standard five-point finite difference
!     approximation on a staggered grid to the Helmholtz equation in
!     spherical coordinates and on the surface of the unit sphere
!     (radius of 1)
!
!             (1/SIN(THETA))(d/dTHETA)(SIN(THETA)(dU/dTHETA)) +
!
!       (1/SIN(THETA)**2)(d/dPHI)(dU/dPHI) + LAMBDA*U = F(THETA,PHI)
!
!     where THETA is colatitude and PHI is longitude.
!
!    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!    * * * * * * * *    Parameter Description     * * * * * * * * * *
!
!            * * * * * *   On Input    * * * * * *
!
!   A,B
!     The range of THETA (colatitude), i.e. A  <=  THETA  <=  B.  A
!     must be less than B and A must be non-negative.  A and B are in
!     radians.  A = 0 corresponds to the north pole and B = PI
!     corresponds to the south pole.
!
!
!                  * * *  IMPORTANT  * * *
!
!     If B is equal to PI, then B must be computed using the statement
!
!     B = PIMACH(DUM)
!
!     This insures that B in the user's program is equal to PI in this
!     program which permits several tests of the input parameters that
!     otherwise would not be possible.
!
!                  * * * * * * * * * * * *
!
!
!
!   M
!     The number of grid points in the interval (A,B).  The grid points
!     in the THETA-direction are given by THETA(I) = A + (I-0.5)DTHETA
!     for I=1,2,...,M where DTHETA =(B-A)/M.  M must be greater than 2.
!
!   MBDCND
!     Indicates the type of boundary conditions at THETA = A and
!     THETA = B.
!
!     = 1  If the solution is specified at THETA = A and THETA = B.
!          (see note 3 below)
!
!     = 2  If the solution is specified at THETA = A and the derivative
!          of the solution with respect to THETA is specified at
!          THETA = B (see notes 2 and 3 below).
!
!     = 3  If the derivative of the solution with respect to THETA is
!          specified at THETA = A (see notes 1, 2 below) and THETA = B.
!
!     = 4  If the derivative of the solution with respect to THETA is
!          specified at THETA = A (see notes 1 and 2 below) and the
!          solution is specified at THETA = B.
!
!     = 5  If the solution is unspecified at THETA = A = 0 and the
!          solution is specified at THETA = B.  (see note 3 below)
!
!     = 6  If the solution is unspecified at THETA = A = 0 and the
!          derivative of the solution with respect to THETA is
!          specified at THETA = B (see note 2 below).
!
!     = 7  If the solution is specified at THETA = A and the
!          solution is unspecified at THETA = B = PI. (see note 3 below)
!
!     = 8  If the derivative of the solution with respect to
!          THETA is specified at THETA = A (see note 1 below)
!          and the solution is unspecified at THETA = B = PI.
!
!     = 9  If the solution is unspecified at THETA = A = 0 and
!          THETA = B = PI.
!
!     NOTES:  1.  If A = 0, do not use MBDCND = 3, 4, or 8,
!                 but instead use MBDCND = 5, 6, or 9.
!
!             2.  If B = PI, do not use MBDCND = 2, 3, or 6,
!                 but instead use MBDCND = 7, 8, or 9.
!
!             3.  When the solution is specified at THETA = 0 and/or
!                 THETA = PI and the other boundary conditions are
!                 combinations of unspecified, normal derivative, or
!                 periodicity a singular system results.  The unique
!                 solution is determined by extrapolation to the
!                 specification of the solution at either THETA = 0 or
!                 THETA = PI.  But in these cases the right side of the
!                 system will be perturbed by the constant PERTRB.
!
!   BDA
!     A one-dimensional array of length N that specifies the boundary
!     values (if any) of the solution at THETA = A.  When
!     MBDCND = 1, 2, or 7,
!
!              BDA(J) = U(A,PHI(J)) ,              J=1,2,...,N.
!
!     When MBDCND = 3, 4, or 8,
!
!              BDA(J) = (d/dTHETA)U(A,PHI(J)) ,    J=1,2,...,N.
!
!     When MBDCND has any other value, BDA is a dummy variable.
!
!   BDB
!     A one-dimensional array of length N that specifies the boundary
!     values of the solution at THETA = B.  When MBDCND = 1,4, or 5,
!
!              BDB(J) = U(B,PHI(J)) ,              J=1,2,...,N.
!
!     When MBDCND = 2,3, or 6,
!
!              BDB(J) = (d/dTHETA)U(B,PHI(J)) ,    J=1,2,...,N.
!
!     When MBDCND has any other value, BDB is a dummy variable.
!
!   C,D
!     The range of PHI (longitude), i.e. C  <=  PHI  <=  D.
!     C must be less than D.  If D-C = 2*PI, periodic boundary
!     conditions are usually prescribed.
!
!   N
!     The number of unknowns in the interval (C,D).  The unknowns in
!     the PHI-direction are given by PHI(J) = C + (J-0.5)DPHI,
!     J=1,2,...,N, where DPHI = (D-C)/N.  N must be greater than 2.
!
!   NBDCND
!     Indicates the type of boundary conditions at PHI = C
!     and PHI = D.
!
!     = 0  If the solution is periodic in PHI, i.e.
!          U(I,J) = U(I,N+J).
!
!     = 1  If the solution is specified at PHI = C and PHI = D
!          (see note below).
!
!     = 2  If the solution is specified at PHI = C and the derivative
!          of the solution with respect to PHI is specified at
!          PHI = D (see note below).
!
!     = 3  If the derivative of the solution with respect to PHI is
!          specified at PHI = C and PHI = D.
!
!     = 4  If the derivative of the solution with respect to PHI is
!          specified at PHI = C and the solution is specified at
!          PHI = D (see note below).
!
!     NOTE:  When NBDCND = 1, 2, or 4, do not use MBDCND = 5, 6, 7, 8,
!     or 9 (the former indicates that the solution is specified at
!     a pole; the latter indicates the solution is unspecified).  Use
!     instead MBDCND = 1 or 2.
!
!   BDC
!     A one dimensional array of length M that specifies the boundary
!     values of the solution at PHI = C.   When NBDCND = 1 or 2,
!
!              BDC(I) = U(THETA(I),C) ,              I=1,2,...,M.
!
!     When NBDCND = 3 or 4,
!
!              BDC(I) = (d/dPHI)U(THETA(I),C),       I=1,2,...,M.
!
!     When NBDCND = 0, BDC is a dummy variable.
!
!   BDD
!     A one-dimensional array of length M that specifies the boundary
!     values of the solution at PHI = D.  When NBDCND = 1 or 4,
!
!              BDD(I) = U(THETA(I),D) ,              I=1,2,...,M.
!
!     When NBDCND = 2 or 3,
!
!              BDD(I) = (d/dPHI)U(THETA(I),D) ,      I=1,2,...,M.
!
!     When NBDCND = 0, BDD is a dummy variable.
!
!   ELMBDA
!     The constant LAMBDA in the Helmholtz equation.  If LAMBDA is
!     greater than 0, a solution may not exist.  However, HSTSSP will
!     attempt to find a solution.
!
!   F
!     A two-dimensional array that specifies the values of the right
!     side of the Helmholtz equation.  For I=1,2,...,M and J=1,2,...,N
!
!              F(I,J) = F(THETA(I),PHI(J)) .
!
!     F must be dimensioned at least M X N.
!
!   IDIMF
!     The row (or first) dimension of the array F as it appears in the
!     program calling HSTSSP.  This parameter is used to specify the
!     variable dimension of F.  IDIMF must be at least M.
!
!   W
!     A one-dimensional array that must be provided by the user for
!     work space.  W may require up to 13M + 4N + M*INT(log2(N))
!     locations.  The actual number of locations used is computed by
!     HSTSSP and is returned in the location W(1).
!
!
!            * * * * * *   On Output   * * * * * *
!
!   F
!     Contains the solution U(I,J) of the finite difference
!     approximation for the grid point (THETA(I),PHI(J)) for
!     I=1,2,...,M, J=1,2,...,N.
!
!   PERTRB
!     If a combination of periodic, derivative, or unspecified
!     boundary conditions is specified for a Poisson equation
!     (LAMBDA = 0), a solution may not exist.  PERTRB is a con-
!     stant, calculated and subtracted from F, which ensures
!     that a solution exists.  HSTSSP then computes this
!     solution, which is a least squares solution to the
!     original approximation.  This solution plus any constant is also
!     a solution; hence, the solution is not unique.  The value of
!     PERTRB should be small compared to the right side F.
!     Otherwise, a solution is obtained to an essentially different
!     problem.  This comparison should always be made to insure that
!     a meaningful solution has been obtained.
!
!   IERROR
!     An error flag that indicates invalid input parameters.
!      Except for numbers 0 and 14, a solution is not attempted.
!
!     =  0  No error
!
!     =  1  A  <  0 or B  >  PI
!
!     =  2  A  >=  B
!
!     =  3  MBDCND  <  1 or MBDCND  >  9
!
!     =  4  C  >=  D
!
!     =  5  N  <=  2
!
!     =  6  NBDCND  <  0 or NBDCND  >  4
!
!     =  7  A  >  0 and MBDCND = 5, 6, or 9
!
!     =  8  A = 0 and MBDCND = 3, 4, or 8
!
!     =  9  B  <  PI and MBDCND  >=  7
!
!     = 10  B = PI and MBDCND = 2,3, or 6
!
!     = 11  MBDCND  >=  5 and NDBCND = 1, 2, or 4
!
!     = 12  IDIMF  <  M
!
!     = 13  M  <=  2
!
!     = 14  LAMBDA  >  0
!
!     Since this is the only means of indicating a possibly
!     incorrect call to HSTSSP, the user should test IERROR after
!     the call.
!
!   W
!     W(1) contains the required length of W.
!
! *Long Description:
!
!    * * * * * * *   Program Specifications    * * * * * * * * * * * *
!
!    Dimension of   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N),
!    Arguments      W(see argument list)
!
!    Latest         June 1, 1977
!    Revision
!
!    Subprograms    HSTSSP,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2,
!    Required       COSGEN,MERGE,TRIX,TRI3,PIMACH
!
!    Special        NONE
!    Conditions
!
!    Common         NONE
!    Blocks
!
!    I/O            NONE
!
!    Precision      Single
!
!    Specialist     Roland Sweet
!
!    Language       FORTRAN
!
!    History        Written by Roland Sweet at NCAR in April, 1977
!
!    Algorithm      This subroutine defines the finite-difference
!                   equations, incorporates boundary data, adjusts the
!                   right side when the system is singular and calls
!                   either POISTG or GENBUN which solves the linear
!                   system of equations.
!
!    Space          8427(decimal) = 20353(octal) locations on the
!    Required       NCAR Control Data 7600
!
!     Timing and        The execution time T on the NCAR Control Data
!     Accuracy       7600 for subroutine HSTSSP is roughly proportional
!                    to M*N*log2(N).  Some typical values are listed in
!                    the table below.
!                       The solution process employed results in a loss
!                    of no more than four significant digits for N and M
!                    as large as 64.  More detailed information about
!                    accuracy can be found in the documentation for
!                    subroutine POISTG which is the routine that
!                    actually solves the finite difference equations.
!
!
!                       M(=N)    MBDCND    NBDCND    T(MSECS)
!                       -----    ------    ------    --------
!
!                        32       1-9       1-4         56
!                        64       1-9       1-4        230
!
!    Portability     American National Standards Institute FORTRAN.
!                    The machine dependent constant PI is defined in
!                    function PIMACH.
!
!    Required       COS
!    Resident
!    Routines
!
!    Reference      Schumann, U. and R. Sweet,'A Direct Method For
!                   The Solution Of Poisson's Equation With Neumann
!                   Boundary Conditions On A Staggered Grid Of
!                   Arbitrary Size,' J. Comp. Phys. 20(1976),
!                   pp. 171-182.
!
!    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!***REFERENCES  U. Schumann and R. Sweet, A direct method for the
!                 solution of Poisson's equation with Neumann boundary
!                 conditions on a staggered grid of arbitrary size,
!                 Journal of Computational Physics 20, (1976),
!                 pp. 171-182.
!***ROUTINES CALLED  GENBUN, PIMACH, POISTG
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  HSTSSP
!
!
  DIMENSION       F(IDIMF,*) ,BDA(*)     ,BDB(*)     ,BDC(*)     , &
                  BDD(*)     ,W(*)
!***FIRST EXECUTABLE STATEMENT  HSTSSP
  IERROR = 0
  PI = PIMACH(DUM)
  if (A < 0. .OR. B > PI) IERROR = 1
  if (A  >=  B) IERROR = 2
  if (MBDCND <= 0 .OR. MBDCND > 9) IERROR = 3
  if (C  >=  D) IERROR = 4
  if (N  <=  2) IERROR = 5
  if (NBDCND < 0 .OR. NBDCND >= 5) IERROR = 6
  if (A > 0. .AND. (MBDCND == 5 .OR. MBDCND == 6 .OR. MBDCND == 9)) &
      IERROR = 7
  if (A == 0. .AND. (MBDCND == 3 .OR. MBDCND == 4 .OR. MBDCND == 8)) &
      IERROR = 8
  if (B < PI .AND. MBDCND >= 7) IERROR = 9
  if (B == PI .AND. (MBDCND == 2 .OR. MBDCND == 3 .OR. MBDCND == 6)) &
      IERROR = 10
  if (MBDCND >= 5 .AND. &
      (NBDCND == 1 .OR. NBDCND == 2 .OR. NBDCND == 4)) IERROR = 11
  if (IDIMF  <  M) IERROR = 12
  if (M  <=  2) IERROR = 13
  if (IERROR  /=  0) RETURN
  DELTAR = (B-A)/M
  DLRSQ = DELTAR**2
  DELTHT = (D-C)/N
  DLTHSQ = DELTHT**2
  NP = NBDCND+1
  ISW = 1
  JSW = 1
  MB = MBDCND
  if (ELMBDA  /=  0.) go to 105
  go to (101,102,105,103,101,105,101,105,105),MBDCND
  101 if (A /= 0. .OR. B /= PI) go to 105
  MB = 9
  go to 104
  102 if (A  /=  0.) go to 105
  MB = 6
  go to 104
  103 if (B  /=  PI) go to 105
  MB = 8
  104 JSW = 2
  105 CONTINUE
!
!     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
!
  IWB = M
  IWC = IWB+M
  IWR = IWC+M
  IWS = IWR+M
  DO 106 I=1,M
     J = IWR+I
     W(J) = SIN(A+(I-0.5)*DELTAR)
     W(I) = SIN((A+(I-1)*DELTAR))/DLRSQ
  106 CONTINUE
  MM1 = M-1
  DO 107 I=1,MM1
     K = IWC+I
     W(K) = W(I+1)
     J = IWR+I
     K = IWB+I
     W(K) = ELMBDA*W(J)-(W(I)+W(I+1))
  107 CONTINUE
  W(IWR) = SIN(B)/DLRSQ
  W(IWC) = ELMBDA*W(IWS)-(W(M)+W(IWR))
  DO 109 I=1,M
     J = IWR+I
     A1 = W(J)
     DO 108 J=1,N
        F(I,J) = A1*F(I,J)
  108    CONTINUE
  109 CONTINUE
!
!     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
!
  go to (110,110,112,112,114,114,110,112,114),MB
  110 A1 = 2.*W(1)
  W(IWB+1) = W(IWB+1)-W(1)
  DO 111 J=1,N
     F(1,J) = F(1,J)-A1*BDA(J)
  111 CONTINUE
  go to 114
  112 A1 = DELTAR*W(1)
  W(IWB+1) = W(IWB+1)+W(1)
  DO 113 J=1,N
     F(1,J) = F(1,J)+A1*BDA(J)
  113 CONTINUE
  114 go to (115,117,117,115,115,117,119,119,119),MB
  115 A1 = 2.*W(IWR)
  W(IWC) = W(IWC)-W(IWR)
  DO 116 J=1,N
     F(M,J) = F(M,J)-A1*BDB(J)
  116 CONTINUE
  go to 119
  117 A1 = DELTAR*W(IWR)
  W(IWC) = W(IWC)+W(IWR)
  DO 118 J=1,N
     F(M,J) = F(M,J)-A1*BDB(J)
  118 CONTINUE
!
!     ENTER BOUNDARY DATA FOR PHI-BOUNDARIES.
!
  119 A1 = 2./DLTHSQ
  go to (129,120,120,122,122),NP
  120 DO 121 I=1,M
     J = IWR+I
     F(I,1) = F(I,1)-A1*BDC(I)/W(J)
  121 CONTINUE
  go to 124
  122 A1 = 1./DELTHT
  DO 123 I=1,M
     J = IWR+I
     F(I,1) = F(I,1)+A1*BDC(I)/W(J)
  123 CONTINUE
  124 A1 = 2./DLTHSQ
  go to (129,125,127,127,125),NP
  125 DO 126 I=1,M
     J = IWR+I
     F(I,N) = F(I,N)-A1*BDD(I)/W(J)
  126 CONTINUE
  go to 129
  127 A1 = 1./DELTHT
  DO 128 I=1,M
     J = IWR+I
     F(I,N) = F(I,N)-A1*BDD(I)/W(J)
  128 CONTINUE
  129 CONTINUE
!
!     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
!     SOLUTION.
!
  PERTRB = 0.
  if (ELMBDA) 139,131,130
  130 IERROR = 14
  go to 139
  131 go to (139,139,132,139,139,132,139,132,132),MB
  132 go to (133,139,139,133,139),NP
  133 CONTINUE
  ISW = 2
  DO 135 J=1,N
     DO 134 I=1,M
        PERTRB = PERTRB+F(I,J)
  134    CONTINUE
  135 CONTINUE
  A1 = N*(COS(A)-COS(B))/(2.*SIN(0.5*DELTAR))
  PERTRB = PERTRB/A1
  DO 137 I=1,M
     J = IWR+I
     A1 = PERTRB*W(J)
     DO 136 J=1,N
        F(I,J) = F(I,J)-A1
  136    CONTINUE
  137 CONTINUE
  A2 = 0.
  A3 = 0.
  DO 138 J=1,N
     A2 = A2+F(1,J)
     A3 = A3+F(M,J)
  138 CONTINUE
  A2 = A2/W(IWR+1)
  A3 = A3/W(IWS)
  139 CONTINUE
!
!     MULTIPLY I-TH EQUATION THROUGH BY  R(I)*DELTHT**2
!
  DO 141 I=1,M
     J = IWR+I
     A1 = DLTHSQ*W(J)
     W(I) = A1*W(I)
     J = IWC+I
     W(J) = A1*W(J)
     J = IWB+I
     W(J) = A1*W(J)
     DO 140 J=1,N
        F(I,J) = A1*F(I,J)
  140    CONTINUE
  141 CONTINUE
  LP = NBDCND
  W(1) = 0.
  W(IWR) = 0.
!
!     call POISTG OR GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
!
  if (NBDCND  ==  0) go to 142
  call POISTG (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
  go to 143
  142 call GENBUN (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
  143 CONTINUE
  W(1) = W(IWR+1)+3*M
  if (ISW /= 2 .OR. JSW /= 2) go to 150
  if (MB  /=  8) go to 145
  A1 = 0.
  DO 144 J=1,N
     A1 = A1+F(M,J)
  144 CONTINUE
  A1 = (A1-DLRSQ*A3/16.)/N
  if (NBDCND  ==  3) A1 = A1+(BDD(M)-BDC(M))/(D-C)
  A1 = BDB(1)-A1
  go to 147
  145 A1 = 0.
  DO 146 J=1,N
     A1 = A1+F(1,J)
  146 CONTINUE
  A1 = (A1-DLRSQ*A2/16.)/N
  if (NBDCND  ==  3) A1 = A1+(BDD(1)-BDC(1))/(D-C)
  A1 = BDA(1)-A1
  147 DO 149 I=1,M
     DO 148 J=1,N
        F(I,J) = F(I,J)+A1
  148    CONTINUE
  149 CONTINUE
  150 CONTINUE
  return
end
