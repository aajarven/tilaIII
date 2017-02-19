subroutine HSTPLR (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, &
     BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
!
!! HSTPLR solves the standard five-point finite difference ...
!            approximation on a staggered grid to the Helmholtz equation
!            in polar coordinates.
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B1A1A
!***TYPE      SINGLE PRECISION (HSTPLR-S)
!***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, POLAR
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!      HSTPLR solves the standard five-point finite difference
!      approximation on a staggered grid to the Helmholtz equation in
!      polar coordinates
!
!      (1/R)(d/DR)(R(dU/DR)) + (1/R**2)(d/dTHETA)(dU/dTHETA)
!
!                      + LAMBDA*U = F(R,THETA)
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!     * * * * * * * *    Parameter Description     * * * * * * * * * *
!
!             * * * * * *   On Input    * * * * * *
!
!    A,B
!      The range of R, i.e. A  <=  R  <=  B.  A must be less than B and
!      A must be non-negative.
!
!    M
!      The number of grid points in the interval (A,B).  The grid points
!      in the R-direction are given by R(I) = A + (I-0.5)DR for
!      I=1,2,...,M where DR =(B-A)/M.  M must be greater than 2.
!
!    MBDCND
!      Indicates the type of boundary conditions at R = A and R = B.
!
!      = 1  If the solution is specified at R = A and R = B.
!
!      = 2  If the solution is specified at R = A and the derivative
!           of the solution with respect to R is specified at R = B.
!           (see note 1 below)
!
!      = 3  If the derivative of the solution with respect to R is
!           specified at R = A (see note 2 below) and R = B.
!
!      = 4  If the derivative of the solution with respect to R is
!           specified at R = A (see note 2 below) and the solution is
!           specified at R = B.
!
!      = 5  If the solution is unspecified at R = A = 0 and the solution
!           is specified at R = B.
!
!      = 6  If the solution is unspecified at R = A = 0 and the
!           derivative of the solution with respect to R is specified at
!           R = B.
!
!      NOTE 1:  If A = 0, MBDCND = 2, and NBDCND = 0 or 3, the system of
!               equations to be solved is singular.  The unique solution
!               is determined by extrapolation to the specification of
!               U(0,THETA(1)).  But in this case the right side of the
!               system will be perturbed by the constant PERTRB.
!
!      NOTE 2:  If A = 0, do not use MBDCND = 3 or 4, but instead use
!               MBDCND = 1,2,5, or 6.
!
!    BDA
!      A one-dimensional array of length N that specifies the boundary
!      values (if any) of the solution at R = A.  When MBDCND = 1 or 2,
!
!               BDA(J) = U(A,THETA(J)) ,          J=1,2,...,N.
!
!      When MBDCND = 3 or 4,
!
!               BDA(J) = (d/dR)U(A,THETA(J)) ,    J=1,2,...,N.
!
!      When MBDCND = 5 or 6, BDA is a dummy variable.
!
!    BDB
!      A one-dimensional array of length N that specifies the boundary
!      values of the solution at R = B.  When MBDCND = 1,4, or 5,
!
!               BDB(J) = U(B,THETA(J)) ,          J=1,2,...,N.
!
!      When MBDCND = 2,3, or 6,
!
!               BDB(J) = (d/dR)U(B,THETA(J)) ,    J=1,2,...,N.
!
!    C,D
!      The range of THETA, i.e. C  <=  THETA  <=  D.  C must be less
!      than D.
!
!    N
!      The number of unknowns in the interval (C,D).  The unknowns in
!      the THETA-direction are given by THETA(J) = C + (J-0.5)DT,
!      J=1,2,...,N, where DT = (D-C)/N.  N must be greater than 2.
!
!    NBDCND
!      Indicates the type of boundary conditions at THETA = C
!      and THETA = D.
!
!      = 0  If the solution is periodic in THETA, i.e.
!           U(I,J) = U(I,N+J).
!
!      = 1  If the solution is specified at THETA = C and THETA = D
!           (see note below).
!
!      = 2  If the solution is specified at THETA = C and the derivative
!           of the solution with respect to THETA is specified at
!           THETA = D (see note below).
!
!      = 3  If the derivative of the solution with respect to THETA is
!           specified at THETA = C and THETA = D.
!
!      = 4  If the derivative of the solution with respect to THETA is
!           specified at THETA = C and the solution is specified at
!           THETA = d (see note below).
!
!      NOTE:  When NBDCND = 1, 2, or 4, do not use MBDCND = 5 or 6 (the
!      former indicates that the solution is specified at R =  0; the
!      latter indicates the solution is unspecified at R = 0).  Use
!      instead MBDCND = 1 or 2.
!
!    BDC
!      A one dimensional array of length M that specifies the boundary
!      values of the solution at THETA = C.   When NBDCND = 1 or 2,
!
!               BDC(I) = U(R(I),C) ,              I=1,2,...,M.
!
!      When NBDCND = 3 or 4,
!
!               BDC(I) = (d/dTHETA)U(R(I),C),     I=1,2,...,M.
!
!      When NBDCND = 0, BDC is a dummy variable.
!
!    BDD
!      A one-dimensional array of length M that specifies the boundary
!      values of the solution at THETA = D.  When NBDCND = 1 or 4,
!
!               BDD(I) = U(R(I),D) ,              I=1,2,...,M.
!
!      When NBDCND = 2 or 3,
!
!               BDD(I) = (d/dTHETA)U(R(I),D) ,    I=1,2,...,M.
!
!      When NBDCND = 0, BDD is a dummy variable.
!
!    ELMBDA
!      The constant LAMBDA in the Helmholtz equation.  If LAMBDA is
!      greater than 0, a solution may not exist.  However, HSTPLR will
!      attempt to find a solution.
!
!    F
!      A two-dimensional array that specifies the values of the right
!      side of the Helmholtz equation.  For I=1,2,...,M and J=1,2,...,N
!
!               F(I,J) = F(R(I),THETA(J)) .
!
!      F must be dimensioned at least M X N.
!
!    IDIMF
!      The row (or first) dimension of the array F as it appears in the
!      program calling HSTPLR.  This parameter is used to specify the
!      variable dimension of F.  IDIMF must be at least M.
!
!    W
!      A one-dimensional array that must be provided by the user for
!      work space.  W may require up to 13M + 4N + M*INT(log2(N))
!      locations.  The actual number of locations used is computed by
!      HSTPLR and is returned in the location W(1).
!
!
!             * * * * * *   On Output   * * * * * *
!
!    F
!      Contains the solution U(I,J) of the finite difference
!      approximation for the grid point (R(I),THETA(J)) for
!      I=1,2,...,M, J=1,2,...,N.
!
!    PERTRB
!      If a combination of periodic, derivative, or unspecified
!      boundary conditions is specified for a Poisson equation
!      (LAMBDA = 0), a solution may not exist.  PERTRB is a con-
!      stant, calculated and subtracted from F, which ensures
!      that a solution exists.  HSTPLR then computes this
!      solution, which is a least squares solution to the
!      original approximation.  This solution plus any constant is also
!      a solution; hence, the solution is not unique.  The value of
!      PERTRB should be small compared to the right side F.
!      Otherwise, a solution is obtained to an essentially different
!      problem.  This comparison should always be made to insure that
!      a meaningful solution has been obtained.
!
!    IERROR
!      An error flag that indicates invalid input parameters.
!      Except for numbers 0 and 11, a solution is not attempted.
!
!      =  0  No error
!
!      =  1  A  <  0
!
!      =  2  A  >=  B
!
!      =  3  MBDCND  <  1 or MBDCND  >  6
!
!      =  4  C  >=  D
!
!      =  5  N  <=  2
!
!      =  6  NBDCND  <  0 or NBDCND  >  4
!
!      =  7  A = 0 and MBDCND = 3 or 4
!
!      =  8  A  >  0 and MBDCND  >=  5
!
!      =  9  MBDCND  >=  5 and NBDCND  /=  0 or 3
!
!      = 10  IDIMF  <  M
!
!      = 11  LAMBDA  >  0
!
!      = 12  M  <=  2
!
!      Since this is the only means of indicating a possibly
!      incorrect call to HSTPLR, the user should test IERROR after
!      the call.
!
!    W
!      W(1) contains the required length of W.
!
! *Long Description:
!
!     * * * * * * *   Program Specifications    * * * * * * * * * * * *
!
!     Dimension of   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N),
!     Arguments      W(see ARGUMENT LIST)
!
!     Latest         June 1, 1977
!     Revision
!
!     Subprograms    HSTPLR,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2,
!     Required       COSGEN,MERGE,TRIX,TRI3,PIMACH
!
!     Special        NONE
!     Conditions
!
!     Common         NONE
!     Blocks
!
!     I/O            NONE
!
!     Precision      Single
!
!     Specialist     Roland Sweet
!
!     Language       FORTRAN
!
!     History        Written by Roland Sweet at NCAR in February, 1977
!
!     Algorithm      This subroutine defines the finite-difference
!                    equations, incorporates boundary data, adjusts the
!                    right side when the system is singular and calls
!                    either POISTG or GENBUN which solves the linear
!                    system of equations.
!
!     Space          8265(decimal) = 20111(octal) LOCATIONS ON THE
!     Required       NCAR Control Data 7600
!
!     Timing and        The execution time T on the NCAR Control Data
!     Accuracy       7600 for subroutine HSTPLR is roughly proportional
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
!                        32       1-6       1-4         56
!                        64       1-6       1-4        230
!
!     Portability    American National Standards Institute Fortran.
!                    The machine dependent constant PI is defined in
!                    function PIMACH.
!
!     Required       COS
!     Resident
!     Routines
!
!     Reference      Schumann, U. and R. Sweet,'A Direct Method For
!                    The Solution Of Poisson's Equation With Neumann
!                    Boundary Conditions On A Staggered Grid of
!                    Arbitrary Size,' J. Comp. Phys. 20(1976),
!                    pp. 171-182.
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!***REFERENCES  U. Schumann and R. Sweet, A direct method for the
!                 solution of Poisson's equation with Neumann boundary
!                 conditions on a staggered grid of arbitrary size,
!                 Journal of Computational Physics 20, (1976),
!                 pp. 171-182.
!***ROUTINES CALLED  GENBUN, POISTG
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  HSTPLR
!
!
  DIMENSION       F(IDIMF,*)
  DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     , &
                  W(*)
!***FIRST EXECUTABLE STATEMENT  HSTPLR
  IERROR = 0
  if (A  <  0.) IERROR = 1
  if (A  >=  B) IERROR = 2
  if (MBDCND <= 0 .OR. MBDCND >= 7) IERROR = 3
  if (C  >=  D) IERROR = 4
  if (N  <=  2) IERROR = 5
  if (NBDCND < 0 .OR. NBDCND >= 5) IERROR = 6
  if (A == 0. .AND. (MBDCND == 3 .OR. MBDCND == 4)) IERROR = 7
  if (A > 0. .AND. MBDCND >= 5) IERROR = 8
  if (MBDCND >= 5 .AND. NBDCND /= 0 .AND. NBDCND /= 3) IERROR = 9
  if (IDIMF  <  M) IERROR = 10
  if (M  <=  2) IERROR = 12
  if (IERROR  /=  0) RETURN
  DELTAR = (B-A)/M
  DLRSQ = DELTAR**2
  DELTHT = (D-C)/N
  DLTHSQ = DELTHT**2
  NP = NBDCND+1
  ISW = 1
  MB = MBDCND
  if (A == 0. .AND. MBDCND == 2) MB = 6
!
!     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
!
  IWB = M
  IWC = IWB+M
  IWR = IWC+M
  DO 101 I=1,M
     J = IWR+I
     W(J) = A+(I-0.5)*DELTAR
     W(I) = (A+(I-1)*DELTAR)/DLRSQ
     K = IWC+I
     W(K) = (A+I*DELTAR)/DLRSQ
     K = IWB+I
     W(K) = (ELMBDA-2./DLRSQ)*W(J)
  101 CONTINUE
  DO 103 I=1,M
     J = IWR+I
     A1 = W(J)
     DO 102 J=1,N
        F(I,J) = A1*F(I,J)
  102    CONTINUE
  103 CONTINUE
!
!     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
!
  go to (104,104,106,106,108,108),MB
  104 A1 = 2.*W(1)
  W(IWB+1) = W(IWB+1)-W(1)
  DO 105 J=1,N
     F(1,J) = F(1,J)-A1*BDA(J)
  105 CONTINUE
  go to 108
  106 A1 = DELTAR*W(1)
  W(IWB+1) = W(IWB+1)+W(1)
  DO 107 J=1,N
     F(1,J) = F(1,J)+A1*BDA(J)
  107 CONTINUE
  108 go to (109,111,111,109,109,111),MB
  109 A1 = 2.*W(IWR)
  W(IWC) = W(IWC)-W(IWR)
  DO 110 J=1,N
     F(M,J) = F(M,J)-A1*BDB(J)
  110 CONTINUE
  go to 113
  111 A1 = DELTAR*W(IWR)
  W(IWC) = W(IWC)+W(IWR)
  DO 112 J=1,N
     F(M,J) = F(M,J)-A1*BDB(J)
  112 CONTINUE
!
!     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
!
  113 A1 = 2./DLTHSQ
  go to (123,114,114,116,116),NP
  114 DO 115 I=1,M
     J = IWR+I
     F(I,1) = F(I,1)-A1*BDC(I)/W(J)
  115 CONTINUE
  go to 118
  116 A1 = 1./DELTHT
  DO 117 I=1,M
     J = IWR+I
     F(I,1) = F(I,1)+A1*BDC(I)/W(J)
  117 CONTINUE
  118 A1 = 2./DLTHSQ
  go to (123,119,121,121,119),NP
  119 DO 120 I=1,M
     J = IWR+I
     F(I,N) = F(I,N)-A1*BDD(I)/W(J)
  120 CONTINUE
  go to 123
  121 A1 = 1./DELTHT
  DO 122 I=1,M
     J = IWR+I
     F(I,N) = F(I,N)-A1*BDD(I)/W(J)
  122 CONTINUE
  123 CONTINUE
!
!     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
!     SOLUTION.
!
  PERTRB = 0.
  if (ELMBDA) 133,125,124
  124 IERROR = 11
  go to 133
  125 go to (133,133,126,133,133,126),MB
  126 go to (127,133,133,127,133),NP
  127 CONTINUE
  ISW = 2
  DO 129 J=1,N
     DO 128 I=1,M
        PERTRB = PERTRB+F(I,J)
  128    CONTINUE
  129 CONTINUE
  PERTRB = PERTRB/(M*N*0.5*(A+B))
  DO 131 I=1,M
     J = IWR+I
     A1 = PERTRB*W(J)
     DO 130 J=1,N
        F(I,J) = F(I,J)-A1
  130    CONTINUE
  131 CONTINUE
  A2 = 0.
  DO 132 J=1,N
     A2 = A2+F(1,J)
  132 CONTINUE
  A2 = A2/W(IWR+1)
  133 CONTINUE
!
!     MULTIPLY I-TH EQUATION THROUGH BY  R(I)*DELTHT**2
!
  DO 135 I=1,M
     J = IWR+I
     A1 = DLTHSQ*W(J)
     W(I) = A1*W(I)
     J = IWC+I
     W(J) = A1*W(J)
     J = IWB+I
     W(J) = A1*W(J)
     DO 134 J=1,N
        F(I,J) = A1*F(I,J)
  134    CONTINUE
  135 CONTINUE
  LP = NBDCND
  W(1) = 0.
  W(IWR) = 0.
!
!     call POISTG OR GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
!
  if (LP  ==  0) go to 136
  call POISTG (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
  go to 137
  136 call GENBUN (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
  137 CONTINUE
  W(1) = W(IWR+1)+3*M
  if (A /= 0. .OR. MBDCND /= 2 .OR. ISW /= 2) go to 141
  A1 = 0.
  DO 138 J=1,N
     A1 = A1+F(1,J)
  138 CONTINUE
  A1 = (A1-DLRSQ*A2/16.)/N
  if (NBDCND  ==  3) A1 = A1+(BDD(1)-BDC(1))/(D-C)
  A1 = BDA(1)-A1
  DO 140 I=1,M
     DO 139 J=1,N
        F(I,J) = F(I,J)+A1
  139    CONTINUE
  140 CONTINUE
  141 CONTINUE
  return
end
