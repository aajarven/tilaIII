subroutine HSTCYL (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, &
     BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
!
!! HSTCYL solves the standard five-point finite difference ...
!            approximation on a staggered grid to the modified
!            Helmholtz equation in cylindrical coordinates.
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B1A1A
!***TYPE      SINGLE PRECISION (HSTCYL-S)
!***KEYWORDS  CYLINDRICAL, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!      HSTCYL solves the standard five-point finite difference
!      approximation on a staggered grid to the modified Helmholtz
!      equation in cylindrical coordinates
!
!          (1/R)(d/dR)(R(dU/dR)) + (d/dZ)(dU/dZ)C
!                      + LAMBDA*(1/R**2)*U = F(R,Z)
!
!      This two-dimensional modified Helmholtz equation results
!      from the Fourier transform of a three-dimensional Poisson
!      equation.
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
!      = 1  If the solution is specified at R = A (see note below) and
!           R = B.
!
!      = 2  If the solution is specified at R = A (see note below) and
!           the derivative of the solution with respect to R is
!           specified at R = B.
!
!      = 3  If the derivative of the solution with respect to R is
!           specified at R = A (see note below) and R = B.
!
!      = 4  If the derivative of the solution with respect to R is
!           specified at R = A (see note below) and the solution is
!           specified at R = B.
!
!      = 5  If the solution is unspecified at R = A = 0 and the solution
!           is specified at R = B.
!
!      = 6  If the solution is unspecified at R = A = 0 and the
!           derivative of the solution with respect to R is specified at
!           R = B.
!
!      NOTE:  If A = 0, do not use MBDCND = 1,2,3, or 4, but instead
!             use MBDCND = 5 or 6.  The resulting approximation gives
!             the only meaningful boundary condition, i.e. dU/dR = 0.
!             (see D. Greenspan, 'Introductory Numerical Analysis Of
!             Elliptic Boundary Value Problems,' Harper and Row, 1965,
!             Chapter 5.)
!
!    BDA
!      A one-dimensional array of length N that specifies the boundary
!      values (if any) of the solution at R = A.  When MBDCND = 1 or 2,
!
!               BDA(J) = U(A,Z(J)) ,          J=1,2,...,N.
!
!      When MBDCND = 3 or 4,
!
!               BDA(J) = (d/dR)U(A,Z(J)) ,    J=1,2,...,N.
!
!      When MBDCND = 5 or 6, BDA is a dummy variable.
!
!    BDB
!      A one-dimensional array of length N that specifies the boundary
!      values of the solution at R = B.  When MBDCND = 1,4, or 5,
!
!               BDB(J) = U(B,Z(J)) ,          J=1,2,...,N.
!
!      When MBDCND = 2,3, or 6,
!
!               BDB(J) = (d/dR)U(B,Z(J)) ,    J=1,2,...,N.
!
!    C,D
!      The range of Z, i.e. C  <=  Z  <=  D.  C must be less
!      than D.
!
!    N
!      The number of unknowns in the interval (C,D).  The unknowns in
!      the Z-direction are given by Z(J) = C + (J-0.5)DZ,
!      J=1,2,...,N, where DZ = (D-C)/N.  N must be greater than 2.
!
!    NBDCND
!      Indicates the type of boundary conditions at Z = C
!      and Z = D.
!
!      = 0  If the solution is periodic in Z, i.e.
!           U(I,J) = U(I,N+J).
!
!      = 1  If the solution is specified at Z = C and Z = D.
!
!      = 2  If the solution is specified at Z = C and the derivative
!           of the solution with respect to Z is specified at
!           Z = D.
!
!      = 3  If the derivative of the solution with respect to Z is
!           specified at Z = C and Z = D.
!
!      = 4  If the derivative of the solution with respect to Z is
!           specified at Z = C and the solution is specified at
!           Z = D.
!
!    BDC
!      A one dimensional array of length M that specifies the boundary
!      values of the solution at Z = C.   When NBDCND = 1 or 2,
!
!               BDC(I) = U(R(I),C) ,              I=1,2,...,M.
!
!      When NBDCND = 3 or 4,
!
!               BDC(I) = (d/dZ)U(R(I),C),         I=1,2,...,M.
!
!      When NBDCND = 0, BDC is a dummy variable.
!
!    BDD
!      A one-dimensional array of length M that specifies the boundary
!      values of the solution at Z = D.  when NBDCND = 1 or 4,
!
!               BDD(I) = U(R(I),D) ,              I=1,2,...,M.
!
!      When NBDCND = 2 or 3,
!
!               BDD(I) = (d/dZ)U(R(I),D) ,        I=1,2,...,M.
!
!      When NBDCND = 0, BDD is a dummy variable.
!
!    ELMBDA
!      The constant LAMBDA in the modified Helmholtz equation.  If
!      LAMBDA is greater than 0, a solution may not exist.  However,
!      HSTCYL will attempt to find a solution.  LAMBDA must be zero
!      when MBDCND = 5 or 6.
!
!    F
!      A two-dimensional array that specifies the values of the right
!      side of the modified Helmholtz equation.  For I=1,2,...,M
!      and J=1,2,...,N
!
!               F(I,J) = F(R(I),Z(J)) .
!
!      F must be dimensioned at least M X N.
!
!    IDIMF
!      The row (or first) dimension of the array F as it appears in the
!      program calling HSTCYL.  This parameter is used to specify the
!      variable dimension of F.  IDIMF must be at least M.
!
!    W
!      A one-dimensional array that must be provided by the user for
!      work space.  W may require up to 13M + 4N + M*INT(log2(N))
!      locations.  The actual number of locations used is computed by
!      HSTCYL and is returned in the location W(1).
!
!
!             * * * * * *   On Output   * * * * * *
!
!    F
!      Contains the solution U(I,J) of the finite difference
!      approximation for the grid point (R(I),Z(J)) for
!      I=1,2,...,M, J=1,2,...,N.
!
!    PERTRB
!      If a combination of periodic, derivative, or unspecified
!      boundary conditions is specified for a Poisson equation
!      (LAMBDA = 0), a solution may not exist.  PERTRB is a con-
!      stant, calculated and subtracted from F, which ensures
!      that a solution exists.  HSTCYL then computes this
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
!      =  7  A = 0 and MBDCND = 1,2,3, or 4
!
!      =  8  A  >  0 and MBDCND  >=  5
!
!      =  9  M  <=  2
!
!      = 10  IDIMF  <  M
!
!      = 11  LAMBDA  >  0
!
!      = 12  A=0, MBDCND  >=  5, ELMBDA  /=  0
!
!      Since this is the only means of indicating a possibly
!      incorrect call to HSTCYL, the user should test IERROR after
!      the call.
!
!    W
!      W(1) contains the required length of W.
!
! *Long Description:
!
!     * * * * * * *   Program Specifications    * * * * * * * * * * * *
!
!     Dimension OF   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N),
!     Arguments      W(see argument list)
!
!     Latest         June 1, 1977
!     Revision
!
!     Subprograms    HSTCYL,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2,
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
!     History        Written by Roland Sweet at NCAR in March, 1977
!
!     Algorithm      This subroutine defines the finite-difference
!                    equations, incorporates boundary data, adjusts the
!                    right side when the system is singular and calls
!                    either POISTG or GENBUN which solves the linear
!                    system of equations.
!
!     Space          8228(decimal) = 20044(octal) locations on the
!     Required       NCAR Control Data 7600
!
!     Timing and        The execution time T on the NCAR Control Data
!     Accuracy       7600 for subroutine HSTCYL is roughly proportional
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
!                    The Solution of Poisson's Equation With Neumann
!                    Boundary Conditions On A Staggered Grid Of
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
!***END PROLOGUE  HSTCYL
!
!
  DIMENSION       F(IDIMF,*) ,BDA(*)     ,BDB(*)     ,BDC(*)     , &
                  BDD(*)     ,W(*)
!***FIRST EXECUTABLE STATEMENT  HSTCYL
  IERROR = 0
  if (A  <  0.) IERROR = 1
  if (A  >=  B) IERROR = 2
  if (MBDCND <= 0 .OR. MBDCND >= 7) IERROR = 3
  if (C  >=  D) IERROR = 4
  if (N  <=  2) IERROR = 5
  if (NBDCND < 0 .OR. NBDCND >= 5) IERROR = 6
  if (A == 0. .AND. MBDCND /= 5 .AND. MBDCND /= 6) IERROR = 7
  if (A > 0. .AND. MBDCND >= 5) IERROR = 8
  if (IDIMF  <  M) IERROR = 10
  if (M  <=  2) IERROR = 9
  if (A == 0. .AND. MBDCND >= 5 .AND. ELMBDA /= 0.) IERROR = 12
  if (IERROR  /=  0) RETURN
  DELTAR = (B-A)/M
  DLRSQ = DELTAR**2
  DELTHT = (D-C)/N
  DLTHSQ = DELTHT**2
  NP = NBDCND+1
!
!     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
!
  IWB = M
  IWC = IWB+M
  IWR = IWC+M
  DO 101 I=1,M
     J = IWR+I
     W(J) = A+(I-0.5)*DELTAR
     W(I) = (A+(I-1)*DELTAR)/(DLRSQ*W(J))
     K = IWC+I
     W(K) = (A+I*DELTAR)/(DLRSQ*W(J))
     K = IWB+I
     W(K) = ELMBDA/W(J)**2-2./DLRSQ
  101 CONTINUE
!
!     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
!
  go to (102,102,104,104,106,106),MBDCND
  102 A1 = 2.*W(1)
  W(IWB+1) = W(IWB+1)-W(1)
  DO 103 J=1,N
     F(1,J) = F(1,J)-A1*BDA(J)
  103 CONTINUE
  go to 106
  104 A1 = DELTAR*W(1)
  W(IWB+1) = W(IWB+1)+W(1)
  DO 105 J=1,N
     F(1,J) = F(1,J)+A1*BDA(J)
  105 CONTINUE
  106 CONTINUE
  go to (107,109,109,107,107,109),MBDCND
  107 W(IWC) = W(IWC)-W(IWR)
  A1 = 2.*W(IWR)
  DO 108 J=1,N
     F(M,J) = F(M,J)-A1*BDB(J)
  108 CONTINUE
  go to 111
  109 W(IWC) = W(IWC)+W(IWR)
  A1 = DELTAR*W(IWR)
  DO 110 J=1,N
     F(M,J) = F(M,J)-A1*BDB(J)
  110 CONTINUE
!
!     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
!
  111 A1 = 2./DLTHSQ
  go to (121,112,112,114,114),NP
  112 DO 113 I=1,M
     F(I,1) = F(I,1)-A1*BDC(I)
  113 CONTINUE
  go to 116
  114 A1 = 1./DELTHT
  DO 115 I=1,M
     F(I,1) = F(I,1)+A1*BDC(I)
  115 CONTINUE
  116 A1 = 2./DLTHSQ
  go to (121,117,119,119,117),NP
  117 DO 118 I=1,M
     F(I,N) = F(I,N)-A1*BDD(I)
  118 CONTINUE
  go to 121
  119 A1 = 1./DELTHT
  DO 120 I=1,M
     F(I,N) = F(I,N)-A1*BDD(I)
  120 CONTINUE
  121 CONTINUE
!
!     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
!     SOLUTION.
!
  PERTRB = 0.
  if (ELMBDA) 130,123,122
  122 IERROR = 11
  go to 130
  123 go to (130,130,124,130,130,124),MBDCND
  124 go to (125,130,130,125,130),NP
  125 CONTINUE
  DO 127 I=1,M
     A1 = 0.
     DO 126 J=1,N
        A1 = A1+F(I,J)
  126    CONTINUE
     J = IWR+I
     PERTRB = PERTRB+A1*W(J)
  127 CONTINUE
  PERTRB = PERTRB/(M*N*0.5*(A+B))
  DO 129 I=1,M
     DO 128 J=1,N
        F(I,J) = F(I,J)-PERTRB
  128    CONTINUE
  129 CONTINUE
  130 CONTINUE
!
!     MULTIPLY I-TH EQUATION THROUGH BY  DELTHT**2
!
  DO 132 I=1,M
     W(I) = W(I)*DLTHSQ
     J = IWC+I
     W(J) = W(J)*DLTHSQ
     J = IWB+I
     W(J) = W(J)*DLTHSQ
     DO 131 J=1,N
        F(I,J) = F(I,J)*DLTHSQ
  131    CONTINUE
  132 CONTINUE
  LP = NBDCND
  W(1) = 0.
  W(IWR) = 0.
!
!     call GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
!
  if (NBDCND  ==  0) go to 133
  call POISTG (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
  go to 134
  133 call GENBUN (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
  134 CONTINUE
  W(1) = W(IWR+1)+3*M
  return
end
