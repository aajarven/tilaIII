subroutine HSTCRT (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, &
     BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
!
!! HSTCRT solves the standard five-point finite difference ...
!            approximation on a staggered grid to the Helmholtz equation
!            in Cartesian coordinates.
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B1A1A
!***TYPE      SINGLE PRECISION (HSTCRT-S)
!***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!      HSTCRT solves the standard five-point finite difference
!      approximation on a staggered grid to the Helmholtz equation in
!      Cartesian coordinates
!
!      (d/dX)(dU/dX) + (d/dY)(dU/dY) + LAMBDA*U = F(X,Y)
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!     * * * * * * * *    Parameter Description     * * * * * * * * * *
!
!             * * * * * *   On Input    * * * * * *
!
!    A,B
!      The range of X, i.e. A  <=  X  <=  B.  A must be less than B.
!
!    M
!      The number of grid points in the interval (A,B).  The grid points
!      in the X-direction are given by X(I) = A + (I-0.5)dX for
!      I=1,2,...,M where dX =(B-A)/M.  M must be greater than 2.
!
!    MBDCND
!      Indicates the type of boundary conditions at X = A and X = B.
!
!      = 0  If the solution is periodic in X,
!           U(M+I,J) = U(I,J).
!
!      = 1  If the solution is specified at X = A and X = B.
!
!      = 2  If the solution is specified at X = A and the derivative
!           of the solution with respect to X is specified at X = B.
!
!      = 3  If the derivative of the solution with respect to X is
!           specified at X = A  and X = B.
!
!      = 4  If the derivative of the solution with respect to X is
!           specified at X = A  and the solution is specified at X = B.
!
!    BDA
!      A one-dimensional array of length N that specifies the boundary
!      values (if any) of the solution at X = A.  When MBDCND = 1 or 2,
!
!               BDA(J) = U(A,Y(J)) ,          J=1,2,...,N.
!
!      When MBDCND = 3 or 4,
!
!               BDA(J) = (d/dX)U(A,Y(J)) ,    J=1,2,...,N.
!
!    BDB
!      A one-dimensional array of length N that specifies the boundary
!      values of the solution at X = B.  When MBDCND = 1 or 4
!
!               BDB(J) = U(B,Y(J)) ,          J=1,2,...,N.
!
!      When MBDCND = 2 or 3
!
!               BDB(J) = (d/dX)U(B,Y(J)) ,    J=1,2,...,N.
!
!    C,D
!      The range of Y, i.e. C  <=  Y  <=  D.  C must be less
!      than D.
!
!    N
!      The number of unknowns in the interval (C,D).  The unknowns in
!      the Y-direction are given by Y(J) = C + (J-0.5)DY,
!      J=1,2,...,N, where DY = (D-C)/N.  N must be greater than 2.
!
!    NBDCND
!      Indicates the type of boundary conditions at Y = C
!      and Y = D.
!
!      = 0  If the solution is periodic in Y, i.e.
!           U(I,J) = U(I,N+J).
!
!      = 1  If the solution is specified at Y = C and Y = D.
!
!      = 2  If the solution is specified at Y = C and the derivative
!           of the solution with respect to Y is specified at Y = D.
!
!      = 3  If the derivative of the solution with respect to Y is
!           specified at Y = C and Y = D.
!
!      = 4  If the derivative of the solution with respect to Y is
!           specified at Y = C and the solution is specified at Y = D.
!
!    BDC
!      A one dimensional array of length M that specifies the boundary
!      values of the solution at Y = C.   When NBDCND = 1 or 2,
!
!               BDC(I) = U(X(I),C) ,              I=1,2,...,M.
!
!      When NBDCND = 3 or 4,
!
!               BDC(I) = (d/dY)U(X(I),C),     I=1,2,...,M.
!
!      When NBDCND = 0, BDC is a dummy variable.
!
!    BDD
!      A one-dimensional array of length M that specifies the boundary
!      values of the solution at Y = D.  When NBDCND = 1 or 4,
!
!               BDD(I) = U(X(I),D) ,              I=1,2,...,M.
!
!      When NBDCND = 2 or 3,
!
!               BDD(I) = (d/dY)U(X(I),D) ,    I=1,2,...,M.
!
!      When NBDCND = 0, BDD is a dummy variable.
!
!    ELMBDA
!      The constant LAMBDA in the Helmholtz equation.  If LAMBDA is
!      greater than 0, a solution may not exist.  However, HSTCRT will
!      attempt to find a solution.
!
!    F
!      A two-dimensional array that specifies the values of the right
!      side of the Helmholtz equation.  For I=1,2,...,M and J=1,2,...,N
!
!               F(I,J) = F(X(I),Y(J)) .
!
!      F must be dimensioned at least M X N.
!
!    IDIMF
!      The row (or first) dimension of the array F as it appears in the
!      program calling HSTCRT.  This parameter is used to specify the
!      variable dimension of F.  IDIMF must be at least M.
!
!    W
!      A one-dimensional array that must be provided by the user for
!      work space.  W may require up to 13M + 4N + M*INT(log2(N))
!      locations.  The actual number of locations used is computed by
!      HSTCRT and is returned in the location W(1).
!
!
!             * * * * * *   On Output   * * * * * *
!
!    F
!      Contains the solution U(I,J) of the finite difference
!      approximation for the grid point (X(I),Y(J)) for
!      I=1,2,...,M, J=1,2,...,N.
!
!    PERTRB
!      If a combination of periodic or derivative boundary conditions is
!      specified for a Poisson equation (LAMBDA = 0), a solution may not
!      exist.  PERTRB is a constant, calculated and subtracted from F,
!      which ensures that a solution exists.  HSTCRT then computes this
!      solution, which is a least squares solution to the original
!      approximation.  This solution plus any constant is also a
!      solution; hence, the solution is not unique.  The value of PERTRB
!      should be small compared to the right side F.  Otherwise, a
!      solution is obtained to an essentially different problem.  This
!      comparison should always be made to insure that a meaningful
!      solution has been obtained.
!
!    IERROR
!      An error flag that indicates invalid input parameters.
!       Except for numbers 0 and  6, a solution is not attempted.
!
!      =  0  No error
!
!      =  1  A  >=  B
!
!      =  2  MBDCND  <  0 or MBDCND  >  4
!
!      =  3  C  >=  D
!
!      =  4  N  <=  2
!
!      =  5  NBDCND  <  0 or NBDCND  >  4
!
!      =  6  LAMBDA  >  0
!
!      =  7  IDIMF  <  M
!
!      =  8  M  <=  2
!
!      Since this is the only means of indicating a possibly
!      incorrect call to HSTCRT, the user should test IERROR after
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
!     Arguments      W(See argument list)
!
!     Latest         June 1, 1977
!     Revision
!
!     Subprograms    HSTCRT,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2,
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
!     History        Written by Roland Sweet at NCAR in January , 1977
!
!     Algorithm      This subroutine defines the finite-difference
!                    equations, incorporates boundary data, adjusts the
!                    right side when the system is singular and calls
!                    either POISTG or GENBUN which solves the linear
!                    system of equations.
!
!     Space          8131(decimal) = 17703(octal) locations on the
!     Required       NCAR Control Data 7600
!
!     Timing and        The execution time T on the NCAR Control Data
!     Accuracy       7600 for subroutine HSTCRT is roughly proportional
!                    to M*N*log2(N).  Some typical values are listed in
!                    the table below.
!                       The solution process employed results in a loss
!                    of no more than FOUR significant digits for N and M
!                    as large as 64.  More detailed information about
!                    accuracy can be found in the documentation for
!                    subroutine POISTG which is the routine that
!                    actually solves the finite difference equations.
!
!
!                       M(=N)    MBDCND    NBDCND    T(MSECS)
!                       -----    ------    ------    --------
!
!                        32       1-4       1-4         56
!                        64       1-4       1-4        230
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
!                    Boundary Conditions On A Staggered Grid Of
!                    Arbitrary Size,' J. COMP. PHYS. 20(1976),
!                    PP. 171-182.
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
!***END PROLOGUE  HSTCRT
!
!
  DIMENSION       F(IDIMF,*) ,BDA(*)     ,BDB(*)     ,BDC(*)     , &
                  BDD(*)     ,W(*)
!***FIRST EXECUTABLE STATEMENT  HSTCRT
  IERROR = 0
  if (A  >=  B) IERROR = 1
  if (MBDCND < 0 .OR. MBDCND > 4) IERROR = 2
  if (C  >=  D) IERROR = 3
  if (N  <=  2) IERROR = 4
  if (NBDCND < 0 .OR. NBDCND > 4) IERROR = 5
  if (IDIMF  <  M) IERROR = 7
  if (M  <=  2) IERROR = 8
  if (IERROR  /=  0) RETURN
  NPEROD = NBDCND
  MPEROD = 0
  if (MBDCND  >  0) MPEROD = 1
  DELTAX = (B-A)/M
  TWDELX = 1./DELTAX
  DELXSQ = 2./DELTAX**2
  DELTAY = (D-C)/N
  TWDELY = 1./DELTAY
  DELYSQ = DELTAY**2
  TWDYSQ = 2./DELYSQ
  NP = NBDCND+1
  MP = MBDCND+1
!
!     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
!
  ID2 = M
  ID3 = ID2+M
  ID4 = ID3+M
  S = (DELTAY/DELTAX)**2
  ST2 = 2.*S
  DO 101 I=1,M
     W(I) = S
     J = ID2+I
     W(J) = -ST2+ELMBDA*DELYSQ
     J = ID3+I
     W(J) = S
  101 CONTINUE
!
!     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
!
  go to (111,102,102,104,104),MP
  102 DO 103 J=1,N
     F(1,J) = F(1,J)-BDA(J)*DELXSQ
  103 CONTINUE
  W(ID2+1) = W(ID2+1)-W(1)
  go to 106
  104 DO 105 J=1,N
     F(1,J) = F(1,J)+BDA(J)*TWDELX
  105 CONTINUE
  W(ID2+1) = W(ID2+1)+W(1)
  106 go to (111,107,109,109,107),MP
  107 DO 108 J=1,N
     F(M,J) = F(M,J)-BDB(J)*DELXSQ
  108 CONTINUE
  W(ID3) = W(ID3)-W(1)
  go to 111
  109 DO 110 J=1,N
     F(M,J) = F(M,J)-BDB(J)*TWDELX
  110 CONTINUE
  W(ID3) = W(ID3)+W(1)
  111 CONTINUE
!
!     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
!
  go to (121,112,112,114,114),NP
  112 DO 113 I=1,M
     F(I,1) = F(I,1)-BDC(I)*TWDYSQ
  113 CONTINUE
  go to 116
  114 DO 115 I=1,M
     F(I,1) = F(I,1)+BDC(I)*TWDELY
  115 CONTINUE
  116 go to (121,117,119,119,117),NP
  117 DO 118 I=1,M
     F(I,N) = F(I,N)-BDD(I)*TWDYSQ
  118 CONTINUE
  go to 121
  119 DO 120 I=1,M
     F(I,N) = F(I,N)-BDD(I)*TWDELY
  120 CONTINUE
  121 CONTINUE
  DO 123 I=1,M
     DO 122 J=1,N
        F(I,J) = F(I,J)*DELYSQ
  122    CONTINUE
  123 CONTINUE
  if (MPEROD  ==  0) go to 124
  W(1) = 0.
  W(ID4) = 0.
  124 CONTINUE
  PERTRB = 0.
  if (ELMBDA) 133,126,125
  125 IERROR = 6
  go to 133
  126 go to (127,133,133,127,133),MP
  127 go to (128,133,133,128,133),NP
!
!     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
!     WILL EXIST.
!
  128 CONTINUE
  S = 0.
  DO 130 J=1,N
     DO 129 I=1,M
        S = S+F(I,J)
  129    CONTINUE
  130 CONTINUE
  PERTRB = S/(M*N)
  DO 132 J=1,N
     DO 131 I=1,M
        F(I,J) = F(I,J)-PERTRB
  131    CONTINUE
  132 CONTINUE
  PERTRB = PERTRB/DELYSQ
!
!     SOLVE THE EQUATION.
!
  133 CONTINUE
  if (NPEROD  ==  0) go to 134
  call POISTG (NPEROD,N,MPEROD,M,W(1),W(ID2+1),W(ID3+1),IDIMF,F, &
               IERR1,W(ID4+1))
  go to 135
  134 CONTINUE
  call GENBUN (NPEROD,N,MPEROD,M,W(1),W(ID2+1),W(ID3+1),IDIMF,F, &
               IERR1,W(ID4+1))
  135 CONTINUE
  W(1) = W(ID4+1)+3*M
  return
end
