subroutine HWSCRT (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, &
     BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
!
!! HWSCRT solves the standard five-point finite difference ...
!            approximation to the Helmholtz equation in Cartesian
!            coordinates.
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B1A1A
!***TYPE      SINGLE PRECISION (HWSCRT-S)
!***KEYWORDS  CARTESIAN, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!     Subroutine HWSCRT solves the standard five-point finite
!     difference approximation to the Helmholtz equation in Cartesian
!     coordinates:
!
!          (d/dX)(dU/dX) + (d/dY)(dU/dY) + LAMBDA*U = F(X,Y).
!
!
!
!     * * * * * * * *    Parameter Description     * * * * * * * * * *
!
!             * * * * * *   On Input    * * * * * *
!
!     A,B
!       The range of X, i.e., A  <=  X  <=  B.  A must be less than B.
!
!     M
!       The number of panels into which the interval (A,B) is
!       subdivided.  Hence, there will be M+1 grid points in the
!       X-direction given by X(I) = A+(I-1)DX for I = 1,2,...,M+1,
!       where DX = (B-A)/M is the panel width. M must be greater than 3.
!
!     MBDCND
!       Indicates the type of boundary conditions at X = A and X = B.
!
!       = 0  If the solution is periodic in X, i.e., U(I,J) = U(M+I,J).
!       = 1  If the solution is specified at X = A and X = B.
!       = 2  If the solution is specified at X = A and the derivative of
!            the solution with respect to X is specified at X = B.
!       = 3  If the derivative of the solution with respect to X is
!            specified at X = A and X = B.
!       = 4  If the derivative of the solution with respect to X is
!            specified at X = A and the solution is specified at X = B.
!
!     BDA
!       A one-dimensional array of length N+1 that specifies the values
!       of the derivative of the solution with respect to X at X = A.
!       When MBDCND = 3 or 4,
!
!            BDA(J) = (d/dX)U(A,Y(J)), J = 1,2,...,N+1  .
!
!       When MBDCND has any other value, BDA is a dummy variable.
!
!     BDB
!       A one-dimensional array of length N+1 that specifies the values
!       of the derivative of the solution with respect to X at X = B.
!       When MBDCND = 2 or 3,
!
!            BDB(J) = (d/dX)U(B,Y(J)), J = 1,2,...,N+1  .
!
!       When MBDCND has any other value BDB is a dummy variable.
!
!     C,D
!       The range of Y, i.e., C  <=  Y  <=  D.  C must be less than D.
!
!     N
!       The number of panels into which the interval (C,D) is
!       subdivided.  Hence, there will be N+1 grid points in the
!       Y-direction given by Y(J) = C+(J-1)DY for J = 1,2,...,N+1, where
!       DY = (D-C)/N is the panel width.  N must be greater than 3.
!
!     NBDCND
!       Indicates the type of boundary conditions at Y = C and Y = D.
!
!       = 0  If the solution is periodic in Y, i.e., U(I,J) = U(I,N+J).
!       = 1  If the solution is specified at Y = C and Y = D.
!       = 2  If the solution is specified at Y = C and the derivative of
!            the solution with respect to Y is specified at Y = D.
!       = 3  If the derivative of the solution with respect to Y is
!            specified at Y = C and Y = D.
!       = 4  If the derivative of the solution with respect to Y is
!            specified at Y = C and the solution is specified at Y = D.
!
!     BDC
!       A one-dimensional array of length M+1 that specifies the values
!       of the derivative of the solution with respect to Y at Y = C.
!       When NBDCND = 3 or 4,
!
!            BDC(I) = (d/dY)U(X(I),C), I = 1,2,...,M+1  .
!
!       When NBDCND has any other value, BDC is a dummy variable.
!
!     BDD
!       A one-dimensional array of length M+1 that specifies the values
!       of the derivative of the solution with respect to Y at Y = D.
!       When NBDCND = 2 or 3,
!
!            BDD(I) = (d/dY)U(X(I),D), I = 1,2,...,M+1  .
!
!       When NBDCND has any other value, BDD is a dummy variable.
!
!     ELMBDA
!       The constant LAMBDA in the Helmholtz equation.  If
!       LAMBDA  >  0, a solution may not exist.  However, HWSCRT will
!       attempt to find a solution.
!
!     F
!       A two-dimensional array which specifies the values of the right
!       side of the Helmholtz equation and boundary values (if any).
!       For I = 2,3,...,M and J = 2,3,...,N
!
!            F(I,J) = F(X(I),Y(J)).
!
!       On the boundaries F is defined by
!
!            MBDCND     F(1,J)        F(M+1,J)
!            ------     ---------     --------
!
!              0        F(A,Y(J))     F(A,Y(J))
!              1        U(A,Y(J))     U(B,Y(J))
!              2        U(A,Y(J))     F(B,Y(J))     J = 1,2,...,N+1
!              3        F(A,Y(J))     F(B,Y(J))
!              4        F(A,Y(J))     U(B,Y(J))
!
!
!            NBDCND     F(I,1)        F(I,N+1)
!            ------     ---------     --------
!
!              0        F(X(I),C)     F(X(I),C)
!              1        U(X(I),C)     U(X(I),D)
!              2        U(X(I),C)     F(X(I),D)     I = 1,2,...,M+1
!              3        F(X(I),C)     F(X(I),D)
!              4        F(X(I),C)     U(X(I),D)
!
!       F must be dimensioned at least (M+1)*(N+1).
!
!       NOTE:
!
!       If the table calls for both the solution U and the right side F
!       at a corner then the solution must be specified.
!
!     IDIMF
!       The row (or first) dimension of the array F as it appears in the
!       program calling HWSCRT.  This parameter is used to specify the
!       variable dimension of F.  IDIMF must be at least M+1  .
!
!     W
!       A one-dimensional array that must be provided by the user for
!       work space.  W may require up to 4*(N+1) +
!       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of
!       locations used is computed by HWSCRT and is returned in location
!       W(1).
!
!
!             * * * * * *   On Output     * * * * * *
!
!     F
!       Contains the solution U(I,J) of the finite difference
!       approximation for the grid point (X(I),Y(J)), I = 1,2,...,M+1,
!       J = 1,2,...,N+1  .
!
!     PERTRB
!       If a combination of periodic or derivative boundary conditions
!       is specified for a Poisson equation (LAMBDA = 0), a solution may
!       not exist.  PERTRB is a constant, calculated and subtracted from
!       F, which ensures that a solution exists.  HWSCRT then computes
!       this solution, which is a least squares solution to the original
!       approximation.  This solution plus any constant is also a
!       solution.  Hence, the solution is not unique.  The value of
!       PERTRB should be small compared to the right side F.  Otherwise,
!       a solution is obtained to an essentially different problem.
!       This comparison should always be made to insure that a
!       meaningful solution has been obtained.
!
!     IERROR
!       An error flag that indicates invalid input parameters.  Except
!       for numbers 0 and 6, a solution is not attempted.
!
!       = 0  No error.
!       = 1  A  >=  B.
!       = 2  MBDCND  <  0 or MBDCND  >  4  .
!       = 3  C  >=  D.
!       = 4  N  <=  3
!       = 5  NBDCND  <  0 or NBDCND  >  4  .
!       = 6  LAMBDA  >  0  .
!       = 7  IDIMF  <  M+1  .
!       = 8  M  <=  3
!
!       Since this is the only means of indicating a possibly incorrect
!       call to HWSCRT, the user should test IERROR after the call.
!
!     W
!       W(1) contains the required length of W.
!
! *Long Description:
!
!     * * * * * * *   Program Specifications    * * * * * * * * * * * *
!
!
!     Dimension of   BDA(N+1),BDB(N+1),BDC(M+1),BDD(M+1),F(IDIMF,N+1),
!     Arguments      W(see argument list)
!
!     Latest         June 1, 1976
!     Revision
!
!     Subprograms    HWSCRT,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,
!     Required       TRIX,TRI3,PIMACH
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
!     History        Standardized September 1, 1973
!                    Revised April 1, 1976
!
!     Algorithm      The routine defines the finite difference
!                    equations, incorporates boundary data, and adjusts
!                    the right side of singular systems and then calls
!                    GENBUN to solve the system.
!
!     Space          13110(octal) = 5704(decimal) locations on the NCAR
!     Required       Control Data 7600
!
!     Timing and        The execution time T on the NCAR Control Data
!     Accuracy       7600 for subroutine HWSCRT is roughly proportional
!                    to M*N*log2(N), but also depends on the input
!                    parameters NBDCND and MBDCND.  Some typical values
!                    are listed in the table below.
!                       The solution process employed results in a loss
!                    of no more than three significant digits for N and
!                    M as large as 64.  More detailed information about
!                    accuracy can be found in the documentation for
!                    subroutine GENBUN which is the routine that
!                    solves the finite difference equations.
!
!
!                       M(=N)    MBDCND    NBDCND    T(MSECS)
!                       -----    ------    ------    --------
!
!                        32        0         0          31
!                        32        1         1          23
!                        32        3         3          36
!                        64        0         0         128
!                        64        1         1          96
!                        64        3         3         142
!
!     Portability    American National Standards Institute FORTRAN.
!                    The machine dependent constant PI is defined in
!                    function PIMACH.
!
!     Reference      Swarztrauber, P. and R. Sweet, 'Efficient FORTRAN
!                    Subprograms for The Solution Of Elliptic Equations'
!                    NCAR TN/IA-109, July, 1975, 138 pp.
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran
!                 subprograms for the solution of elliptic equations,
!                 NCAR TN/IA-109, July 1975, 138 pp.
!***ROUTINES CALLED  GENBUN
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  HWSCRT
!
!
  DIMENSION       F(IDIMF,*)
  DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     , &
                  W(*)
!***FIRST EXECUTABLE STATEMENT  HWSCRT
  IERROR = 0
  if (A  >=  B) IERROR = 1
  if (MBDCND < 0 .OR. MBDCND > 4) IERROR = 2
  if (C  >=  D) IERROR = 3
  if (N  <=  3) IERROR = 4
  if (NBDCND < 0 .OR. NBDCND > 4) IERROR = 5
  if (IDIMF  <  M+1) IERROR = 7
  if (M  <=  3) IERROR = 8
  if (IERROR  /=  0) RETURN
  NPEROD = NBDCND
  MPEROD = 0
  if (MBDCND  >  0) MPEROD = 1
  DELTAX = (B-A)/M
  TWDELX = 2./DELTAX
  DELXSQ = 1./DELTAX**2
  DELTAY = (D-C)/N
  TWDELY = 2./DELTAY
  DELYSQ = 1./DELTAY**2
  NP = NBDCND+1
  NP1 = N+1
  MP = MBDCND+1
  MP1 = M+1
  NSTART = 1
  NSTOP = N
  NSKIP = 1
  go to (104,101,102,103,104),NP
  101 NSTART = 2
  go to 104
  102 NSTART = 2
  103 NSTOP = NP1
  NSKIP = 2
  104 NUNK = NSTOP-NSTART+1
!
!     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
!
  MSTART = 1
  MSTOP = M
  MSKIP = 1
  go to (117,105,106,109,110),MP
  105 MSTART = 2
  go to 107
  106 MSTART = 2
  MSTOP = MP1
  MSKIP = 2
  107 DO 108 J=NSTART,NSTOP
     F(2,J) = F(2,J)-F(1,J)*DELXSQ
  108 CONTINUE
  go to 112
  109 MSTOP = MP1
  MSKIP = 2
  110 DO 111 J=NSTART,NSTOP
     F(1,J) = F(1,J)+BDA(J)*TWDELX
  111 CONTINUE
  112 go to (113,115),MSKIP
  113 DO 114 J=NSTART,NSTOP
     F(M,J) = F(M,J)-F(MP1,J)*DELXSQ
  114 CONTINUE
  go to 117
  115 DO 116 J=NSTART,NSTOP
     F(MP1,J) = F(MP1,J)-BDB(J)*TWDELX
  116 CONTINUE
  117 MUNK = MSTOP-MSTART+1
!
!     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
!
  go to (127,118,118,120,120),NP
  118 DO 119 I=MSTART,MSTOP
     F(I,2) = F(I,2)-F(I,1)*DELYSQ
  119 CONTINUE
  go to 122
  120 DO 121 I=MSTART,MSTOP
     F(I,1) = F(I,1)+BDC(I)*TWDELY
  121 CONTINUE
  122 go to (123,125),NSKIP
  123 DO 124 I=MSTART,MSTOP
     F(I,N) = F(I,N)-F(I,NP1)*DELYSQ
  124 CONTINUE
  go to 127
  125 DO 126 I=MSTART,MSTOP
     F(I,NP1) = F(I,NP1)-BDD(I)*TWDELY
  126 CONTINUE
!
!    MULTIPLY RIGHT SIDE BY DELTAY**2.
!
  127 DELYSQ = DELTAY*DELTAY
  DO 129 I=MSTART,MSTOP
     DO 128 J=NSTART,NSTOP
        F(I,J) = F(I,J)*DELYSQ
  128    CONTINUE
  129 CONTINUE
!
!     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
!
  ID2 = MUNK
  ID3 = ID2+MUNK
  ID4 = ID3+MUNK
  S = DELYSQ*DELXSQ
  ST2 = 2.*S
  DO 130 I=1,MUNK
     W(I) = S
     J = ID2+I
     W(J) = -ST2+ELMBDA*DELYSQ
     J = ID3+I
     W(J) = S
  130 CONTINUE
  if (MP  ==  1) go to 131
  W(1) = 0.
  W(ID4) = 0.
  131 CONTINUE
  go to (135,135,132,133,134),MP
  132 W(ID2) = ST2
  go to 135
  133 W(ID2) = ST2
  134 W(ID3+1) = ST2
  135 CONTINUE
  PERTRB = 0.
  if (ELMBDA) 144,137,136
  136 IERROR = 6
  go to 144
  137 if ((NBDCND == 0 .OR. NBDCND == 3) .AND. &
      (MBDCND == 0 .OR. MBDCND == 3)) go to 138
  go to 144
!
!     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
!     WILL EXIST.
!
  138 A1 = 1.
  A2 = 1.
  if (NBDCND  ==  3) A2 = 2.
  if (MBDCND  ==  3) A1 = 2.
  S1 = 0.
  MSP1 = MSTART+1
  MSTM1 = MSTOP-1
  NSP1 = NSTART+1
  NSTM1 = NSTOP-1
  DO 140 J=NSP1,NSTM1
     S = 0.
     DO 139 I=MSP1,MSTM1
        S = S+F(I,J)
  139    CONTINUE
     S1 = S1+S*A1+F(MSTART,J)+F(MSTOP,J)
  140 CONTINUE
  S1 = A2*S1
  S = 0.
  DO 141 I=MSP1,MSTM1
     S = S+F(I,NSTART)+F(I,NSTOP)
  141 CONTINUE
  S1 = S1+S*A1+F(MSTART,NSTART)+F(MSTART,NSTOP)+F(MSTOP,NSTART)+ &
       F(MSTOP,NSTOP)
  S = (2.+(NUNK-2)*A2)*(2.+(MUNK-2)*A1)
  PERTRB = S1/S
  DO 143 J=NSTART,NSTOP
     DO 142 I=MSTART,MSTOP
        F(I,J) = F(I,J)-PERTRB
  142    CONTINUE
  143 CONTINUE
  PERTRB = PERTRB/DELYSQ
!
!     SOLVE THE EQUATION.
!
  144 call GENBUN (NPEROD,NUNK,MPEROD,MUNK,W(1),W(ID2+1),W(ID3+1), &
               IDIMF,F(MSTART,NSTART),IERR1,W(ID4+1))
  W(1) = W(ID4+1)+3*MUNK
!
!     FILL IN IDENTICAL VALUES WHEN HAVE PERIODIC BOUNDARY CONDITIONS.
!
  if (NBDCND  /=  0) go to 146
  DO 145 I=MSTART,MSTOP
     F(I,NP1) = F(I,1)
  145 CONTINUE
  146 if (MBDCND  /=  0) go to 148
  DO 147 J=NSTART,NSTOP
     F(MP1,J) = F(1,J)
  147 CONTINUE
  if (NBDCND  ==  0) F(MP1,NP1) = F(1,NP1)
  148 CONTINUE
  return
end
