subroutine HWSCYL (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, &
     BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
!
!! HWSCYL solves a standard finite difference approximation ...
!            to the Helmholtz equation in cylindrical coordinates.
!
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B1A1A
!***TYPE      SINGLE PRECISION (HWSCYL-S)
!***KEYWORDS  CYLINDRICAL, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!     Subroutine HWSCYL solves a finite difference approximation to the
!     Helmholtz equation in cylindrical coordinates:
!
!          (1/R)(d/dR)(R(dU/dR)) + (d/dZ)(dU/dZ)
!
!                                + (LAMBDA/R**2)U = F(R,Z)
!
!     This modified Helmholtz equation results from the Fourier
!     transform of the three-dimensional Poisson equation.
!
!     * * * * * * * *    Parameter Description     * * * * * * * * * *
!
!             * * * * * *   On Input    * * * * * *
!
!     A,B
!       The range of R, i.e., A  <=  R  <=  B.  A must be less than B
!       and A must be non-negative.
!
!     M
!       The number of panels into which the interval (A,B) is
!       subdivided.  Hence, there will be M+1 grid points in the
!       R-direction given by R(I) = A+(I-1)DR, for I = 1,2,...,M+1,
!       where DR = (B-A)/M is the panel width. M must be greater than 3.
!
!     MBDCND
!       Indicates the type of boundary conditions at R = A and R = B.
!
!       = 1  If the solution is specified at R = A and R = B.
!       = 2  If the solution is specified at R = A and the derivative of
!            the solution with respect to R is specified at R = B.
!       = 3  If the derivative of the solution with respect to R is
!            specified at R = A (see note below) and R = B.
!       = 4  If the derivative of the solution with respect to R is
!            specified at R = A (see note below) and the solution is
!            specified at R = B.
!       = 5  If the solution is unspecified at R = A = 0 and the
!            solution is specified at R = B.
!       = 6  If the solution is unspecified at R = A = 0 and the
!            derivative of the solution with respect to R is specified
!            at R = B.
!
!       NOTE:  If A = 0, do not use MBDCND = 3 or 4, but instead use
!              MBDCND = 1,2,5, or 6  .
!
!     BDA
!       A one-dimensional array of length N+1 that specifies the values
!       of the derivative of the solution with respect to R at R = A.
!       When MBDCND = 3 or 4,
!
!            BDA(J) = (d/dR)U(A,Z(J)), J = 1,2,...,N+1  .
!
!       When MBDCND has any other value, BDA is a dummy variable.
!
!     BDB
!       A one-dimensional array of length N+1 that specifies the values
!       of the derivative of the solution with respect to R at R = B.
!       When MBDCND = 2,3, or 6,
!
!            BDB(J) = (d/dR)U(B,Z(J)), J = 1,2,...,N+1  .
!
!       When MBDCND has any other value, BDB is a dummy variable.
!
!     C,D
!       The range of Z, i.e., C  <=  Z  <=  D.  C must be less than D.
!
!     N
!       The number of panels into which the interval (C,D) is
!       subdivided.  Hence, there will be N+1 grid points in the
!       Z-direction given by Z(J) = C+(J-1)DZ, for J = 1,2,...,N+1,
!       where DZ = (D-C)/N is the panel width. N must be greater than 3.
!
!     NBDCND
!       Indicates the type of boundary conditions at Z = C and Z = D.
!
!       = 0  If the solution is periodic in Z, i.e., U(I,1) = U(I,N+1).
!       = 1  If the solution is specified at Z = C and Z = D.
!       = 2  If the solution is specified at Z = C and the derivative of
!            the solution with respect to Z is specified at Z = D.
!       = 3  If the derivative of the solution with respect to Z is
!            specified at Z = C and Z = D.
!       = 4  If the derivative of the solution with respect to Z is
!            specified at Z = C and the solution is specified at Z = D.
!
!     BDC
!       A one-dimensional array of length M+1 that specifies the values
!       of the derivative of the solution with respect to Z at Z = C.
!       When NBDCND = 3 or 4,
!
!            BDC(I) = (d/dZ)U(R(I),C), I = 1,2,...,M+1  .
!
!       When NBDCND has any other value, BDC is a dummy variable.
!
!     BDD
!       A one-dimensional array of length M+1 that specifies the values
!       of the derivative of the solution with respect to Z at Z = D.
!       When NBDCND = 2 or 3,
!
!            BDD(I) = (d/dZ)U(R(I),D), I = 1,2,...,M+1  .
!
!       When NBDCND has any other value, BDD is a dummy variable.
!
!     ELMBDA
!       The constant LAMBDA in the Helmholtz equation.  If
!       LAMBDA  >  0, a solution may not exist.  However, HWSCYL will
!       attempt to find a solution.  LAMBDA must be zero when
!       MBDCND = 5 or 6  .
!
!     F
!       A two-dimensional array that specifies the values of the right
!       side of the Helmholtz equation and boundary data (if any).  For
!       I = 2,3,...,M and J = 2,3,...,N
!
!            F(I,J) = F(R(I),Z(J)).
!
!       On the boundaries F is defined by
!
!            MBDCND   F(1,J)            F(M+1,J)
!            ------   ---------         ---------
!
!              1      U(A,Z(J))         U(B,Z(J))
!              2      U(A,Z(J))         F(B,Z(J))
!              3      F(A,Z(J))         F(B,Z(J))   J = 1,2,...,N+1
!              4      F(A,Z(J))         U(B,Z(J))
!              5      F(0,Z(J))         U(B,Z(J))
!              6      F(0,Z(J))         F(B,Z(J))
!
!            NBDCND   F(I,1)            F(I,N+1)
!            ------   ---------         ---------
!
!              0      F(R(I),C)         F(R(I),C)
!              1      U(R(I),C)         U(R(I),D)
!              2      U(R(I),C)         F(R(I),D)   I = 1,2,...,M+1
!              3      F(R(I),C)         F(R(I),D)
!              4      F(R(I),C)         U(R(I),D)
!
!       F must be dimensioned at least (M+1)*(N+1).
!
!       NOTE
!
!       If the table calls for both the solution U and the right side F
!       at a corner then the solution must be specified.
!
!     IDIMF
!       The row (or first) dimension of the array F as it appears in the
!       program calling HWSCYL.  This parameter is used to specify the
!       variable dimension of F.  IDIMF must be at least M+1  .
!
!     W
!       A one-dimensional array that must be provided by the user for
!       work space.  W may require up to 4*(N+1) +
!       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of
!       locations used is computed by HWSCYL and is returned in location
!       W(1).
!
!
!             * * * * * *   On Output     * * * * * *
!
!     F
!       Contains the solution U(I,J) of the finite difference
!       approximation for the grid point (R(I),Z(J)), I = 1,2,...,M+1,
!       J = 1,2,...,N+1  .
!
!     PERTRB
!       If one specifies a combination of periodic, derivative, and
!       unspecified boundary conditions for a Poisson equation
!       (LAMBDA = 0), a solution may not exist.  PERTRB is a constant,
!       calculated and subtracted from F, which ensures that a solution
!       exists.  HWSCYL then computes this solution, which is a least
!       squares solution to the original approximation.  This solution
!       plus any constant is also a solution.  Hence, the solution is
!       not unique.  The value of PERTRB should be small compared to the
!       right side F.  Otherwise, a solution is obtained to an
!       essentially different problem.  This comparison should always
!       be made to insure that a meaningful solution has been obtained.
!
!     IERROR
!       An error flag which indicates invalid input parameters.  Except
!       for numbers 0 and 11, a solution is not attempted.
!
!       =  0  No error.
!       =  1  A  <  0  .
!       =  2  A  >=  B.
!       =  3  MBDCND  <  1 or MBDCND  >  6  .
!       =  4  C  >=  D.
!       =  5  N  <=  3
!       =  6  NBDCND  <  0 or NBDCND  >  4  .
!       =  7  A = 0, MBDCND = 3 or 4  .
!       =  8  A  >  0, MBDCND  >=  5  .
!       =  9  A = 0, LAMBDA  /=  0, MBDCND  >=  5  .
!       = 10  IDIMF  <  M+1  .
!       = 11  LAMBDA  >  0  .
!       = 12  M  <=  3
!
!       Since this is the only means of indicating a possibly incorrect
!       call to HWSCYL, the user should test IERROR after the call.
!
!     W
!       W(1) contains the required length of W.
!
! *Long Description:
!
!     * * * * * * *   Program Specifications    * * * * * * * * * * * *
!
!     Dimension of   BDA(N+1),BDB(N+1),BDC(M+1),BDD(M+1),F(IDIMF,N+1),
!     Arguments      W(see argument list)
!
!     Latest         June 1, 1976
!     Revision
!
!     Subprograms    HWSCYL,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,
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
!     Space          5818(decimal) = 13272(octal) locations on the NCAR
!     Required       Control Data 7600
!
!     Timing and        The execution time T on the NCAR Control Data
!     Accuracy       7600 for subroutine HWSCYL is roughly proportional
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
!                        32        1         0          31
!                        32        1         1          23
!                        32        3         3          36
!                        64        1         0         128
!                        64        1         1          96
!                        64        3         3         142
!
!     Portability    American National Standards Institute FORTRAN.
!                    The machine dependent constant PI is defined in
!                    function PIMACH.
!
!     Required       COS
!     Resident
!     Routines
!
!     Reference      Swarztrauber, P. and R. Sweet, 'Efficient FORTRAN
!                    Subprograms for the Solution of Elliptic Equations'
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
!***END PROLOGUE  HWSCYL
!
!
  DIMENSION       F(IDIMF,*)
  DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     , &
                  W(*)
!***FIRST EXECUTABLE STATEMENT  HWSCYL
  IERROR = 0
  if (A  <  0.) IERROR = 1
  if (A  >=  B) IERROR = 2
  if (MBDCND <= 0 .OR. MBDCND >= 7) IERROR = 3
  if (C  >=  D) IERROR = 4
  if (N  <=  3) IERROR = 5
  if (NBDCND <= -1 .OR. NBDCND >= 5) IERROR = 6
  if (A == 0. .AND. (MBDCND == 3 .OR. MBDCND == 4)) IERROR = 7
  if (A > 0. .AND. MBDCND >= 5) IERROR = 8
  if (A == 0. .AND. ELMBDA /= 0. .AND. MBDCND >= 5) IERROR = 9
  if (IDIMF  <  M+1) IERROR = 10
  if (M  <=  3) IERROR = 12
  if (IERROR  /=  0) RETURN
  MP1 = M+1
  DELTAR = (B-A)/M
  DLRBY2 = DELTAR/2.
  DLRSQ = DELTAR**2
  NP1 = N+1
  DELTHT = (D-C)/N
  DLTHSQ = DELTHT**2
  NP = NBDCND+1
!
!     DEFINE RANGE OF INDICES I AND J FOR UNKNOWNS U(I,J).
!
  MSTART = 2
  MSTOP = M
  go to (104,103,102,101,101,102),MBDCND
  101 MSTART = 1
  go to 104
  102 MSTART = 1
  103 MSTOP = MP1
  104 MUNK = MSTOP-MSTART+1
  NSTART = 1
  NSTOP = N
  go to (108,105,106,107,108),NP
  105 NSTART = 2
  go to 108
  106 NSTART = 2
  107 NSTOP = NP1
  108 NUNK = NSTOP-NSTART+1
!
!     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
!
  ID2 = MUNK
  ID3 = ID2+MUNK
  ID4 = ID3+MUNK
  ID5 = ID4+MUNK
  ID6 = ID5+MUNK
  ISTART = 1
  A1 = 2./DLRSQ
  IJ = 0
  if (MBDCND == 3 .OR. MBDCND == 4) IJ = 1
  if (MBDCND  <=  4) go to 109
  W(1) = 0.
  W(ID2+1) = -2.*A1
  W(ID3+1) = 2.*A1
  ISTART = 2
  IJ = 1
  109 DO 110 I=ISTART,MUNK
     R = A+(I-IJ)*DELTAR
     J = ID5+I
     W(J) = R
     J = ID6+I
     W(J) = 1./R**2
     W(I) = (R-DLRBY2)/(R*DLRSQ)
     J = ID3+I
     W(J) = (R+DLRBY2)/(R*DLRSQ)
     K = ID6+I
     J = ID2+I
     W(J) = -A1+ELMBDA*W(K)
  110 CONTINUE
  go to (114,111,112,113,114,112),MBDCND
  111 W(ID2) = A1
  go to 114
  112 W(ID2) = A1
  113 W(ID3+1) = A1*ISTART
  114 CONTINUE
!
!     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
!
  go to (115,115,117,117,119,119),MBDCND
  115 A1 = W(1)
  DO 116 J=NSTART,NSTOP
     F(2,J) = F(2,J)-A1*F(1,J)
  116 CONTINUE
  go to 119
  117 A1 = 2.*DELTAR*W(1)
  DO 118 J=NSTART,NSTOP
     F(1,J) = F(1,J)+A1*BDA(J)
  118 CONTINUE
  119 go to (120,122,122,120,120,122),MBDCND
  120 A1 = W(ID4)
  DO 121 J=NSTART,NSTOP
     F(M,J) = F(M,J)-A1*F(MP1,J)
  121 CONTINUE
  go to 124
  122 A1 = 2.*DELTAR*W(ID4)
  DO 123 J=NSTART,NSTOP
     F(MP1,J) = F(MP1,J)-A1*BDB(J)
  123 CONTINUE
!
!     ENTER BOUNDARY DATA FOR Z-BOUNDARIES.
!
  124 A1 = 1./DLTHSQ
  L = ID5-MSTART+1
  go to (134,125,125,127,127),NP
  125 DO 126 I=MSTART,MSTOP
     F(I,2) = F(I,2)-A1*F(I,1)
  126 CONTINUE
  go to 129
  127 A1 = 2./DELTHT
  DO 128 I=MSTART,MSTOP
     F(I,1) = F(I,1)+A1*BDC(I)
  128 CONTINUE
  129 A1 = 1./DLTHSQ
  go to (134,130,132,132,130),NP
  130 DO 131 I=MSTART,MSTOP
     F(I,N) = F(I,N)-A1*F(I,NP1)
  131 CONTINUE
  go to 134
  132 A1 = 2./DELTHT
  DO 133 I=MSTART,MSTOP
     F(I,NP1) = F(I,NP1)-A1*BDD(I)
  133 CONTINUE
  134 CONTINUE
!
!     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
!     SOLUTION.
!
  PERTRB = 0.
  if (ELMBDA) 146,136,135
  135 IERROR = 11
  go to 146
  136 W(ID5+1) = .5*(W(ID5+2)-DLRBY2)
  go to (146,146,138,146,146,137),MBDCND
  137 W(ID5+1) = .5*W(ID5+1)
  138 go to (140,146,146,139,146),NP
  139 A2 = 2.
  go to 141
  140 A2 = 1.
  141 K = ID5+MUNK
  W(K) = .5*(W(K-1)+DLRBY2)
  S = 0.
  DO 143 I=MSTART,MSTOP
     S1 = 0.
     NSP1 = NSTART+1
     NSTM1 = NSTOP-1
     DO 142 J=NSP1,NSTM1
        S1 = S1+F(I,J)
  142    CONTINUE
     K = I+L
     S = S+(A2*S1+F(I,NSTART)+F(I,NSTOP))*W(K)
  143 CONTINUE
  S2 = M*A+(.75+(M-1)*(M+1))*DLRBY2
  if (MBDCND  ==  3) S2 = S2+.25*DLRBY2
  S1 = (2.+A2*(NUNK-2))*S2
  PERTRB = S/S1
  DO 145 I=MSTART,MSTOP
     DO 144 J=NSTART,NSTOP
        F(I,J) = F(I,J)-PERTRB
  144    CONTINUE
  145 CONTINUE
  146 CONTINUE
!
!     MULTIPLY I-TH EQUATION THROUGH BY DELTHT**2 TO PUT EQUATION INTO
!     CORRECT FORM FOR SUBROUTINE GENBUN.
!
  DO 148 I=MSTART,MSTOP
     K = I-MSTART+1
     W(K) = W(K)*DLTHSQ
     J = ID2+K
     W(J) = W(J)*DLTHSQ
     J = ID3+K
     W(J) = W(J)*DLTHSQ
     DO 147 J=NSTART,NSTOP
        F(I,J) = F(I,J)*DLTHSQ
  147    CONTINUE
  148 CONTINUE
  W(1) = 0.
  W(ID4) = 0.
!
!     call GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
!
  call GENBUN (NBDCND,NUNK,1,MUNK,W(1),W(ID2+1),W(ID3+1),IDIMF, &
               F(MSTART,NSTART),IERR1,W(ID4+1))
  W(1) = W(ID4+1)+3*MUNK
  if (NBDCND  /=  0) go to 150
  DO 149 I=MSTART,MSTOP
     F(I,NP1) = F(I,1)
  149 CONTINUE
  150 CONTINUE
  return
end
