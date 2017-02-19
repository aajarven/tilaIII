subroutine HWSPLR (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, &
     BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
!
!! HWSPLR solves a finite difference approximation to the Helmholtz ...
!            equation in polar coordinates.
!
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B1A1A
!***TYPE      SINGLE PRECISION (HWSPLR-S)
!***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, POLAR
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!     Subroutine HWSPLR solves a finite difference approximation to the
!     Helmholtz equation in polar coordinates:
!
!          (1/R)(d/dR)(R(dU/dR)) + (1/R**2)(d/dTHETA)(dU/dTHETA)
!
!                                + LAMBDA*U = F(R,THETA).
!
!
!
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
!       Indicates the type of boundary condition at R = A and R = B.
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
!            BDA(J) = (d/dR)U(A,THETA(J)), J = 1,2,...,N+1  .
!
!       When MBDCND has any other value, BDA is a dummy variable.
!
!     BDB
!       A one-dimensional array of length N+1 that specifies the values
!       of the derivative of the solution with respect to R at R = B.
!       When MBDCND = 2,3, or 6,
!
!            BDB(J) = (d/dR)U(B,THETA(J)), J = 1,2,...,N+1  .
!
!       When MBDCND has any other value, BDB is a dummy variable.
!
!     C,D
!       The range of THETA, i.e., C  <=  THETA  <=  D.  C must be less
!       than D.
!
!     N
!       The number of panels into which the interval (C,D) is
!       subdivided.  Hence, there will be N+1 grid points in the
!       THETA-direction given by THETA(J) = C+(J-1)DTHETA for
!       J = 1,2,...,N+1, where DTHETA = (D-C)/N is the panel width.  N
!       must be greater than 3.
!
!     NBDCND
!       Indicates the type of boundary conditions at THETA = C and
!       at THETA = D.
!
!       = 0  If the solution is periodic in THETA, i.e.,
!            U(I,J) = U(I,N+J).
!       = 1  If the solution is specified at THETA = C and THETA = D
!            (see note below).
!       = 2  If the solution is specified at THETA = C and the
!            derivative of the solution with respect to THETA is
!            specified at THETA = D (see note below).
!       = 4  If the derivative of the solution with respect to THETA is
!            specified at THETA = C and the solution is specified at
!            THETA = D (see note below).
!
!       NOTE:  When NBDCND = 1,2, or 4, do not use MBDCND = 5 or 6
!              (the former indicates that the solution is specified at
!              R = 0, the latter indicates the solution is unspecified
!              at R = 0).  Use instead MBDCND = 1 or 2  .
!
!     BDC
!       A one-dimensional array of length M+1 that specifies the values
!       of the derivative of the solution with respect to THETA at
!       THETA = C.  When NBDCND = 3 or 4,
!
!            BDC(I) = (d/dTHETA)U(R(I),C), I = 1,2,...,M+1  .
!
!       When NBDCND has any other value, BDC is a dummy variable.
!
!     BDD
!       A one-dimensional array of length M+1 that specifies the values
!       of the derivative of the solution with respect to THETA at
!       THETA = D.  When NBDCND = 2 or 3,
!
!            BDD(I) = (d/dTHETA)U(R(I),D), I = 1,2,...,M+1  .
!
!       When NBDCND has any other value, BDD is a dummy variable.
!
!     ELMBDA
!       The constant LAMBDA in the Helmholtz equation.  If
!       LAMBDA  <  0, a solution may not exist.  However, HWSPLR will
!       attempt to find a solution.
!
!     F
!       A two-dimensional array that specifies the values of the right
!       side of the Helmholtz equation and boundary values (if any).
!       For I = 2,3,...,M and J = 2,3,...,N
!
!            F(I,J) = F(R(I),THETA(J)).
!
!       On the boundaries F is defined by
!
!            MBDCND   F(1,J)            F(M+1,J)
!            ------   -------------     -------------
!
!              1      U(A,THETA(J))     U(B,THETA(J))
!              2      U(A,THETA(J))     F(B,THETA(J))
!              3      F(A,THETA(J))     F(B,THETA(J))
!              4      F(A,THETA(J))     U(B,THETA(J))   J = 1,2,...,N+1
!              5      F(0,0)            U(B,THETA(J))
!              6      F(0,0)            F(B,THETA(J))
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
!
!     IDIMF
!       The row (or first) dimension of the array F as it appears in the
!       program calling HWSPLR.  This parameter is used to specify the
!       variable dimension of F.  IDIMF must be at least M+1  .
!
!     W
!       A one-dimensional array that must be provided by the user for
!       work space.  W may require up to 4*(N+1) +
!       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of
!       locations used is computed by HWSPLR and is returned in location
!       W(1).
!
!
!             * * * * * *   On Output     * * * * * *
!
!     F
!       Contains the solution U(I,J) of the finite difference
!       approximation for the grid point (R(I),THETA(J)),
!       I = 1,2,...,M+1, J = 1,2,...,N+1  .
!
!     PERTRB
!       If a combination of periodic, derivative, or unspecified
!       boundary conditions is specified for a Poisson equation
!       (LAMBDA = 0), a solution may not exist.  PERTRB is a constant,
!       calculated and subtracted from F, which ensures that a solution
!       exists.  HWSPLR then computes this solution, which is a least
!       squares solution to the original approximation.  This solution
!       plus any constant is also a solution.  Hence, the solution is
!       not unique.  PERTRB should be small compared to the right side.
!       Otherwise, a solution is obtained to an essentially different
!       problem.  This comparison should always be made to insure that a
!       meaningful solution has been obtained.
!
!     IERROR
!       An error flag that indicates invalid input parameters.  Except
!       for numbers 0 and 11, a solution is not attempted.
!
!       =  0  No error.
!       =  1  A  <  0  .
!       =  2  A  >=  B.
!       =  3  MBDCND  <  1 or MBDCND  >  6  .
!       =  4  C  >=  D.
!       =  5  N  <=  3
!       =  6  NBDCND  <  0 or  >  4  .
!       =  7  A = 0, MBDCND = 3 or 4  .
!       =  8  A  >  0, MBDCND  >=  5  .
!       =  9  MBDCND  >=  5, NBDCND  /=  0 and NBDCND  /=  3  .
!       = 10  IDIMF  <  M+1  .
!       = 11  LAMBDA  >  0  .
!       = 12  M  <=  3
!
!       Since this is the only means of indicating a possibly incorrect
!       call to HWSPLR, the user should test IERROR after the call.
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
!     Subprograms    HWSPLR,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,
!     Required       TRIX,TRI3,PIMACH
!
!     Special        None
!     Conditions
!
!     Common         NONE
!     Blocks
!
!     I/O
!
!     Precision      Single
!
!     Specialist     Roland Sweet
!
!     Language       FORTRAN
!
!     History        Standardized April 1, 1973
!                    Revised January 1, 1976
!
!     Algorithm      The routine defines the finite difference
!                    equations, incorporates boundary data, and adjusts
!                    the right side of singular systems and then calls
!                    GENBUN to solve the system.
!
!     Space          13430(octal) = 5912(decimal)  locations on the NCAR
!     Required       Control Data 7600
!
!     Timing and        The execution time T on the NCAR Control Data
!     Accuracy       7600 for subroutine HWSPLR is roughly proportional
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
!                    Subprograms For The Solution Of Elliptic Equations'
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
!***END PROLOGUE  HWSPLR
!
!
  DIMENSION       F(IDIMF,*)
  DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     , &
                  W(*)
!***FIRST EXECUTABLE STATEMENT  HWSPLR
  IERROR = 0
  if (A  <  0.) IERROR = 1
  if (A  >=  B) IERROR = 2
  if (MBDCND <= 0 .OR. MBDCND >= 7) IERROR = 3
  if (C  >=  D) IERROR = 4
  if (N  <=  3) IERROR = 5
  if (NBDCND <= -1 .OR. NBDCND >= 5) IERROR = 6
  if (A == 0. .AND. (MBDCND == 3 .OR. MBDCND == 4)) IERROR = 7
  if (A > 0. .AND. MBDCND >= 5) IERROR = 8
  if (MBDCND >= 5 .AND. NBDCND /= 0 .AND. NBDCND /= 3) IERROR = 9
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
  MSTOP = MP1
  go to (101,105,102,103,104,105),MBDCND
  101 MSTOP = M
  go to 105
  102 MSTART = 1
  go to 105
  103 MSTART = 1
  104 MSTOP = M
  105 MUNK = MSTOP-MSTART+1
  NSTART = 1
  NSTOP = N
  go to (109,106,107,108,109),NP
  106 NSTART = 2
  go to 109
  107 NSTART = 2
  108 NSTOP = NP1
  109 NUNK = NSTOP-NSTART+1
!
!     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
!
  ID2 = MUNK
  ID3 = ID2+MUNK
  ID4 = ID3+MUNK
  ID5 = ID4+MUNK
  ID6 = ID5+MUNK
  A1 = 2./DLRSQ
  IJ = 0
  if (MBDCND == 3 .OR. MBDCND == 4) IJ = 1
  DO 110 I=1,MUNK
     R = A+(I-IJ)*DELTAR
     J = ID5+I
     W(J) = R
     J = ID6+I
     W(J) = 1./R**2
     W(I) = (R-DLRBY2)/(R*DLRSQ)
     J = ID3+I
     W(J) = (R+DLRBY2)/(R*DLRSQ)
     J = ID2+I
     W(J) = -A1+ELMBDA
  110 CONTINUE
  go to (114,111,112,113,114,111),MBDCND
  111 W(ID2) = A1
  go to 114
  112 W(ID2) = A1
  113 W(ID3+1) = A1
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
!     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
!
  124 A1 = 1./DLTHSQ
  L = ID5-MSTART+1
  LP = ID6-MSTART+1
  go to (134,125,125,127,127),NP
  125 DO 126 I=MSTART,MSTOP
     J = I+LP
     F(I,2) = F(I,2)-A1*W(J)*F(I,1)
  126 CONTINUE
  go to 129
  127 A1 = 2./DELTHT
  DO 128 I=MSTART,MSTOP
     J = I+LP
     F(I,1) = F(I,1)+A1*W(J)*BDC(I)
  128 CONTINUE
  129 A1 = 1./DLTHSQ
  go to (134,130,132,132,130),NP
  130 DO 131 I=MSTART,MSTOP
     J = I+LP
     F(I,N) = F(I,N)-A1*W(J)*F(I,NP1)
  131 CONTINUE
  go to 134
  132 A1 = 2./DELTHT
  DO 133 I=MSTART,MSTOP
     J = I+LP
     F(I,NP1) = F(I,NP1)-A1*W(J)*BDD(I)
  133 CONTINUE
  134 CONTINUE
!
!     ADJUST RIGHT SIDE OF EQUATION FOR UNKNOWN AT POLE WHEN HAVE
!     DERIVATIVE SPECIFIED BOUNDARY CONDITIONS.
!
  if (MBDCND >= 5 .AND. NBDCND == 3) &
      F(1,1) = F(1,1)-(BDD(2)-BDC(2))*4./(N*DELTHT*DLRSQ)
!
!     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
!     SOLUTION.
!
  PERTRB = 0.
  if (ELMBDA) 144,136,135
  135 IERROR = 11
  go to 144
  136 if (NBDCND /= 0 .AND. NBDCND /= 3) go to 144
  S2 = 0.
  go to (144,144,137,144,144,138),MBDCND
  137 W(ID5+1) = .5*(W(ID5+2)-DLRBY2)
  S2 = .25*DELTAR
  138 A2 = 2.
  if (NBDCND  ==  0) A2 = 1.
  J = ID5+MUNK
  W(J) = .5*(W(J-1)+DLRBY2)
  S = 0.
  DO 140 I=MSTART,MSTOP
     S1 = 0.
     IJ = NSTART+1
     K = NSTOP-1
     DO 139 J=IJ,K
        S1 = S1+F(I,J)
  139    CONTINUE
     J = I+L
     S = S+(A2*S1+F(I,NSTART)+F(I,NSTOP))*W(J)
  140 CONTINUE
  S2 = M*A+DELTAR*((M-1)*(M+1)*.5+.25)+S2
  S1 = (2.+A2*(NUNK-2))*S2
  if (MBDCND  ==  3) go to 141
  S2 = N*A2*DELTAR/8.
  S = S+F(1,1)*S2
  S1 = S1+S2
  141 CONTINUE
  PERTRB = S/S1
  DO 143 I=MSTART,MSTOP
     DO 142 J=NSTART,NSTOP
        F(I,J) = F(I,J)-PERTRB
  142    CONTINUE
  143 CONTINUE
  144 CONTINUE
!
!     MULTIPLY I-TH EQUATION THROUGH BY (R(I)*DELTHT)**2.
!
  DO 146 I=MSTART,MSTOP
     K = I-MSTART+1
     J = I+LP
     A1 = DLTHSQ/W(J)
     W(K) = A1*W(K)
     J = ID2+K
     W(J) = A1*W(J)
     J = ID3+K
     W(J) = A1*W(J)
     DO 145 J=NSTART,NSTOP
        F(I,J) = A1*F(I,J)
  145    CONTINUE
  146 CONTINUE
  W(1) = 0.
  W(ID4) = 0.
!
!     call GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
!
  call GENBUN (NBDCND,NUNK,1,MUNK,W(1),W(ID2+1),W(ID3+1),IDIMF, &
               F(MSTART,NSTART),IERR1,W(ID4+1))
  IWSTOR = W(ID4+1)+3*MUNK
  go to (157,157,157,157,148,147),MBDCND
!
!     ADJUST THE SOLUTION AS NECESSARY FOR THE PROBLEMS WHERE A = 0.
!
  147 if (ELMBDA  /=  0.) go to 148
  YPOLE = 0.
  go to 155
  148 CONTINUE
  J = ID5+MUNK
  W(J) = W(ID2)/W(ID3)
  DO 149 IP=3,MUNK
     I = MUNK-IP+2
     J = ID5+I
     LP = ID2+I
     K = ID3+I
     W(J) = W(I)/(W(LP)-W(K)*W(J+1))
  149 CONTINUE
  W(ID5+1) = -.5*DLTHSQ/(W(ID2+1)-W(ID3+1)*W(ID5+2))
  DO 150 I=2,MUNK
     J = ID5+I
     W(J) = -W(J)*W(J-1)
  150 CONTINUE
  S = 0.
  DO 151 J=NSTART,NSTOP
     S = S+F(2,J)
  151 CONTINUE
  A2 = NUNK
  if (NBDCND  ==  0) go to 152
  S = S-.5*(F(2,NSTART)+F(2,NSTOP))
  A2 = A2-1.
  152 YPOLE = (.25*DLRSQ*F(1,1)-S/A2)/(W(ID5+1)-1.+ELMBDA*DLRSQ*.25)
  DO 154 I=MSTART,MSTOP
     K = L+I
     DO 153 J=NSTART,NSTOP
        F(I,J) = F(I,J)+YPOLE*W(K)
  153    CONTINUE
  154 CONTINUE
  155 DO 156 J=1,NP1
     F(1,J) = YPOLE
  156 CONTINUE
  157 CONTINUE
  if (NBDCND  /=  0) go to 159
  DO 158 I=MSTART,MSTOP
     F(I,NP1) = F(I,1)
  158 CONTINUE
  159 CONTINUE
  W(1) = IWSTOR
  return
end
