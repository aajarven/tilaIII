subroutine POISTG (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, &
     IERROR, W)
!
!! POISTG solves a block tridiagonal system of linear equations ...
!            that results from a staggered grid finite difference
!            approximation to 2-D elliptic PDE's.
!
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B4B
!***TYPE      SINGLE PRECISION (POISTG-S)
!***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, TRIDIAGONAL
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!     Subroutine POISTG solves the linear system of equations
!
!       A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
!       + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
!
!       for I=1,2,...,M and J=1,2,...,N.
!
!     The indices I+1 and I-1 are evaluated modulo M, i.e.
!     X(0,J) = X(M,J) and X(M+1,J) = X(1,J), and X(I,0) may be equal to
!     X(I,1) or -X(I,1) and X(I,N+1) may be equal to X(I,N) or -X(I,N)
!     depending on an input parameter.
!
!
!     * * * * * * * *    Parameter Description     * * * * * * * * * *
!
!             * * * * * *   On Input    * * * * * *
!
!   NPEROD
!     Indicates the values which X(I,0) and X(I,N+1) are assumed
!     to have.
!     = 1 If X(I,0) = -X(I,1) and X(I,N+1) = -X(I,N)
!     = 2 If X(I,0) = -X(I,1) and X(I,N+1) =  X(I,N)
!     = 3 If X(I,0) =  X(I,1) and X(I,N+1) =  X(I,N)
!     = 4 If X(I,0) =  X(I,1) and X(I,N+1) = -X(I,N)
!
!   N
!     The number of unknowns in the J-direction.  N must
!     be greater than 2.
!
!   MPEROD
!     = 0 If A(1) and C(M) are not zero
!     = 1 If A(1) = C(M) = 0
!
!   M
!     The number of unknowns in the I-direction.  M must
!     be greater than 2.
!
!   A,B,C
!     One-dimensional arrays of length M that specify the coefficients
!     in the linear equations given above.  If MPEROD = 0 the array
!     elements must not depend on the index I, but must be constant.
!     Specifically, the subroutine checks the following condition
!
!           A(I) = C(1)
!           B(I) = B(1)
!           C(I) = C(1)
!
!     for I = 1, 2, ..., M.
!
!   IDIMY
!     The row (or first) dimension of the two-dimensional array Y as
!     it appears in the program calling POISTG.  This parameter is
!     used to specify the variable dimension of Y.  IDIMY must be at
!     least M.
!
!   Y
!     A two-dimensional array that specifies the values of the
!     right side of the linear system of equations given above.
!     Y must be dimensioned at least M X N.
!
!   W
!     A one-dimensional work array that must be provided by the user
!     for work space.  W may require up to 9M + 4N + M(INT(log2(N)))
!     locations.  The actual number of locations used is computed by
!     POISTG and returned in location W(1).
!
!
!             * * * * * *   On Output     * * * * * *
!
!   Y
!     Contains the solution X.
!
!   IERROR
!     An error flag that indicates invalid input parameters.  Except
!     for number zero, a solution is not attempted.
!     = 0  No error
!     = 1  If M  <=  2
!     = 2  If N  <=  2
!     = 3  IDIMY  <  M
!     = 4  If NPEROD  <  1 or NPEROD  >  4
!     = 5  If MPEROD  <  0 or MPEROD  >  1
!     = 6  If MPEROD = 0 and
!          A(I)  /=  C(1) or B(I)  /=  B(1) or C(I)  /=  C(1)
!          for some I = 1, 2, ..., M.
!       = 7 If MPEROD  ==  1 .AND. (A(1) /= 0 .OR. C(M) /= 0)
!
!   W
!     W(1) contains the required length of W.
!
! *Long Description:
!
!     * * * * * * *   Program Specifications    * * * * * * * * * * * *
!
!     Dimension of   A(M),B(M),C(M),Y(IDIMY,N),
!     Arguments      W(see argument list)
!
!     Latest         June 1, 1977
!     Revision
!
!     Subprograms    POISTG,POSTG2,COSGEN,MERGE,TRIX,TRI3,PIMACH
!     Required
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
!     History        Written by Roland Sweet in 1973
!                    Revised by Roland Sweet in 1977
!
!
!     Space          3297(decimal) = 6341(octal) locations on the
!     Required       NCAR Control Data 7600
!
!     Timing and        The execution time T on the NCAR Control Data
!     Accuracy       7600 for subroutine POISTG is roughly proportional
!                    to M*N*log2(N).  Some typical values are listed
!                    in the table below.  More comprehensive timing
!                    charts may be found in the reference.
!                       To measure the accuracy of the algorithm a
!                    uniform random number generator was used to create
!                    a solution array X for the system given in the
!                    'PURPOSE ' with
!
!                       A(I) = C(I) = -0.5*B(I) = 1,       I=1,2,...,M
!
!                    and, when MPEROD = 1
!
!                       A(1) = C(M) = 0
!                       B(1) = B(M) =-1.
!
!                    The solution X was substituted into the given sys-
!                    tem and, using double precision, a right side Y was
!                    computed.  Using this array Y subroutine POISTG was
!                    called to produce an approximate solution Z.  Then
!                    the relative error, defined as
!
!                       E = MAX(ABS(Z(I,J)-X(I,J)))/MAX(ABS(X(I,J)))
!
!                    where the two maxima are taken over all I=1,2,...,M
!                    and J=1,2,...,N, was computed.  The value of E is
!                    given in the table below for some typical values of
!                    M and N.
!
!
!                       M (=N)    MPEROD    NPEROD    T(MSECS)    E
!                       ------    ------    ------    --------  ------
!
!                         31        0-1       1-4        45     9.E-13
!                         31        1         1          21     4.E-13
!                         31        1         3          41     3.E-13
!                         32        0-1       1-4        51     3.E-12
!                         32        1         1          32     3.E-13
!                         32        1         3          48     1.E-13
!                         33        0-1       1-4        42     1.E-12
!                         33        1         1          30     4.E-13
!                         33        1         3          34     1.E-13
!                         63        0-1       1-4       186     3.E-12
!                         63        1         1          91     1.E-12
!                         63        1         3         173     2.E-13
!                         64        0-1       1-4       209     4.E-12
!                         64        1         1         128     1.E-12
!                         64        1         3         199     6.E-13
!                         65        0-1       1-4       143     2.E-13
!                         65        1         1         160     1.E-11
!                         65        1         3         138     4.E-13
!
!     Portability    American National Standards Institute FORTRAN.
!                    The machine dependent constant PI is defined in
!                    function PIMACH.
!
!     Required       COS
!     Resident
!     Routines
!
!     Reference      Schumann, U. and R. Sweet,'A Direct Method for
!                    the Solution of Poisson's Equation With Neumann
!                    Boundary Conditions on a Staggered Grid of
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
!***ROUTINES CALLED  POSTG2
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  POISTG
!
!
  DIMENSION       Y(IDIMY,*)
  DIMENSION       W(*)       ,B(*)       ,A(*)       ,C(*)
!***FIRST EXECUTABLE STATEMENT  POISTG
  IERROR = 0
  if (M  <=  2) IERROR = 1
  if (N  <=  2) IERROR = 2
  if (IDIMY  <  M) IERROR = 3
  if (NPEROD < 1 .OR. NPEROD > 4) IERROR = 4
  if (MPEROD < 0 .OR. MPEROD > 1) IERROR = 5
  if (MPEROD  ==  1) go to 103
  DO 101 I=1,M
     if (A(I)  /=  C(1)) go to 102
     if (C(I)  /=  C(1)) go to 102
     if (B(I)  /=  B(1)) go to 102
  101 CONTINUE
  go to 104
  102 IERROR = 6
  return
  103 if (A(1) /= 0. .OR. C(M) /= 0.) IERROR = 7
  104 if (IERROR  /=  0) RETURN
  IWBA = M+1
  IWBB = IWBA+M
  IWBC = IWBB+M
  IWB2 = IWBC+M
  IWB3 = IWB2+M
  IWW1 = IWB3+M
  IWW2 = IWW1+M
  IWW3 = IWW2+M
  IWD = IWW3+M
  IWTCOS = IWD+M
  IWP = IWTCOS+4*N
  DO 106 I=1,M
     K = IWBA+I-1
     W(K) = -A(I)
     K = IWBC+I-1
     W(K) = -C(I)
     K = IWBB+I-1
     W(K) = 2.-B(I)
     DO 105 J=1,N
        Y(I,J) = -Y(I,J)
  105    CONTINUE
  106 CONTINUE
  NP = NPEROD
  MP = MPEROD+1
  go to (110,107),MP
  107 CONTINUE
  go to (108,108,108,119),NPEROD
  108 CONTINUE
  call POSTG2 (NP,N,M,W(IWBA),W(IWBB),W(IWBC),IDIMY,Y,W,W(IWB2), &
               W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS), &
               W(IWP))
  IPSTOR = W(IWW1)
  IREV = 2
  if (NPEROD  ==  4) go to 120
  109 CONTINUE
  go to (123,129),MP
  110 CONTINUE
!
!     REORDER UNKNOWNS WHEN MP =0
!
  MH = (M+1)/2
  MHM1 = MH-1
  MODD = 1
  if (MH*2  ==  M) MODD = 2
  DO 115 J=1,N
     DO 111 I=1,MHM1
        MHPI = MH+I
        MHMI = MH-I
        W(I) = Y(MHMI,J)-Y(MHPI,J)
        W(MHPI) = Y(MHMI,J)+Y(MHPI,J)
  111    CONTINUE
     W(MH) = 2.*Y(MH,J)
     go to (113,112),MODD
  112    W(M) = 2.*Y(M,J)
  113    CONTINUE
     DO 114 I=1,M
        Y(I,J) = W(I)
  114    CONTINUE
  115 CONTINUE
  K = IWBC+MHM1-1
  I = IWBA+MHM1
  W(K) = 0.
  W(I) = 0.
  W(K+1) = 2.*W(K+1)
  go to (116,117),MODD
  116 CONTINUE
  K = IWBB+MHM1-1
  W(K) = W(K)-W(I-1)
  W(IWBC-1) = W(IWBC-1)+W(IWBB-1)
  go to 118
  117 W(IWBB-1) = W(K+1)
  118 CONTINUE
  go to 107
  119 CONTINUE
!
!     REVERSE COLUMNS WHEN NPEROD = 4.
!
  IREV = 1
  NBY2 = N/2
  NP = 2
  120 DO 122 J=1,NBY2
     MSKIP = N+1-J
     DO 121 I=1,M
        A1 = Y(I,J)
        Y(I,J) = Y(I,MSKIP)
        Y(I,MSKIP) = A1
  121    CONTINUE
  122 CONTINUE
  go to (108,109),IREV
  123 CONTINUE
  DO 128 J=1,N
     DO 124 I=1,MHM1
        MHMI = MH-I
        MHPI = MH+I
        W(MHMI) = .5*(Y(MHPI,J)+Y(I,J))
        W(MHPI) = .5*(Y(MHPI,J)-Y(I,J))
  124    CONTINUE
     W(MH) = .5*Y(MH,J)
     go to (126,125),MODD
  125    W(M) = .5*Y(M,J)
  126    CONTINUE
     DO 127 I=1,M
        Y(I,J) = W(I)
  127    CONTINUE
  128 CONTINUE
  129 CONTINUE
!
!     return STORAGE REQUIREMENTS FOR W ARRAY.
!
  W(1) = IPSTOR+IWP-1
  return
end
