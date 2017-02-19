subroutine GENBUN (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, &
     IERROR, W)
!
!! GENBUN solves by a cyclic reduction algorithm the linear system ...
!  of equations that results from a finite difference approximation to
!  certain 2-d elliptic PDE's on a centered grid.
!
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B4B
!***TYPE      SINGLE PRECISION (GENBUN-S, CMGNBN-C)
!***KEYWORDS  ELLIPTIC, FISHPACK, PDE, TRIDIAGONAL
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!     Subroutine GENBUN solves the linear system of equations
!
!          A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
!
!          + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
!
!               for I = 1,2,...,M  and  J = 1,2,...,N.
!
!     The indices I+1 and I-1 are evaluated modulo M, i.e.,
!     X(0,J) = X(M,J) and X(M+1,J) = X(1,J), and X(I,0) may be equal to
!     0, X(I,2), or X(I,N) and X(I,N+1) may be equal to 0, X(I,N-1), or
!     X(I,1) depending on an input parameter.
!
!
!     * * * * * * * *    Parameter Description     * * * * * * * * * *
!
!             * * * * * *   On Input    * * * * * *
!
!     NPEROD
!       Indicates the values that X(I,0) and X(I,N+1) are assumed to
!       have.
!
!       = 0  If X(I,0) = X(I,N) and X(I,N+1) = X(I,1).
!       = 1  If X(I,0) = X(I,N+1) = 0  .
!       = 2  If X(I,0) = 0 and X(I,N+1) = X(I,N-1).
!       = 3  If X(I,0) = X(I,2) and X(I,N+1) = X(I,N-1).
!       = 4  If X(I,0) = X(I,2) and X(I,N+1) = 0.
!
!     N
!       The number of unknowns in the J-direction.  N must be greater
!       than 2.
!
!     MPEROD
!       = 0 if A(1) and C(M) are not zero.
!       = 1 if A(1) = C(M) = 0.
!
!     M
!       The number of unknowns in the I-direction.  M must be greater
!       than 2.
!
!     A,B,C
!       One-dimensional arrays of length M that specify the
!       coefficients in the linear equations given above.  If MPEROD = 0
!       the array elements must not depend upon the index I, but must be
!       constant.  Specifically, the subroutine checks the following
!       condition
!
!             A(I) = C(1)
!             C(I) = C(1)
!             B(I) = B(1)
!
!       for I=1,2,...,M.
!
!     IDIMY
!       The row (or first) dimension of the two-dimensional array Y as
!       it appears in the program calling GENBUN.  This parameter is
!       used to specify the variable dimension of Y.  IDIMY must be at
!       least M.
!
!     Y
!       A two-dimensional array that specifies the values of the right
!       side of the linear system of equations given above.  Y must be
!       dimensioned at least M*N.
!
!     W
!       A one-dimensional array that must be provided by the user for
!       work space.  W may require up to 4*N + (10 + INT(log2(N)))*M
!       locations.  The actual number of locations used is computed by
!       GENBUN and is returned in location W(1).
!
!
!             * * * * * *   On Output     * * * * * *
!
!     Y
!       Contains the solution X.
!
!     IERROR
!       An error flag that indicates invalid input parameters.  Except
!       for number zero, a solution is not attempted.
!
!       = 0  No error.
!       = 1  M  <=  2
!       = 2  N  <=  2
!       = 3  IDIMY  <  M
!       = 4  NPEROD  <  0 or NPEROD  >  4
!       = 5  MPEROD  <  0 or MPEROD  >  1
!       = 6  A(I)  /=  C(1) or C(I)  /=  C(1) or B(I)  /=  B(1) for
!            some I=1,2,...,M.
!       = 7  A(1)  /=  0 or C(M)  /=  0 and MPEROD = 1
!
!     W
!       W(1) contains the required length of W.
!
! *Long Description:
!
!     * * * * * * *   Program Specifications    * * * * * * * * * * * *
!
!     Dimension of   A(M),B(M),C(M),Y(IDIMY,N),W(see parameter list)
!     Arguments
!
!     Latest         June 1, 1976
!     Revision
!
!     Subprograms    GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,TRIX,TRI3,
!     Required       PIMACH
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
!     History        Standardized April 1, 1973
!                    Revised August 20,1973
!                    Revised January 1, 1976
!
!     Algorithm      The linear system is solved by a cyclic reduction
!                    algorithm described in the reference.
!
!     Space          4944(decimal) = 11520(octal) locations on the NCAR
!     Required       Control Data 7600.
!
!     Timing and        The execution time T on the NCAR Control Data
!     Accuracy       7600 for subroutine GENBUN is roughly proportional
!                    to M*N*log2(N), but also depends on the input
!                    parameter NPEROD.  Some typical values are listed
!                    in the table below.  More comprehensive timing
!                    charts may be found in the reference.
!                       To measure the accuracy of the algorithm a
!                    uniform random number generator was used to create
!                    a solution array X for the system given in the
!                    'PURPOSE' with
!
!                       A(I) = C(I) = -0.5*B(I) = 1,       I=1,2,...,M
!
!                    and, when MPEROD = 1
!
!                       A(1) = C(M) = 0
!                       A(M) = C(1) = 2.
!
!                    The solution X was substituted into the given sys-
!                    tem and, using double precision, a right side Y was
!                    computed.  Using this array Y subroutine GENBUN was
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
!                         31        0         0          36     6.E-14
!                         31        1         1          21     4.E-13
!                         31        1         3          41     3.E-13
!                         32        0         0          29     9.E-14
!                         32        1         1          32     3.E-13
!                         32        1         3          48     1.E-13
!                         33        0         0          36     9.E-14
!                         33        1         1          30     4.E-13
!                         33        1         3          34     1.E-13
!                         63        0         0         150     1.E-13
!                         63        1         1          91     1.E-12
!                         63        1         3         173     2.E-13
!                         64        0         0         122     1.E-13
!                         64        1         1         128     1.E-12
!                         64        1         3         199     6.E-13
!                         65        0         0         143     2.E-13
!                         65        1         1         120     1.E-12
!                         65        1         3         138     4.E-13
!
!     Portability    American National Standards Institute Fortran.
!                    The machine dependent constant PI is defined in
!                    function PIMACH.
!
!     Required       COS
!     Resident
!     Routines
!
!     Reference      Sweet, R., 'A Cyclic Reduction Algorithm For
!                    Solving Block Tridiagonal Systems Of Arbitrary
!                    Dimensions,' SIAM J. on Numer. Anal.,
!                    14(Sept., 1977), PP. 706-720.
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!***REFERENCES  R. Sweet, A cyclic reduction algorithm for solving
!                 block tridiagonal systems of arbitrary dimensions,
!                 SIAM Journal on Numerical Analysis 14, (September
!                 1977), pp. 706-720.
!***ROUTINES CALLED  POISD2, POISN2, POISP2
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  GENBUN
!
!
  DIMENSION       Y(IDIMY,*)
  DIMENSION       W(*)       ,B(*)       ,A(*)       ,C(*)
!***FIRST EXECUTABLE STATEMENT  GENBUN
  IERROR = 0
  if (M  <=  2) IERROR = 1
  if (N  <=  2) IERROR = 2
  if (IDIMY  <  M) IERROR = 3
  if (NPEROD < 0 .OR. NPEROD > 4) IERROR = 4
  if (MPEROD < 0 .OR. MPEROD > 1) IERROR = 5
  if (MPEROD  ==  1) go to 102

  DO I=2,M
     if (A(I)  /=  C(1)) go to 103
     if (C(I)  /=  C(1)) go to 103
     if (B(I)  /=  B(1)) go to 103
  end do

  go to 104
  102 if (A(1) /= 0. .OR. C(M) /= 0.) IERROR = 7
  go to 104
  103 IERROR = 6
  104 if (IERROR  /=  0) RETURN
  MP1 = M+1
  IWBA = MP1
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
  MP = MPEROD+1
  NP = NPEROD+1
  go to (114,107),MP
  107 go to (108,109,110,111,123),NP
  108 call POISP2 (M,N,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2), &
               W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS), &
               W(IWP))
  go to 112
  109 call POISD2 (M,N,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWW1), &
               W(IWD),W(IWTCOS),W(IWP))
  go to 112
  110 call POISN2 (M,N,1,2,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2), &
               W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS), &
               W(IWP))
  go to 112
  111 call POISN2 (M,N,1,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2), &
               W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS), &
               W(IWP))
  112 IPSTOR = W(IWW1)
  IREV = 2
  if (NPEROD  ==  4) go to 124
  113 go to (127,133),MP
  114 CONTINUE
!
!     REORDER UNKNOWNS WHEN MP =0
!
  MH = (M+1)/2
  MHM1 = MH-1
  MODD = 1
  if (MH*2  ==  M) MODD = 2
  DO 119 J=1,N
     DO 115 I=1,MHM1
        MHPI = MH+I
        MHMI = MH-I
        W(I) = Y(MHMI,J)-Y(MHPI,J)
        W(MHPI) = Y(MHMI,J)+Y(MHPI,J)
  115    CONTINUE
     W(MH) = 2.*Y(MH,J)
     go to (117,116),MODD
  116    W(M) = 2.*Y(M,J)
  117    CONTINUE
     DO 118 I=1,M
        Y(I,J) = W(I)
  118    CONTINUE
  119 CONTINUE
  K = IWBC+MHM1-1
  I = IWBA+MHM1
  W(K) = 0.
  W(I) = 0.
  W(K+1) = 2.*W(K+1)
  go to (120,121),MODD
  120 CONTINUE
  K = IWBB+MHM1-1
  W(K) = W(K)-W(I-1)
  W(IWBC-1) = W(IWBC-1)+W(IWBB-1)
  go to 122
  121 W(IWBB-1) = W(K+1)
  122 CONTINUE
  go to 107
!
!     REVERSE COLUMNS WHEN NPEROD = 4.
!
  123 IREV = 1
  NBY2 = N/2
  124 DO 126 J=1,NBY2
     MSKIP = N+1-J
     DO 125 I=1,M
        A1 = Y(I,J)
        Y(I,J) = Y(I,MSKIP)
        Y(I,MSKIP) = A1
  125    CONTINUE
  126 CONTINUE
  go to (110,113),IREV
  127 CONTINUE
  DO 132 J=1,N
     DO 128 I=1,MHM1
        MHMI = MH-I
        MHPI = MH+I
        W(MHMI) = .5*(Y(MHPI,J)+Y(I,J))
        W(MHPI) = .5*(Y(MHPI,J)-Y(I,J))
  128    CONTINUE
     W(MH) = .5*Y(MH,J)
     go to (130,129),MODD
  129    W(M) = .5*Y(M,J)
  130    CONTINUE
     DO 131 I=1,M
        Y(I,J) = W(I)
  131    CONTINUE
  132 CONTINUE
  133 CONTINUE
!
!     return STORAGE REQUIREMENTS FOR W ARRAY.
!
  W(1) = IPSTOR+IWP-1
  return
end
