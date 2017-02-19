subroutine CMGNBN (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, &
     IERROR, W)
!
!! CMGNBN solves a complex block tridiagonal linear system of ...
!  equations by a cyclic reduction algorithm.
!
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B4B
!***TYPE      COMPLEX (GENBUN-S, CMGNBN-C)
!***KEYWORDS  CYCLIC REDUCTION, ELLIPTIC PDE, FISHPACK,
!             TRIDIAGONAL LINEAR SYSTEM
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!     Subroutine CMGNBN solves the complex linear system of equations
!
!          A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
!
!          + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
!
!               For I = 1,2,...,M  and  J = 1,2,...,N.
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
!       = 0 If A(1) and C(M) are not zero
!       = 1 If A(1) = C(M) = 0
!
!     M
!       The number of unknowns in the I-direction.  N must be greater
!       than 2.
!
!     A,B,C
!       One-dimensional complex arrays of length M that specify the
!       coefficients in the linear equations given above.  If MPEROD = 0
!       the array elements must not depend upon the index I, but must be
!       constant.  Specifically, the subroutine checks the following
!       condition
!
!             A(I) = C(1)
!             C(I) = C(1)
!             B(I) = B(1)
!
!       For I=1,2,...,M.
!
!     IDIMY
!       The row (or first) dimension of the two-dimensional array Y as
!       it appears in the program calling CMGNBN.  This parameter is
!       used to specify the variable dimension of Y.  IDIMY must be at
!       least M.
!
!     Y
!       A two-dimensional complex array that specifies the values of the
!       right side of the linear system of equations given above.  Y
!       must be dimensioned at least M*N.
!
!     W
!       A one-dimensional complex array that must be provided by the
!       user for work space.  W may require up to 4*N +
!       (10 + INT(log2(N)))*M LOCATIONS.  The actual number of locations
!       used is computed by CMGNBN and is returned in location W(1).
!
!
!             * * * * * *   On Output     * * * * * *
!
!     Y
!       Contains the solution X.
!
!     IERROR
!       An error flag which indicates invalid input parameters.  Except
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
!     Latest         June 1979
!     Revision
!
!     Subprograms    CMGNBN,CMPOSD,CMPOSN,CMPOSP,CMPCSG,CMPMRG,
!     Required       CMPTRX,CMPTR3,PIMACH
!
!     Special        None
!     Conditions
!
!     Common         None
!     Blocks
!
!     I/O            None
!
!     Precision      Single
!
!     Specialist     Roland Sweet
!
!     Language       FORTRAN
!
!     History        Written by Roland Sweet at NCAR in June, 1977
!
!     Algorithm      The linear system is solved by a cyclic reduction
!                    algorithm described in the reference.
!
!     Space          4944(DECIMAL) = 11520(octal) locations on the NCAR
!     Required       Control Data 7600
!
!     Timing and      The execution time T on the NCAR Control Data
!     Accuracy       7600 for subroutine CMGNBN is roughly proportional
!                    to M*N*log2(N), but also depends on the input
!                    parameter NPEROD.  Some typical values are listed
!                    in the table below.
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
!                    tem and a right side Y was computed.  Using this
!                    array Y subroutine CMGNBN was called to produce an
!                    approximate solution Z.  Then the relative error,
!                    defined as
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
!                         31        0         0          77     1.E-12
!                         31        1         1          45     4.E-13
!                         31        1         3          91     2.E-12
!                         32        0         0          59     7.E-14
!                         32        1         1          65     5.E-13
!                         32        1         3          97     2.E-13
!                         33        0         0          80     6.E-13
!                         33        1         1          67     5.E-13
!                         33        1         3          76     3.E-12
!                         63        0         0         350     5.E-12
!                         63        1         1         215     6.E-13
!                         63        1         3         412     1.E-11
!                         64        0         0         264     1.E-13
!                         64        1         1         287     3.E-12
!                         64        1         3         421     3.E-13
!                         65        0         0         338     2.E-12
!                         65        1         1         292     5.E-13
!                         65        1         3         329     1.E-11
!
!     Portability    American National Standards Institute Fortran.
!                    The machine dependent constant PI is defined in
!                    function PIMACH.
!
!     Required       COS
!     Resident
!     Routines
!
!     Reference      Sweet, R., 'A Cyclic Reduction Algorithm for
!                    Solving Block Tridiagonal Systems Of Arbitrary
!                    Dimensions,' SIAM J. on Numer. Anal.,
!                    14(SEPT., 1977), PP. 706-720.
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!***REFERENCES  R. Sweet, A cyclic reduction algorithm for solving
!                 block tridiagonal systems of arbitrary dimensions,
!                 SIAM Journal on Numerical Analysis 14, (September
!                 1977), pp. 706-720.
!***ROUTINES CALLED  CMPOSD, CMPOSN, CMPOSP
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CMGNBN
!
!
  COMPLEX         A          ,B          ,C          ,Y          , &
                  W          ,A1
  DIMENSION       Y(IDIMY,*)
  DIMENSION       W(*)       ,B(*)       ,A(*)       ,C(*)
!***FIRST EXECUTABLE STATEMENT  CMGNBN
  IERROR = 0
  if (M  <=  2) IERROR = 1
  if (N  <=  2) IERROR = 2
  if (IDIMY  <  M) IERROR = 3
  if (NPEROD < 0 .OR. NPEROD > 4) IERROR = 4
  if (MPEROD < 0 .OR. MPEROD > 1) IERROR = 5
  if (MPEROD  ==  1) go to 102
  DO 101 I=2,M
     if (ABS(A(I)-C(1))  /=  0.) go to 103
     if (ABS(C(I)-C(1))  /=  0.) go to 103
     if (ABS(B(I)-B(1))  /=  0.) go to 103
  101 CONTINUE
  go to 104
  102 if (ABS(A(1)) /= 0. .AND. ABS(C(M)) /= 0.) IERROR = 7
  go to 104
  103 IERROR = 6
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
  MP = MPEROD+1
  NP = NPEROD+1
  go to (114,107),MP
  107 go to (108,109,110,111,123),NP
  108 call CMPOSP (M,N,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2), &
               W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS), &
               W(IWP))
  go to 112
  109 call CMPOSD (M,N,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWW1), &
               W(IWD),W(IWTCOS),W(IWP))
  go to 112
  110 call CMPOSN (M,N,1,2,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2), &
               W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS), &
               W(IWP))
  go to 112
  111 call CMPOSN (M,N,1,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2), &
               W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS), &
               W(IWP))
  112 IPSTOR = REAL(W(IWW1))
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
  W(K) = (0.,0.)
  W(I) = (0.,0.)
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
!     REVERSE COLUMNS WHEN NPEROD = 4
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
  W(1) = CMPLX(REAL(IPSTOR+IWP-1),0.)
  return
end
