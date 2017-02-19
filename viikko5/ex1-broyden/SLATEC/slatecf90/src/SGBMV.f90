subroutine SGBMV (TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX, &
     BETA, Y, INCY)
!
!! SGBMV multiplies a real vector by a real general band matrix.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B4
!***TYPE      SINGLE PRECISION (SGBMV-S, DGBMV-D, CGBMV-C)
!***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J. J., (ANL)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!  SGBMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n band matrix, with kl sub-diagonals and ku super-diagonals.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  KL     - INTEGER.
!           On entry, KL specifies the number of sub-diagonals of the
!           matrix A. KL must satisfy  0 .le. KL.
!           Unchanged on exit.
!
!  KU     - INTEGER.
!           On entry, KU specifies the number of super-diagonals of the
!           matrix A. KU must satisfy  0 .le. KU.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry, the leading ( kl + ku + 1 ) by n part of the
!           array A must contain the matrix of coefficients, supplied
!           column by column, with the leading diagonal of the matrix in
!           row ( ku + 1 ) of the array, the first super-diagonal
!           starting at position 2 in row ku, the first sub-diagonal
!           starting at position 1 in row ( ku + 2 ), and so on.
!           Elements in the array A that do not correspond to elements
!           in the band matrix (such as the top left ku by ku triangle)
!           are not referenced.
!           The following program segment will transfer a band matrix
!           from conventional full matrix storage to band storage:
!
!                 DO 20, J = 1, N
!                    K = KU + 1 - J
!                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )
!                       A( K + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( kl + ku + 1 ).
!           Unchanged on exit.
!
!  X      - REAL             array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - REAL            .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - REAL             array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry, the incremented array Y must contain the
!           vector y. On exit, Y is overwritten by the updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!***REFERENCES  Dongarra, J. J., Du Croz, J., Hammarling, S., and
!                 Hanson, R. J.  An extended set of Fortran basic linear
!                 algebra subprograms.  ACM TOMS, Vol. 14, No. 1,
!                 pp. 1-17, March 1988.
!***ROUTINES CALLED  LSAME, XERBLA
!***REVISION HISTORY  (YYMMDD)
!   861022  DATE WRITTEN
!   910605  Modified to meet SLATEC prologue standards.  Only comment
!           lines were modified.  (BKS)
!***END PROLOGUE  SGBMV
!     .. Scalar Arguments ..
  REAL               ALPHA, BETA
  INTEGER            INCX, INCY, KL, KU, LDA, M, N
  CHARACTER*1        TRANS
!     .. Array Arguments ..
  REAL               A( LDA, * ), X( * ), Y( * )
  REAL               ONE         , ZERO
  PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     .. Local Scalars ..
  REAL               TEMP
  INTEGER            I, INFO, IX, IY, J, JX, JY, K, KUP1, KX, KY, &
                     LENX, LENY
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!***FIRST EXECUTABLE STATEMENT  SGBMV
!
!     Test the input parameters.
!
  INFO = 0
  if     ( .NOT.LSAME( TRANS, 'N' ).AND. &
           .NOT.LSAME( TRANS, 'T' ).AND. &
           .NOT.LSAME( TRANS, 'C' )      )THEN
     INFO = 1
  ELSE if (  M < 0 )THEN
     INFO = 2
  ELSE if (  N < 0 )THEN
     INFO = 3
  ELSE if (  KL < 0 )THEN
     INFO = 4
  ELSE if (  KU < 0 )THEN
     INFO = 5
  ELSE if (  LDA < ( KL + KU + 1 ) )THEN
     INFO = 8
  ELSE if (  INCX == 0 )THEN
     INFO = 10
  ELSE if (  INCY == 0 )THEN
     INFO = 13
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'SGBMV ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if (  ( M == 0 ).OR.( N == 0 ).OR. &
      ( ( ALPHA == ZERO ).AND.( BETA == ONE ) ) ) &
     return
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
  if (  LSAME( TRANS, 'N' ) )THEN
     LENX = N
     LENY = M
  ELSE
     LENX = M
     LENY = N
  end if
  if (  INCX > 0 )THEN
     KX = 1
  ELSE
     KX = 1 - ( LENX - 1 )*INCX
  end if
  if (  INCY > 0 )THEN
     KY = 1
  ELSE
     KY = 1 - ( LENY - 1 )*INCY
  end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the band part of A.
!
!     First form  y := beta*y.
!
  if (  BETA /= ONE )THEN
     if (  INCY == 1 )THEN
        if (  BETA == ZERO )THEN
           DO 10, I = 1, LENY
              Y( I ) = ZERO
   10          CONTINUE
        ELSE
           DO 20, I = 1, LENY
              Y( I ) = BETA*Y( I )
   20          CONTINUE
        end if
     ELSE
        IY = KY
        if (  BETA == ZERO )THEN
           DO 30, I = 1, LENY
              Y( IY ) = ZERO
              IY      = IY   + INCY
   30          CONTINUE
        ELSE
           DO 40, I = 1, LENY
              Y( IY ) = BETA*Y( IY )
              IY      = IY           + INCY
   40          CONTINUE
        end if
     end if
  end if
  if (  ALPHA == ZERO ) &
     return
  KUP1 = KU + 1
  if (  LSAME( TRANS, 'N' ) )THEN
!
!        Form  y := alpha*A*x + y.
!
     JX = KX
     if (  INCY == 1 )THEN
        DO 60, J = 1, N
           if (  X( JX ) /= ZERO )THEN
              TEMP = ALPHA*X( JX )
              K    = KUP1 - J
              DO 50, I = MAX( 1, J - KU ), MIN( M, J + KL )
                 Y( I ) = Y( I ) + TEMP*A( K + I, J )
   50             CONTINUE
           end if
           JX = JX + INCX
   60       CONTINUE
     ELSE
        DO 80, J = 1, N
           if (  X( JX ) /= ZERO )THEN
              TEMP = ALPHA*X( JX )
              IY   = KY
              K    = KUP1 - J
              DO 70, I = MAX( 1, J - KU ), MIN( M, J + KL )
                 Y( IY ) = Y( IY ) + TEMP*A( K + I, J )
                 IY      = IY      + INCY
   70             CONTINUE
           end if
           JX = JX + INCX
           if (  J > KU ) &
              KY = KY + INCY
   80       CONTINUE
     end if
  ELSE
!
!        Form  y := alpha*A'*x + y.
!
     JY = KY
     if (  INCX == 1 )THEN
        DO 100, J = 1, N
           TEMP = ZERO
           K    = KUP1 - J
           DO 90, I = MAX( 1, J - KU ), MIN( M, J + KL )
              TEMP = TEMP + A( K + I, J )*X( I )
   90          CONTINUE
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
  100       CONTINUE
     ELSE
        DO 120, J = 1, N
           TEMP = ZERO
           IX   = KX
           K    = KUP1 - J
           DO 110, I = MAX( 1, J - KU ), MIN( M, J + KL )
              TEMP = TEMP + A( K + I, J )*X( IX )
              IX   = IX   + INCX
  110          CONTINUE
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
           if (  J > KU ) &
              KX = KX + INCX
  120       CONTINUE
     end if
  end if
!
  return
!
!     End of SGBMV .
!
end
