subroutine CGEMV (TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, &
     INCY)
!
!! CGEMV multiplies a complex vector by a complex general matrix.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B4
!***TYPE      COMPLEX (SGEMV-S, DGEMV-D, CGEMV-C)
!***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J. J., (ANL)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!  CGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
!
!     y := alpha*conjg( A' )*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
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
!              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
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
!  ALPHA  - COMPLEX         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - COMPLEX          array of DIMENSION at least
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
!  BETA   - COMPLEX         .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - COMPLEX          array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
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
!***END PROLOGUE  CGEMV
!     .. Scalar Arguments ..
  COMPLEX            ALPHA, BETA
  INTEGER            INCX, INCY, LDA, M, N
  CHARACTER*1        TRANS
!     .. Array Arguments ..
  COMPLEX            A( LDA, * ), X( * ), Y( * )
!     .. Parameters ..
  COMPLEX            ONE
  PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
  COMPLEX            ZERO
  PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     .. Local Scalars ..
  COMPLEX            TEMP
  INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
  LOGICAL            NOCONJ
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          CONJG, MAX
!***FIRST EXECUTABLE STATEMENT  CGEMV
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
  ELSE if (  LDA < MAX( 1, M ) )THEN
     INFO = 6
  ELSE if (  INCX == 0 )THEN
     INFO = 8
  ELSE if (  INCY == 0 )THEN
     INFO = 11
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'CGEMV ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if (  ( M == 0 ).OR.( N == 0 ).OR. &
      ( ( ALPHA == ZERO ).AND.( BETA == ONE ) ) ) &
     return
!
  NOCONJ = LSAME( TRANS, 'T' )
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
!     accessed sequentially with one pass through A.
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
  if (  LSAME( TRANS, 'N' ) )THEN
!
!        Form  y := alpha*A*x + y.
!
     JX = KX
     if (  INCY == 1 )THEN
        DO 60, J = 1, N
           if (  X( JX ) /= ZERO )THEN
              TEMP = ALPHA*X( JX )
              DO 50, I = 1, M
                 Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
           end if
           JX = JX + INCX
   60       CONTINUE
     ELSE
        DO 80, J = 1, N
           if (  X( JX ) /= ZERO )THEN
              TEMP = ALPHA*X( JX )
              IY   = KY
              DO 70, I = 1, M
                 Y( IY ) = Y( IY ) + TEMP*A( I, J )
                 IY      = IY      + INCY
   70             CONTINUE
           end if
           JX = JX + INCX
   80       CONTINUE
     end if
  ELSE
!
!        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
!
     JY = KY
     if (  INCX == 1 )THEN
        DO 110, J = 1, N
           TEMP = ZERO
           if (  NOCONJ )THEN
              DO 90, I = 1, M
                 TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
           ELSE
              DO 100, I = 1, M
                 TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  100             CONTINUE
           end if
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
  110       CONTINUE
     ELSE
        DO 140, J = 1, N
           TEMP = ZERO
           IX   = KX
           if (  NOCONJ )THEN
              DO 120, I = 1, M
                 TEMP = TEMP + A( I, J )*X( IX )
                 IX   = IX   + INCX
  120             CONTINUE
           ELSE
              DO 130, I = 1, M
                 TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
                 IX   = IX   + INCX
  130             CONTINUE
           end if
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
  140       CONTINUE
     end if
  end if
!
  return
!
!     End of CGEMV .
!
end