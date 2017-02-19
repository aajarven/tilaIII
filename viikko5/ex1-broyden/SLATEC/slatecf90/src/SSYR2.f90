subroutine SSYR2 (UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA)
!
!! SSYR2 performs symmetric rank 2 update of a real symmetric matrix.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B4
!***TYPE      SINGLE PRECISION (SSYR2-S, DSYR2-D, CSYR2-C)
!***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J. J., (ANL)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!  SSYR2  performs the symmetric rank 2 operation
!
!     A := alpha*x*y' + alpha*y*x' + A,
!
!  where alpha is a scalar, x and y are n element vectors and A is an n
!  by n symmetric matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - REAL             array of dimension at least
!           ( 1 + ( n - 1)*abs( INCX)).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - REAL             array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
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
!***END PROLOGUE  SSYR2
!     .. Scalar Arguments ..
  REAL               ALPHA
  INTEGER            INCX, INCY, LDA, N
  CHARACTER*1        UPLO
!     .. Array Arguments ..
  REAL               A( LDA, * ), X( * ), Y( * )
!     .. Parameters ..
  REAL               ZERO
  PARAMETER        ( ZERO = 0.0E+0 )
!     .. Local Scalars ..
  REAL               TEMP1, TEMP2
  INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!***FIRST EXECUTABLE STATEMENT  SSYR2
!
!     Test the input parameters.
!
  INFO = 0
  if     ( .NOT.LSAME( UPLO, 'U' ).AND. &
           .NOT.LSAME( UPLO, 'L' )      )THEN
     INFO = 1
  ELSE if (  N < 0 )THEN
     INFO = 2
  ELSE if (  INCX == 0 )THEN
     INFO = 5
  ELSE if (  INCY == 0 )THEN
     INFO = 7
  ELSE if (  LDA < MAX( 1, N ) )THEN
     INFO = 9
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'SSYR2 ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if (  ( N == 0 ).OR.( ALPHA == ZERO ) ) &
     return
!
!     Set up the start points in X and Y if the increments are not both
!     unity.
!
  if (  ( INCX /= 1 ).OR.( INCY /= 1 ) )THEN
     if (  INCX > 0 )THEN
        KX = 1
     ELSE
        KX = 1 - ( N - 1 )*INCX
     end if
     if (  INCY > 0 )THEN
        KY = 1
     ELSE
        KY = 1 - ( N - 1 )*INCY
     end if
     JX = KX
     JY = KY
  end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
  if (  LSAME( UPLO, 'U' ) )THEN
!
!        Form  A  when A is stored in the upper triangle.
!
     if (  ( INCX == 1 ).AND.( INCY == 1 ) )THEN
        DO 20, J = 1, N
           if (  ( X( J ) /= ZERO ).OR.( Y( J ) /= ZERO ) )THEN
              TEMP1 = ALPHA*Y( J )
              TEMP2 = ALPHA*X( J )
              DO 10, I = 1, J
                 A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   10             CONTINUE
           end if
   20       CONTINUE
     ELSE
        DO 40, J = 1, N
           if (  ( X( JX ) /= ZERO ).OR.( Y( JY ) /= ZERO ) )THEN
              TEMP1 = ALPHA*Y( JY )
              TEMP2 = ALPHA*X( JX )
              IX    = KX
              IY    = KY
              DO 30, I = 1, J
                 A( I, J ) = A( I, J ) + X( IX )*TEMP1 &
                                       + Y( IY )*TEMP2
                 IX        = IX        + INCX
                 IY        = IY        + INCY
   30             CONTINUE
           end if
           JX = JX + INCX
           JY = JY + INCY
   40       CONTINUE
     end if
  ELSE
!
!        Form  A  when A is stored in the lower triangle.
!
     if (  ( INCX == 1 ).AND.( INCY == 1 ) )THEN
        DO 60, J = 1, N
           if (  ( X( J ) /= ZERO ).OR.( Y( J ) /= ZERO ) )THEN
              TEMP1 = ALPHA*Y( J )
              TEMP2 = ALPHA*X( J )
              DO 50, I = J, N
                 A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   50             CONTINUE
           end if
   60       CONTINUE
     ELSE
        DO 80, J = 1, N
           if (  ( X( JX ) /= ZERO ).OR.( Y( JY ) /= ZERO ) )THEN
              TEMP1 = ALPHA*Y( JY )
              TEMP2 = ALPHA*X( JX )
              IX    = JX
              IY    = JY
              DO 70, I = J, N
                 A( I, J ) = A( I, J ) + X( IX )*TEMP1 &
                                       + Y( IY )*TEMP2
                 IX        = IX        + INCX
                 IY        = IY        + INCY
   70             CONTINUE
           end if
           JX = JX + INCX
           JY = JY + INCY
   80       CONTINUE
     end if
  end if
!
  return
!
!     End of SSYR2 .
!
end
