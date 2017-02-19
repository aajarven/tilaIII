subroutine CHER2 (UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA)
!
!! CHER2 performs Hermitian rank 2 update of a complex Hermitian matrix.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B4
!***TYPE      COMPLEX (SHER2-S, DHER2-D, CHER2-C)
!***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J. J., (ANL)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!  CHER2  performs the hermitian rank 2 operation
!
!     A := alpha*x*conjg( y') + conjg( alpha)*y*conjg( x') + A,
!
!  where alpha is a scalar, x and y are n element vectors and A is an n
!  by n hermitian matrix.
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
!  ALPHA  - COMPLEX         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - COMPLEX          array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - COMPLEX          array of dimension at least
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
!  A      - COMPLEX          array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the hermitian matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the hermitian matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!           Note that the imaginary parts of the diagonal elements need
!           not be set, they are assumed to be zero, and on exit they
!           are set to zero.
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
!***END PROLOGUE  CHER2
!     .. Scalar Arguments ..
  COMPLEX            ALPHA
  INTEGER            INCX, INCY, LDA, N
  CHARACTER*1        UPLO
!     .. Array Arguments ..
  COMPLEX            A( LDA, * ), X( * ), Y( * )
!     .. Parameters ..
  COMPLEX            ZERO
  PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     .. Local Scalars ..
  COMPLEX            TEMP1, TEMP2
  INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          CONJG, MAX, REAL
!***FIRST EXECUTABLE STATEMENT  CHER2
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
     call XERBLA( 'CHER2 ', INFO )
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
              TEMP1 = ALPHA*CONJG( Y( J ) )
              TEMP2 = CONJG( ALPHA*X( J ) )
              DO 10, I = 1, J - 1
                 A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   10             CONTINUE
              A( J, J ) = REAL( A( J, J ) ) + &
                          REAL( X( J )*TEMP1 + Y( J )*TEMP2 )
           ELSE
              A( J, J ) = REAL( A( J, J ) )
           end if
   20       CONTINUE
     ELSE
        DO 40, J = 1, N
           if (  ( X( JX ) /= ZERO ).OR.( Y( JY ) /= ZERO ) )THEN
              TEMP1 = ALPHA*CONJG( Y( JY ) )
              TEMP2 = CONJG( ALPHA*X( JX ) )
              IX    = KX
              IY    = KY
              DO 30, I = 1, J - 1
                 A( I, J ) = A( I, J ) + X( IX )*TEMP1 &
                                       + Y( IY )*TEMP2
                 IX        = IX        + INCX
                 IY        = IY        + INCY
   30             CONTINUE
              A( J, J ) = REAL( A( J, J ) ) + &
                          REAL( X( JX )*TEMP1 + Y( JY )*TEMP2 )
           ELSE
              A( J, J ) = REAL( A( J, J ) )
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
              TEMP1     = ALPHA*CONJG( Y( J ) )
              TEMP2     = CONJG( ALPHA*X( J ) )
              A( J, J ) = REAL( A( J, J ) ) + &
                          REAL( X( J )*TEMP1 + Y( J )*TEMP2 )
              DO 50, I = J + 1, N
                 A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   50             CONTINUE
           ELSE
              A( J, J ) = REAL( A( J, J ) )
           end if
   60       CONTINUE
     ELSE
        DO 80, J = 1, N
           if (  ( X( JX ) /= ZERO ).OR.( Y( JY ) /= ZERO ) )THEN
              TEMP1     = ALPHA*CONJG( Y( JY ) )
              TEMP2     = CONJG( ALPHA*X( JX ) )
              A( J, J ) = REAL( A( J, J ) ) + &
                          REAL( X( JX )*TEMP1 + Y( JY )*TEMP2 )
              IX        = JX
              IY        = JY
              DO 70, I = J + 1, N
                 IX        = IX        + INCX
                 IY        = IY        + INCY
                 A( I, J ) = A( I, J ) + X( IX )*TEMP1 &
                                       + Y( IY )*TEMP2
   70             CONTINUE
           ELSE
              A( J, J ) = REAL( A( J, J ) )
           end if
           JX = JX + INCX
           JY = JY + INCY
   80       CONTINUE
     end if
  end if
!
  return
!
!     End of CHER2 .
!
end
