subroutine DSYR (UPLO, N, ALPHA, X, INCX, A, LDA)
!
!! DSYR performs A = alpha*x*x' + A.
!
!***PURPOSE  Perform the symmetric rank 1 operation.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B4
!***TYPE      DOUBLE PRECISION (DSYR-D)
!***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J. J., (ANL)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!  DSYR   performs the symmetric rank 1 operation
!
!     A := alpha*x*x' + A,
!
!  where alpha is a real scalar, x is an n element vector and A is an
!  n by n symmetric matrix.
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
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
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
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
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
!***END PROLOGUE  DSYR
!     .. Scalar Arguments ..
  DOUBLE PRECISION   ALPHA
  INTEGER            INCX, LDA, N
  CHARACTER*1        UPLO
!     .. Array Arguments ..
  DOUBLE PRECISION   A( LDA, * ), X( * )
!     .. Parameters ..
  DOUBLE PRECISION   ZERO
  PARAMETER        ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
  DOUBLE PRECISION   TEMP
  INTEGER            I, INFO, IX, J, JX, KX
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!***FIRST EXECUTABLE STATEMENT  DSYR
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
  ELSE if (  LDA < MAX( 1, N ) )THEN
     INFO = 7
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'DSYR  ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if (  ( N == 0 ).OR.( ALPHA == ZERO ) ) &
     return
!
!     Set the start point in X if the increment is not unity.
!
  if (  INCX <= 0 )THEN
     KX = 1 - ( N - 1 )*INCX
  ELSE if (  INCX /= 1 )THEN
     KX = 1
  end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
  if (  LSAME( UPLO, 'U' ) )THEN
!
!        Form  A  when A is stored in upper triangle.
!
     if (  INCX == 1 )THEN
        DO 20, J = 1, N
           if (  X( J ) /= ZERO )THEN
              TEMP = ALPHA*X( J )
              DO 10, I = 1, J
                 A( I, J ) = A( I, J ) + X( I )*TEMP
   10             CONTINUE
           end if
   20       CONTINUE
     ELSE
        JX = KX
        DO 40, J = 1, N
           if (  X( JX ) /= ZERO )THEN
              TEMP = ALPHA*X( JX )
              IX   = KX
              DO 30, I = 1, J
                 A( I, J ) = A( I, J ) + X( IX )*TEMP
                 IX        = IX        + INCX
   30             CONTINUE
           end if
           JX = JX + INCX
   40       CONTINUE
     end if
  ELSE
!
!        Form  A  when A is stored in lower triangle.
!
     if (  INCX == 1 )THEN
        DO 60, J = 1, N
           if (  X( J ) /= ZERO )THEN
              TEMP = ALPHA*X( J )
              DO 50, I = J, N
                 A( I, J ) = A( I, J ) + X( I )*TEMP
   50             CONTINUE
           end if
   60       CONTINUE
     ELSE
        JX = KX
        DO 80, J = 1, N
           if (  X( JX ) /= ZERO )THEN
              TEMP = ALPHA*X( JX )
              IX   = JX
              DO 70, I = J, N
                 A( I, J ) = A( I, J ) + X( IX )*TEMP
                 IX        = IX        + INCX
   70             CONTINUE
           end if
           JX = JX + INCX
   80       CONTINUE
     end if
  end if
!
  return
end