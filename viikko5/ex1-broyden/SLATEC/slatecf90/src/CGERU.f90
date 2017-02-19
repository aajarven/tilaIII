subroutine CGERU (M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
!
!! CGERU performs unconjugated rank 1 update of a complex general matrix.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B4
!***TYPE      COMPLEX (SGERU-S, DGERU-D, CGERU-C)
!***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J. J., (ANL)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!  CGERU  performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
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
!  X      - COMPLEX          array of dimension at least
!           ( 1 + ( m - 1)*abs( INCX)).
!           Before entry, the incremented array X must contain the m
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
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!***REFERENCES  Dongarra, J. J., Du Croz, J., Hammarling, S., and
!                 Hanson, R. J.  An extended set of Fortran basic linear
!                 algebra subprograms.  ACM TOMS, Vol. 14, No. 1,
!                 pp. 1-17, March 1988.
!***ROUTINES CALLED  XERBLA
!***REVISION HISTORY  (YYMMDD)
!   861022  DATE WRITTEN
!   910605  Modified to meet SLATEC prologue standards.  Only comment
!           lines were modified.  (BKS)
!***END PROLOGUE  CGERU
!     .. Scalar Arguments ..
  COMPLEX            ALPHA
  INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
  COMPLEX            A( LDA, * ), X( * ), Y( * )
!     .. Parameters ..
  COMPLEX            ZERO
  PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     .. Local Scalars ..
  COMPLEX            TEMP
  INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!***FIRST EXECUTABLE STATEMENT  CGERU
!
!     Test the input parameters.
!
  INFO = 0
  if     ( M < 0 )THEN
     INFO = 1
  ELSE if (  N < 0 )THEN
     INFO = 2
  ELSE if (  INCX == 0 )THEN
     INFO = 5
  ELSE if (  INCY == 0 )THEN
     INFO = 7
  ELSE if (  LDA < MAX( 1, M ) )THEN
     INFO = 9
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'CGERU ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if (  ( M == 0 ).OR.( N == 0 ).OR.( ALPHA == ZERO ) ) &
     return
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
  if (  INCY > 0 )THEN
     JY = 1
  ELSE
     JY = 1 - ( N - 1 )*INCY
  end if
  if (  INCX == 1 )THEN
     DO 20, J = 1, N
        if (  Y( JY ) /= ZERO )THEN
           TEMP = ALPHA*Y( JY )
           DO 10, I = 1, M
              A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
        end if
        JY = JY + INCY
   20    CONTINUE
  ELSE
     if (  INCX > 0 )THEN
        KX = 1
     ELSE
        KX = 1 - ( M - 1 )*INCX
     end if
     DO 40, J = 1, N
        if (  Y( JY ) /= ZERO )THEN
           TEMP = ALPHA*Y( JY )
           IX   = KX
           DO 30, I = 1, M
              A( I, J ) = A( I, J ) + X( IX )*TEMP
              IX        = IX        + INCX
   30          CONTINUE
        end if
        JY = JY + INCY
   40    CONTINUE
  end if
!
  return
!
!     End of CGERU .
!
end
