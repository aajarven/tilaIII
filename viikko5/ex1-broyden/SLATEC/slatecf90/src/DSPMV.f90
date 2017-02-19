subroutine DSPMV (UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY)
!
!! DSPMV performs the matrix-vector operation y := alpha*A*x + beta*y.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B4
!***TYPE      DOUBLE PRECISION (SSPMV-S, DSPMV-D, CSPMV-C)
!***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J. J., (ANL)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!  DSPMV  performs the matrix-vector operation
!
!     y := alpha*A*x + beta*y,
!
!  where alpha and beta are scalars, x and y are n element vectors and
!  A is an n by n symmetric matrix, supplied in packed form.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the matrix A is supplied in the packed
!           array AP as follows:
!
!              UPLO = 'U' or 'u'   The upper triangular part of A is
!                                  supplied in AP.
!
!              UPLO = 'L' or 'l'   The lower triangular part of A is
!                                  supplied in AP.
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
!  AP     - DOUBLE PRECISION array of DIMENSION at least
!           ( ( n*( n + 1))/2).
!           Before entry with UPLO = 'U' or 'u', the array AP must
!           contain the upper triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP( 1 )
!           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!           and a( 2, 2 ) respectively, and so on.
!           Before entry with UPLO = 'L' or 'l', the array AP must
!           contain the lower triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP( 1 )
!           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!           and a( 3, 1 ) respectively, and so on.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
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
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y. On exit, Y is overwritten by the updated
!           vector y.
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
!***END PROLOGUE  DSPMV
!     .. Scalar Arguments ..
  DOUBLE PRECISION   ALPHA, BETA
  INTEGER            INCX, INCY, N
  CHARACTER*1        UPLO
!     .. Array Arguments ..
  DOUBLE PRECISION   AP( * ), X( * ), Y( * )
!     .. Parameters ..
  DOUBLE PRECISION   ONE         , ZERO
  PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     .. Local Scalars ..
  DOUBLE PRECISION   TEMP1, TEMP2
  INTEGER            I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!***FIRST EXECUTABLE STATEMENT  DSPMV
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
     INFO = 6
  ELSE if (  INCY == 0 )THEN
     INFO = 9
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'DSPMV ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if (  ( N == 0 ).OR.( ( ALPHA == ZERO ).AND.( BETA == ONE ) ) ) &
     return
!
!     Set up the start points in  X  and  Y.
!
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
!
!     Start the operations. In this version the elements of the array AP
!     are accessed sequentially with one pass through AP.
!
!     First form  y := beta*y.
!
  if (  BETA /= ONE )THEN
     if (  INCY == 1 )THEN
        if (  BETA == ZERO )THEN
           DO 10, I = 1, N
              Y( I ) = ZERO
   10          CONTINUE
        ELSE
           DO 20, I = 1, N
              Y( I ) = BETA*Y( I )
   20          CONTINUE
        end if
     ELSE
        IY = KY
        if (  BETA == ZERO )THEN
           DO 30, I = 1, N
              Y( IY ) = ZERO
              IY      = IY   + INCY
   30          CONTINUE
        ELSE
           DO 40, I = 1, N
              Y( IY ) = BETA*Y( IY )
              IY      = IY           + INCY
   40          CONTINUE
        end if
     end if
  end if
  if (  ALPHA == ZERO ) &
     return
  KK = 1
  if (  LSAME( UPLO, 'U' ) )THEN
!
!        Form  y  when AP contains the upper triangle.
!
     if (  ( INCX == 1 ).AND.( INCY == 1 ) )THEN
        DO 60, J = 1, N
           TEMP1 = ALPHA*X( J )
           TEMP2 = ZERO
           K     = KK
           DO 50, I = 1, J - 1
              Y( I ) = Y( I ) + TEMP1*AP( K )
              TEMP2  = TEMP2  + AP( K )*X( I )
              K      = K      + 1
   50          CONTINUE
           Y( J ) = Y( J ) + TEMP1*AP( KK + J - 1 ) + ALPHA*TEMP2
           KK     = KK     + J
   60       CONTINUE
     ELSE
        JX = KX
        JY = KY
        DO 80, J = 1, N
           TEMP1 = ALPHA*X( JX )
           TEMP2 = ZERO
           IX    = KX
           IY    = KY
           DO 70, K = KK, KK + J - 2
              Y( IY ) = Y( IY ) + TEMP1*AP( K )
              TEMP2   = TEMP2   + AP( K )*X( IX )
              IX      = IX      + INCX
              IY      = IY      + INCY
   70          CONTINUE
           Y( JY ) = Y( JY ) + TEMP1*AP( KK + J - 1 ) + ALPHA*TEMP2
           JX      = JX      + INCX
           JY      = JY      + INCY
           KK      = KK      + J
   80       CONTINUE
     end if
  ELSE
!
!        Form  y  when AP contains the lower triangle.
!
     if (  ( INCX == 1 ).AND.( INCY == 1 ) )THEN
        DO 100, J = 1, N
           TEMP1  = ALPHA*X( J )
           TEMP2  = ZERO
           Y( J ) = Y( J )       + TEMP1*AP( KK )
           K      = KK           + 1
           DO 90, I = J + 1, N
              Y( I ) = Y( I ) + TEMP1*AP( K )
              TEMP2  = TEMP2  + AP( K )*X( I )
              K      = K      + 1
   90          CONTINUE
           Y( J ) = Y( J ) + ALPHA*TEMP2
           KK     = KK     + ( N - J + 1 )
  100       CONTINUE
     ELSE
        JX = KX
        JY = KY
        DO 120, J = 1, N
           TEMP1   = ALPHA*X( JX )
           TEMP2   = ZERO
           Y( JY ) = Y( JY )       + TEMP1*AP( KK )
           IX      = JX
           IY      = JY
           DO 110, K = KK + 1, KK + N - J
              IX      = IX      + INCX
              IY      = IY      + INCY
              Y( IY ) = Y( IY ) + TEMP1*AP( K )
              TEMP2   = TEMP2   + AP( K )*X( IX )
  110          CONTINUE
           Y( JY ) = Y( JY ) + ALPHA*TEMP2
           JX      = JX      + INCX
           JY      = JY      + INCY
           KK      = KK      + ( N - J + 1 )
  120       CONTINUE
     end if
  end if
!
  return
!
!     End of DSPMV .
!
end
