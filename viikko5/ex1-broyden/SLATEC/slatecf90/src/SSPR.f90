subroutine SSPR (UPLO, N, ALPHA, X, INCX, AP)
!
!! SSPR performs the symmetric rank 1 operation A = A + alpha*x*x'.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B4
!***TYPE      SINGLE PRECISION (SSPR-S)
!***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J. J., (ANL)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!  SSPR    performs the symmetric rank 1 operation
!
!     A := alpha*x*x' + A,
!
!  where alpha is a real scalar, x is an n element vector and A is an
!  n by n symmetric matrix, supplied in packed form.
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
!  AP     - REAL             array of DIMENSION at least
!           ( ( n*( n + 1 ) )/2 ).
!           Before entry with  UPLO = 'U' or 'u', the array AP must
!           contain the upper triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP( 1 )
!           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!           and a( 2, 2 ) respectively, and so on. On exit, the array
!           AP is overwritten by the upper triangular part of the
!           updated matrix.
!           Before entry with UPLO = 'L' or 'l', the array AP must
!           contain the lower triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP( 1 )
!           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!           and a( 3, 1 ) respectively, and so on. On exit, the array
!           AP is overwritten by the lower triangular part of the
!           updated matrix.
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
!***END PROLOGUE  SSPR
!     .. Scalar Arguments ..
  REAL               ALPHA
  INTEGER            INCX, N
  CHARACTER*1        UPLO
!     .. Array Arguments ..
  REAL               AP( * ), X( * )
!     .. Parameters ..
  REAL               ZERO
  PARAMETER        ( ZERO = 0.0E+0 )
!     .. Local Scalars ..
  REAL               TEMP
  INTEGER            I, INFO, IX, J, JX, K, KK, KX
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!***FIRST EXECUTABLE STATEMENT  SSPR
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
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'SSPR  ', INFO )
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
!     Start the operations. In this version the elements of the array AP
!     are accessed sequentially with one pass through AP.
!
  KK = 1
  if (  LSAME( UPLO, 'U' ) )THEN
!
!        Form  A  when upper triangle is stored in AP.
!
     if (  INCX == 1 )THEN
        DO 20, J = 1, N
           if (  X( J ) /= ZERO )THEN
              TEMP = ALPHA*X( J )
              K    = KK
              DO 10, I = 1, J
                 AP( K ) = AP( K ) + X( I )*TEMP
                 K       = K       + 1
   10             CONTINUE
           end if
           KK = KK + J
   20       CONTINUE
     ELSE
        JX = KX
        DO 40, J = 1, N
           if (  X( JX ) /= ZERO )THEN
              TEMP = ALPHA*X( JX )
              IX   = KX
              DO 30, K = KK, KK + J - 1
                 AP( K ) = AP( K ) + X( IX )*TEMP
                 IX      = IX      + INCX
   30             CONTINUE
           end if
           JX = JX + INCX
           KK = KK + J
   40       CONTINUE
     end if
  ELSE
!
!        Form  A  when lower triangle is stored in AP.
!
     if (  INCX == 1 )THEN
        DO 60, J = 1, N
           if (  X( J ) /= ZERO )THEN
              TEMP = ALPHA*X( J )
              K    = KK
              DO 50, I = J, N
                 AP( K ) = AP( K ) + X( I )*TEMP
                 K       = K       + 1
   50             CONTINUE
           end if
           KK = KK + N - J + 1
   60       CONTINUE
     ELSE
        JX = KX
        DO 80, J = 1, N
           if (  X( JX ) /= ZERO )THEN
              TEMP = ALPHA*X( JX )
              IX   = JX
              DO 70, K = KK, KK + N - J
                 AP( K ) = AP( K ) + X( IX )*TEMP
                 IX      = IX      + INCX
   70             CONTINUE
           end if
           JX = JX + INCX
           KK = KK + N - J + 1
   80       CONTINUE
     end if
  end if
!
  return
!
!     End of SSPR  .
!
end
