subroutine DTRSV (UPLO, TRANS, DIAG, N, A, LDA, X, INCX)
!
!! DTRSV solves A*x=b or A'*x=b where A is triangular.
!
!***PURPOSE  Solve one of the systems of equations.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B4
!***TYPE      DOUBLE PRECISION (STRSV-S, DTRSV-D, CTRSV-C)
!***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J. J., (ANL)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!  DTRSV  solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular matrix.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   A*x = b.
!
!              TRANS = 'T' or 't'   A'*x = b.
!
!              TRANS = 'C' or 'c'   A'*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
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
!***END PROLOGUE  DTRSV
!     .. Scalar Arguments ..
  INTEGER            INCX, LDA, N
  CHARACTER*1        DIAG, TRANS, UPLO
!     .. Array Arguments ..
  DOUBLE PRECISION   A( LDA, * ), X( * )
!     .. Parameters ..
  DOUBLE PRECISION   ZERO
  PARAMETER        ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
  DOUBLE PRECISION   TEMP
  INTEGER            I, INFO, IX, J, JX, KX
  LOGICAL            NOUNIT
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!***FIRST EXECUTABLE STATEMENT  DTRSV
!
!     Test the input parameters.
!
  INFO = 0
  if     ( .NOT.LSAME( UPLO , 'U' ).AND. &
           .NOT.LSAME( UPLO , 'L' )      )THEN
     INFO = 1
  ELSE if (  .NOT.LSAME( TRANS, 'N' ).AND. &
           .NOT.LSAME( TRANS, 'T' ).AND. &
           .NOT.LSAME( TRANS, 'C' )      )THEN
     INFO = 2
  ELSE if (  .NOT.LSAME( DIAG , 'U' ).AND. &
           .NOT.LSAME( DIAG , 'N' )      )THEN
     INFO = 3
  ELSE if (  N < 0 )THEN
     INFO = 4
  ELSE if (  LDA < MAX( 1, N ) )THEN
     INFO = 6
  ELSE if (  INCX == 0 )THEN
     INFO = 8
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'DTRSV ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if (  N == 0 ) &
     return
!
  NOUNIT = LSAME( DIAG, 'N' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
  if (  INCX <= 0 )THEN
     KX = 1 - ( N - 1 )*INCX
  ELSE if (  INCX /= 1 )THEN
     KX = 1
  end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
  if (  LSAME( TRANS, 'N' ) )THEN
!
!        Form  x := inv( A )*x.
!
     if (  LSAME( UPLO, 'U' ) )THEN
        if (  INCX == 1 )THEN
           DO 20, J = N, 1, -1
              if (  X( J ) /= ZERO )THEN
                 if (  NOUNIT ) &
                    X( J ) = X( J )/A( J, J )
                 TEMP = X( J )
                 DO 10, I = J - 1, 1, -1
                    X( I ) = X( I ) - TEMP*A( I, J )
   10                CONTINUE
              end if
   20          CONTINUE
        ELSE
           JX = KX + ( N - 1 )*INCX
           DO 40, J = N, 1, -1
              if (  X( JX ) /= ZERO )THEN
                 if (  NOUNIT ) &
                    X( JX ) = X( JX )/A( J, J )
                 TEMP = X( JX )
                 IX   = JX
                 DO 30, I = J - 1, 1, -1
                    IX      = IX      - INCX
                    X( IX ) = X( IX ) - TEMP*A( I, J )
   30                CONTINUE
              end if
              JX = JX - INCX
   40          CONTINUE
        end if
     ELSE
        if (  INCX == 1 )THEN
           DO 60, J = 1, N
              if (  X( J ) /= ZERO )THEN
                 if (  NOUNIT ) &
                    X( J ) = X( J )/A( J, J )
                 TEMP = X( J )
                 DO 50, I = J + 1, N
                    X( I ) = X( I ) - TEMP*A( I, J )
   50                CONTINUE
              end if
   60          CONTINUE
        ELSE
           JX = KX
           DO 80, J = 1, N
              if (  X( JX ) /= ZERO )THEN
                 if (  NOUNIT ) &
                    X( JX ) = X( JX )/A( J, J )
                 TEMP = X( JX )
                 IX   = JX
                 DO 70, I = J + 1, N
                    IX      = IX      + INCX
                    X( IX ) = X( IX ) - TEMP*A( I, J )
   70                CONTINUE
              end if
              JX = JX + INCX
   80          CONTINUE
        end if
     end if
  ELSE
!
!        Form  x := inv( A' )*x.
!
     if (  LSAME( UPLO, 'U' ) )THEN
        if (  INCX == 1 )THEN
           DO 100, J = 1, N
              TEMP = X( J )
              DO 90, I = 1, J - 1
                 TEMP = TEMP - A( I, J )*X( I )
   90             CONTINUE
              if (  NOUNIT ) &
                 TEMP = TEMP/A( J, J )
              X( J ) = TEMP
  100          CONTINUE
        ELSE
           JX = KX
           DO 120, J = 1, N
              TEMP = X( JX )
              IX   = KX
              DO 110, I = 1, J - 1
                 TEMP = TEMP - A( I, J )*X( IX )
                 IX   = IX   + INCX
  110             CONTINUE
              if (  NOUNIT ) &
                 TEMP = TEMP/A( J, J )
              X( JX ) = TEMP
              JX      = JX   + INCX
  120          CONTINUE
        end if
     ELSE
        if (  INCX == 1 )THEN
           DO 140, J = N, 1, -1
              TEMP = X( J )
              DO 130, I = N, J + 1, -1
                 TEMP = TEMP - A( I, J )*X( I )
  130             CONTINUE
              if (  NOUNIT ) &
                 TEMP = TEMP/A( J, J )
              X( J ) = TEMP
  140          CONTINUE
        ELSE
           KX = KX + ( N - 1 )*INCX
           JX = KX
           DO 160, J = N, 1, -1
              TEMP = X( JX )
              IX   = KX
              DO 150, I = N, J + 1, -1
                 TEMP = TEMP - A( I, J )*X( IX )
                 IX   = IX   - INCX
  150             CONTINUE
              if (  NOUNIT ) &
                 TEMP = TEMP/A( J, J )
              X( JX ) = TEMP
              JX      = JX   - INCX
  160          CONTINUE
        end if
     end if
  end if

  return
end