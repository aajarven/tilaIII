subroutine DTBSV (UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX)
!
!! DTBSV solves a triangular band linear system.
!
!***PURPOSE  Solve one of the systems of equations.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B4
!***TYPE      DOUBLE PRECISION (STBSV-S, DTBSV-D, CTBSV-C)
!***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J. J., (ANL)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!  DTBSV  solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular band matrix, with ( k + 1)
!  diagonals.
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
!  K      - INTEGER.
!           On entry with UPLO = 'U' or 'u', K specifies the number of
!           super-diagonals of the matrix A.
!           On entry with UPLO = 'L' or 'l', K specifies the number of
!           sub-diagonals of the matrix A.
!           K must satisfy  0 .le. K.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer an upper
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = K + 1 - J
!                    DO 10, I = MAX( 1, J - K ), J
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer a lower
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = 1 - J
!                    DO 10, I = J, MIN( N, J + K )
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Note that when DIAG = 'U' or 'u' the elements of the array A
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
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
!***END PROLOGUE  DTBSV
!     .. Scalar Arguments ..
  INTEGER            INCX, K, LDA, N
  CHARACTER*1        DIAG, TRANS, UPLO
!     .. Array Arguments ..
  DOUBLE PRECISION   A( LDA, * ), X( * )
!     .. Parameters ..
  DOUBLE PRECISION   ZERO
  PARAMETER        ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
  DOUBLE PRECISION   TEMP
  INTEGER            I, INFO, IX, J, JX, KPLUS1, KX, L
  LOGICAL            NOUNIT
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!***FIRST EXECUTABLE STATEMENT  DTBSV
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
  ELSE if (  K < 0 )THEN
     INFO = 5
  ELSE if (  LDA < ( K + 1 ) )THEN
     INFO = 7
  ELSE if (  INCX == 0 )THEN
     INFO = 9
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'DTBSV ', INFO )
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
!     accessed by sequentially with one pass through A.
!
  if (  LSAME( TRANS, 'N' ) )THEN
!
!        Form  x := inv( A )*x.
!
     if (  LSAME( UPLO, 'U' ) )THEN
        KPLUS1 = K + 1
        if (  INCX == 1 )THEN
           DO 20, J = N, 1, -1
              if (  X( J ) /= ZERO )THEN
                 L = KPLUS1 - J
                 if (  NOUNIT ) &
                    X( J ) = X( J )/A( KPLUS1, J )
                 TEMP = X( J )
                 DO 10, I = J - 1, MAX( 1, J - K ), -1
                    X( I ) = X( I ) - TEMP*A( L + I, J )
   10                CONTINUE
              end if
   20          CONTINUE
        ELSE
           KX = KX + ( N - 1 )*INCX
           JX = KX
           DO 40, J = N, 1, -1
              KX = KX - INCX
              if (  X( JX ) /= ZERO )THEN
                 IX = KX
                 L  = KPLUS1 - J
                 if (  NOUNIT ) &
                    X( JX ) = X( JX )/A( KPLUS1, J )
                 TEMP = X( JX )
                 DO 30, I = J - 1, MAX( 1, J - K ), -1
                    X( IX ) = X( IX ) - TEMP*A( L + I, J )
                    IX      = IX      - INCX
   30                CONTINUE
              end if
              JX = JX - INCX
   40          CONTINUE
        end if
     ELSE
        if (  INCX == 1 )THEN
           DO 60, J = 1, N
              if (  X( J ) /= ZERO )THEN
                 L = 1 - J
                 if (  NOUNIT ) &
                    X( J ) = X( J )/A( 1, J )
                 TEMP = X( J )
                 DO 50, I = J + 1, MIN( N, J + K )
                    X( I ) = X( I ) - TEMP*A( L + I, J )
   50                CONTINUE
              end if
   60          CONTINUE
        ELSE
           JX = KX
           DO 80, J = 1, N
              KX = KX + INCX
              if (  X( JX ) /= ZERO )THEN
                 IX = KX
                 L  = 1  - J
                 if (  NOUNIT ) &
                    X( JX ) = X( JX )/A( 1, J )
                 TEMP = X( JX )
                 DO 70, I = J + 1, MIN( N, J + K )
                    X( IX ) = X( IX ) - TEMP*A( L + I, J )
                    IX      = IX      + INCX
   70                CONTINUE
              end if
              JX = JX + INCX
   80          CONTINUE
        end if
     end if
  ELSE
!
!        Form  x := inv( A')*x.
!
     if (  LSAME( UPLO, 'U' ) )THEN
        KPLUS1 = K + 1
        if (  INCX == 1 )THEN
           DO 100, J = 1, N
              TEMP = X( J )
              L    = KPLUS1 - J
              DO 90, I = MAX( 1, J - K ), J - 1
                 TEMP = TEMP - A( L + I, J )*X( I )
   90             CONTINUE
              if (  NOUNIT ) &
                 TEMP = TEMP/A( KPLUS1, J )
              X( J ) = TEMP
  100          CONTINUE
        ELSE
           JX = KX
           DO 120, J = 1, N
              TEMP = X( JX )
              IX   = KX
              L    = KPLUS1  - J
              DO 110, I = MAX( 1, J - K ), J - 1
                 TEMP = TEMP - A( L + I, J )*X( IX )
                 IX   = IX   + INCX
  110             CONTINUE
              if (  NOUNIT ) &
                 TEMP = TEMP/A( KPLUS1, J )
              X( JX ) = TEMP
              JX      = JX   + INCX
              if (  J > K ) &
                 KX = KX + INCX
  120          CONTINUE
        end if
     ELSE
        if (  INCX == 1 )THEN
           DO 140, J = N, 1, -1
              TEMP = X( J )
              L    = 1      - J
              DO 130, I = MIN( N, J + K ), J + 1, -1
                 TEMP = TEMP - A( L + I, J )*X( I )
  130             CONTINUE
              if (  NOUNIT ) &
                 TEMP = TEMP/A( 1, J )
              X( J ) = TEMP
  140          CONTINUE
        ELSE
           KX = KX + ( N - 1 )*INCX
           JX = KX
           DO 160, J = N, 1, -1
              TEMP = X( JX )
              IX   = KX
              L    = 1       - J
              DO 150, I = MIN( N, J + K ), J + 1, -1
                 TEMP = TEMP - A( L + I, J )*X( IX )
                 IX   = IX   - INCX
  150             CONTINUE
              if (  NOUNIT ) &
                 TEMP = TEMP/A( 1, J )
              X( JX ) = TEMP
              JX      = JX   - INCX
              if (  ( N - J ) >= K ) &
                 KX = KX - INCX
  160          CONTINUE
        end if
     end if
  end if
!
  return
!
!     End of DTBSV .
!
end
