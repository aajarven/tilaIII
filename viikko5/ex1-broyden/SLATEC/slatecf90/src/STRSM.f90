subroutine STRSM (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     B, LDB)
!
!! STRSM solves a triangular system of linear equations with multiple RHS.
!
!***PURPOSE  Solve a real triangular system of equations with multiple
!            right-hand sides.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B6
!***TYPE      SINGLE PRECISION (STRSM-S, DTRSM-D, CTRSM-C)
!***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S. (NAG)
!***DESCRIPTION
!
!  STRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - REAL             array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!***REFERENCES  Dongarra, J., Du Croz, J., Duff, I., and Hammarling, S.
!                 A set of level 3 basic linear algebra subprograms.
!                 ACM TOMS, Vol. 16, No. 1, pp. 1-17, March 1990.
!***ROUTINES CALLED  LSAME, XERBLA
!***REVISION HISTORY  (YYMMDD)
!   890208  DATE WRITTEN
!   910605  Modified to meet SLATEC prologue standards.  Only comment
!           lines were modified.  (BKS)
!***END PROLOGUE  STRSM
!     .. Scalar Arguments ..
  CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
  INTEGER            M, N, LDA, LDB
  REAL               ALPHA
!     .. Array Arguments ..
  REAL               A( LDA, * ), B( LDB, * )
!
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     .. Local Scalars ..
  LOGICAL            LSIDE, NOUNIT, UPPER
  INTEGER            I, INFO, J, K, NROWA
  REAL               TEMP
!     .. Parameters ..
  REAL               ONE         , ZERO
  PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!***FIRST EXECUTABLE STATEMENT  STRSM
!
!     Test the input parameters.
!
  LSIDE  = LSAME( SIDE  , 'L' )
  if (  LSIDE )THEN
     NROWA = M
  ELSE
     NROWA = N
  end if
  NOUNIT = LSAME( DIAG  , 'N' )
  UPPER  = LSAME( UPLO  , 'U' )
!
  INFO   = 0
  if (       ( .NOT.LSIDE                ).AND. &
           ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
     INFO = 1
  ELSE if (  ( .NOT.UPPER                ).AND. &
           ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
     INFO = 2
  ELSE if (  ( .NOT.LSAME( TRANSA, 'N' ) ).AND. &
           ( .NOT.LSAME( TRANSA, 'T' ) ).AND. &
           ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
     INFO = 3
  ELSE if (  ( .NOT.LSAME( DIAG  , 'U' ) ).AND. &
           ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
     INFO = 4
  ELSE if (  M   < 0               )THEN
     INFO = 5
  ELSE if (  N   < 0               )THEN
     INFO = 6
  ELSE if (  LDA < MAX( 1, NROWA ) )THEN
     INFO = 9
  ELSE if (  LDB < MAX( 1, M     ) )THEN
     INFO = 11
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'STRSM ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if (  N == 0 ) &
     return
!
!     And when  alpha.eq.zero.
!
  if (  ALPHA == ZERO )THEN
     DO 20, J = 1, N
        DO 10, I = 1, M
           B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
     return
  end if
!
!     Start the operations.
!
  if (  LSIDE )THEN
     if (  LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
        if (  UPPER )THEN
           DO 60, J = 1, N
              if (  ALPHA /= ONE )THEN
                 DO 30, I = 1, M
                    B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
              end if
              DO 50, K = M, 1, -1
                 if (  B( K, J ) /= ZERO )THEN
                    if (  NOUNIT ) &
                       B( K, J ) = B( K, J )/A( K, K )
                    DO 40, I = 1, K - 1
                       B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                 end if
   50             CONTINUE
   60          CONTINUE
        ELSE
           DO 100, J = 1, N
              if (  ALPHA /= ONE )THEN
                 DO 70, I = 1, M
                    B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
              end if
              DO 90 K = 1, M
                 if (  B( K, J ) /= ZERO )THEN
                    if (  NOUNIT ) &
                       B( K, J ) = B( K, J )/A( K, K )
                    DO 80, I = K + 1, M
                       B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                 end if
   90             CONTINUE
  100          CONTINUE
        end if
     ELSE
!
!           Form  B := alpha*inv( A' )*B.
!
        if (  UPPER )THEN
           DO 130, J = 1, N
              DO 120, I = 1, M
                 TEMP = ALPHA*B( I, J )
                 DO 110, K = 1, I - 1
                    TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                 if (  NOUNIT ) &
                    TEMP = TEMP/A( I, I )
                 B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
        ELSE
           DO 160, J = 1, N
              DO 150, I = M, 1, -1
                 TEMP = ALPHA*B( I, J )
                 DO 140, K = I + 1, M
                    TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                 if (  NOUNIT ) &
                    TEMP = TEMP/A( I, I )
                 B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
        end if
     end if
  ELSE
     if (  LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
        if (  UPPER )THEN
           DO 210, J = 1, N
              if (  ALPHA /= ONE )THEN
                 DO 170, I = 1, M
                    B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
              end if
              DO 190, K = 1, J - 1
                 if (  A( K, J ) /= ZERO )THEN
                    DO 180, I = 1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                 end if
  190             CONTINUE
              if (  NOUNIT )THEN
                 TEMP = ONE/A( J, J )
                 DO 200, I = 1, M
                    B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
              end if
  210          CONTINUE
        ELSE
           DO 260, J = N, 1, -1
              if (  ALPHA /= ONE )THEN
                 DO 220, I = 1, M
                    B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
              end if
              DO 240, K = J + 1, N
                 if (  A( K, J ) /= ZERO )THEN
                    DO 230, I = 1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                 end if
  240             CONTINUE
              if (  NOUNIT )THEN
                 TEMP = ONE/A( J, J )
                 DO 250, I = 1, M
                   B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
              end if
  260          CONTINUE
        end if
     ELSE
!
!           Form  B := alpha*B*inv( A' ).
!
        if (  UPPER )THEN
           DO 310, K = N, 1, -1
              if (  NOUNIT )THEN
                 TEMP = ONE/A( K, K )
                 DO 270, I = 1, M
                    B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
              end if
              DO 290, J = 1, K - 1
                 if (  A( J, K ) /= ZERO )THEN
                    TEMP = A( J, K )
                    DO 280, I = 1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                 end if
  290             CONTINUE
              if (  ALPHA /= ONE )THEN
                 DO 300, I = 1, M
                    B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
              end if
  310          CONTINUE
        ELSE
           DO 360, K = 1, N
              if (  NOUNIT )THEN
                 TEMP = ONE/A( K, K )
                 DO 320, I = 1, M
                    B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
              end if
              DO 340, J = K + 1, N
                 if (  A( J, K ) /= ZERO )THEN
                    TEMP = A( J, K )
                    DO 330, I = 1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                 end if
  340             CONTINUE
              if (  ALPHA /= ONE )THEN
                 DO 350, I = 1, M
                    B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
              end if
  360          CONTINUE
        end if
     end if
  end if
!
  return
!
!     End of STRSM .
!
end
