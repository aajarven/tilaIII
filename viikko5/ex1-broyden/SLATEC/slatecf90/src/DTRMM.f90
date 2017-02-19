subroutine DTRMM (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     B, LDB)
!
!! DTRMM performs B = alpha*op(A)*B or B = alpha*B*op(A), A triangular.
!
!***PURPOSE  Perform one of the matrix-matrix operations.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B6
!***TYPE      DOUBLE PRECISION (STRMM-S, DTRMM-D, CTRMM-C)
!***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S. (NAG)
!***DESCRIPTION
!
!  DTRMM  performs one of the matrix-matrix operations
!
!     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!
!  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry,  SIDE specifies whether  op( A ) multiplies B from
!           the left or right as follows:
!
!              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!
!              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
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
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
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
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain the matrix  B,  and  on exit  is overwritten  by the
!           transformed matrix.
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
!***END PROLOGUE  DTRMM
!     .. Scalar Arguments ..
  CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
  INTEGER            M, N, LDA, LDB
  DOUBLE PRECISION   ALPHA
!     .. Array Arguments ..
  DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
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
  DOUBLE PRECISION   TEMP
!     .. Parameters ..
  DOUBLE PRECISION   ONE         , ZERO
  PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!***FIRST EXECUTABLE STATEMENT  DTRMM
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
     call XERBLA( 'DTRMM ', INFO )
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
!           Form  B := alpha*A*B.
!
        if (  UPPER )THEN
           DO 50, J = 1, N
              DO 40, K = 1, M
                 if (  B( K, J ) /= ZERO )THEN
                    TEMP = ALPHA*B( K, J )
                    DO 30, I = 1, K - 1
                       B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                    if (  NOUNIT ) &
                       TEMP = TEMP*A( K, K )
                    B( K, J ) = TEMP
                 end if
   40             CONTINUE
   50          CONTINUE
        ELSE
           DO 80, J = 1, N
              DO 70 K = M, 1, -1
                 if (  B( K, J ) /= ZERO )THEN
                    TEMP      = ALPHA*B( K, J )
                    B( K, J ) = TEMP
                    if (  NOUNIT ) &
                       B( K, J ) = B( K, J )*A( K, K )
                    DO 60, I = K + 1, M
                       B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                 end if
   70             CONTINUE
   80          CONTINUE
        end if
     ELSE
!
!           Form  B := alpha*B*A'.
!
        if (  UPPER )THEN
           DO 110, J = 1, N
              DO 100, I = M, 1, -1
                 TEMP = B( I, J )
                 if (  NOUNIT ) &
                    TEMP = TEMP*A( I, I )
                 DO 90, K = 1, I - 1
                    TEMP = TEMP + A( K, I )*B( K, J )
   90                CONTINUE
                 B( I, J ) = ALPHA*TEMP
  100             CONTINUE
  110          CONTINUE
        ELSE
           DO 140, J = 1, N
              DO 130, I = 1, M
                 TEMP = B( I, J )
                 if (  NOUNIT ) &
                    TEMP = TEMP*A( I, I )
                 DO 120, K = I + 1, M
                    TEMP = TEMP + A( K, I )*B( K, J )
  120                CONTINUE
                 B( I, J ) = ALPHA*TEMP
  130             CONTINUE
  140          CONTINUE
        end if
     end if
  ELSE
     if (  LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*A.
!
        if (  UPPER )THEN
           DO 180, J = N, 1, -1
              TEMP = ALPHA
              if (  NOUNIT ) &
                 TEMP = TEMP*A( J, J )
              DO 150, I = 1, M
                 B( I, J ) = TEMP*B( I, J )
  150             CONTINUE
              DO 170, K = 1, J - 1
                 if (  A( K, J ) /= ZERO )THEN
                    TEMP = ALPHA*A( K, J )
                    DO 160, I = 1, M
                       B( I, J ) = B( I, J ) + TEMP*B( I, K )
  160                   CONTINUE
                 end if
  170             CONTINUE
  180          CONTINUE
        ELSE
           DO 220, J = 1, N
              TEMP = ALPHA
              if (  NOUNIT ) &
                 TEMP = TEMP*A( J, J )
              DO 190, I = 1, M
                 B( I, J ) = TEMP*B( I, J )
  190             CONTINUE
              DO 210, K = J + 1, N
                 if (  A( K, J ) /= ZERO )THEN
                    TEMP = ALPHA*A( K, J )
                    DO 200, I = 1, M
                       B( I, J ) = B( I, J ) + TEMP*B( I, K )
  200                   CONTINUE
                 end if
  210             CONTINUE
  220          CONTINUE
        end if
     ELSE
!
!           Form  B := alpha*B*A'.
!
        if (  UPPER )THEN
           DO 260, K = 1, N
              DO 240, J = 1, K - 1
                 if (  A( J, K ) /= ZERO )THEN
                    TEMP = ALPHA*A( J, K )
                    DO 230, I = 1, M
                       B( I, J ) = B( I, J ) + TEMP*B( I, K )
  230                   CONTINUE
                 end if
  240             CONTINUE
              TEMP = ALPHA
              if (  NOUNIT ) &
                 TEMP = TEMP*A( K, K )
              if (  TEMP /= ONE )THEN
                 DO 250, I = 1, M
                    B( I, K ) = TEMP*B( I, K )
  250                CONTINUE
              end if
  260          CONTINUE
        ELSE
           DO 300, K = N, 1, -1
              DO 280, J = K + 1, N
                 if (  A( J, K ) /= ZERO )THEN
                    TEMP = ALPHA*A( J, K )
                    DO 270, I = 1, M
                       B( I, J ) = B( I, J ) + TEMP*B( I, K )
  270                   CONTINUE
                 end if
  280             CONTINUE
              TEMP = ALPHA
              if (  NOUNIT ) &
                 TEMP = TEMP*A( K, K )
              if (  TEMP /= ONE )THEN
                 DO 290, I = 1, M
                    B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
              end if
  300          CONTINUE
        end if
     end if
  end if
!
  return
!
!     End of DTRMM .
!
end
