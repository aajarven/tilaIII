subroutine CTRMM (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     B, LDB)
!
!! CTRMM multiplies a complex general matrix by a complex triangular matrix.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B6
!***TYPE      COMPLEX (STRMM-S, DTRMM-D, CTRMM-C)
!***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S. (NAG)
!***DESCRIPTION
!
!  CTRMM  performs one of the matrix-matrix operations
!
!     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
!
!  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
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
!              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
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
!  ALPHA  - COMPLEX         .
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, k ), where k is m
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
!  B      - COMPLEX          array of DIMENSION ( LDB, n ).
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
!***END PROLOGUE  CTRMM
!     .. Scalar Arguments ..
  CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
  INTEGER            M, N, LDA, LDB
  COMPLEX            ALPHA
!     .. Array Arguments ..
  COMPLEX            A( LDA, * ), B( LDB, * )
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          CONJG, MAX
!     .. Local Scalars ..
  LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
  INTEGER            I, INFO, J, K, NROWA
  COMPLEX            TEMP
!     .. Parameters ..
  COMPLEX            ONE
  PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
  COMPLEX            ZERO
  PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!***FIRST EXECUTABLE STATEMENT  CTRMM
!
!     Test the input parameters.
!
  LSIDE  = LSAME( SIDE  , 'L' )
  if (  LSIDE )THEN
     NROWA = M
  ELSE
     NROWA = N
  end if
  NOCONJ = LSAME( TRANSA, 'T' )
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
     call XERBLA( 'CTRMM ', INFO )
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
!           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
!
        if (  UPPER )THEN
           DO 120, J = 1, N
              DO 110, I = M, 1, -1
                 TEMP = B( I, J )
                 if (  NOCONJ )THEN
                    if (  NOUNIT ) &
                       TEMP = TEMP*A( I, I )
                    DO 90, K = 1, I - 1
                       TEMP = TEMP + A( K, I )*B( K, J )
   90                   CONTINUE
                 ELSE
                    if (  NOUNIT ) &
                       TEMP = TEMP*CONJG( A( I, I ) )
                    DO 100, K = 1, I - 1
                       TEMP = TEMP + CONJG( A( K, I ) )*B( K, J )
  100                   CONTINUE
                 end if
                 B( I, J ) = ALPHA*TEMP
  110             CONTINUE
  120          CONTINUE
        ELSE
           DO 160, J = 1, N
              DO 150, I = 1, M
                 TEMP = B( I, J )
                 if (  NOCONJ )THEN
                    if (  NOUNIT ) &
                       TEMP = TEMP*A( I, I )
                    DO 130, K = I + 1, M
                       TEMP = TEMP + A( K, I )*B( K, J )
  130                   CONTINUE
                 ELSE
                    if (  NOUNIT ) &
                       TEMP = TEMP*CONJG( A( I, I ) )
                    DO 140, K = I + 1, M
                       TEMP = TEMP + CONJG( A( K, I ) )*B( K, J )
  140                   CONTINUE
                 end if
                 B( I, J ) = ALPHA*TEMP
  150             CONTINUE
  160          CONTINUE
        end if
     end if
  ELSE
     if (  LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*A.
!
        if (  UPPER )THEN
           DO 200, J = N, 1, -1
              TEMP = ALPHA
              if (  NOUNIT ) &
                 TEMP = TEMP*A( J, J )
              DO 170, I = 1, M
                 B( I, J ) = TEMP*B( I, J )
  170             CONTINUE
              DO 190, K = 1, J - 1
                 if (  A( K, J ) /= ZERO )THEN
                    TEMP = ALPHA*A( K, J )
                    DO 180, I = 1, M
                       B( I, J ) = B( I, J ) + TEMP*B( I, K )
  180                   CONTINUE
                 end if
  190             CONTINUE
  200          CONTINUE
        ELSE
           DO 240, J = 1, N
              TEMP = ALPHA
              if (  NOUNIT ) &
                 TEMP = TEMP*A( J, J )
              DO 210, I = 1, M
                 B( I, J ) = TEMP*B( I, J )
  210             CONTINUE
              DO 230, K = J + 1, N
                 if (  A( K, J ) /= ZERO )THEN
                    TEMP = ALPHA*A( K, J )
                    DO 220, I = 1, M
                       B( I, J ) = B( I, J ) + TEMP*B( I, K )
  220                   CONTINUE
                 end if
  230             CONTINUE
  240          CONTINUE
        end if
     ELSE
!
!           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
!
        if (  UPPER )THEN
           DO 280, K = 1, N
              DO 260, J = 1, K - 1
                 if (  A( J, K ) /= ZERO )THEN
                    if (  NOCONJ )THEN
                       TEMP = ALPHA*A( J, K )
                    ELSE
                       TEMP = ALPHA*CONJG( A( J, K ) )
                    end if
                    DO 250, I = 1, M
                       B( I, J ) = B( I, J ) + TEMP*B( I, K )
  250                   CONTINUE
                 end if
  260             CONTINUE
              TEMP = ALPHA
              if (  NOUNIT )THEN
                 if (  NOCONJ )THEN
                    TEMP = TEMP*A( K, K )
                 ELSE
                    TEMP = TEMP*CONJG( A( K, K ) )
                 end if
              end if
              if (  TEMP /= ONE )THEN
                 DO 270, I = 1, M
                    B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
              end if
  280          CONTINUE
        ELSE
           DO 320, K = N, 1, -1
              DO 300, J = K + 1, N
                 if (  A( J, K ) /= ZERO )THEN
                    if (  NOCONJ )THEN
                       TEMP = ALPHA*A( J, K )
                    ELSE
                       TEMP = ALPHA*CONJG( A( J, K ) )
                    end if
                    DO 290, I = 1, M
                       B( I, J ) = B( I, J ) + TEMP*B( I, K )
  290                   CONTINUE
                 end if
  300             CONTINUE
              TEMP = ALPHA
              if (  NOUNIT )THEN
                 if (  NOCONJ )THEN
                    TEMP = TEMP*A( K, K )
                 ELSE
                    TEMP = TEMP*CONJG( A( K, K ) )
                 end if
              end if
              if (  TEMP /= ONE )THEN
                 DO 310, I = 1, M
                    B( I, K ) = TEMP*B( I, K )
  310                CONTINUE
              end if
  320          CONTINUE
        end if
     end if
  end if
!
  return
!
!     End of CTRMM .
!
end
