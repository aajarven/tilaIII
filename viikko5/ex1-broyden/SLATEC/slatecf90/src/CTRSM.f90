subroutine CTRSM (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     B, LDB)
!
!! CTRSM solves a complex triangular system of equations with ...
!            multiple right-hand sides.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B6
!***TYPE      COMPLEX (STRSM-S, DTRSM-D, CTRSM-C)
!***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S. (NAG)
!***DESCRIPTION
!
!  CTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
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
!***END PROLOGUE  CTRSM
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
!***FIRST EXECUTABLE STATEMENT  CTRSM
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
     call XERBLA( 'CTRSM ', INFO )
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
!           Form  B := alpha*inv( A' )*B
!           or    B := alpha*inv( conjg( A' ) )*B.
!
        if (  UPPER )THEN
           DO 140, J = 1, N
              DO 130, I = 1, M
                 TEMP = ALPHA*B( I, J )
                 if (  NOCONJ )THEN
                    DO 110, K = 1, I - 1
                       TEMP = TEMP - A( K, I )*B( K, J )
  110                   CONTINUE
                    if (  NOUNIT ) &
                       TEMP = TEMP/A( I, I )
                 ELSE
                    DO 120, K = 1, I - 1
                       TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  120                   CONTINUE
                    if (  NOUNIT ) &
                       TEMP = TEMP/CONJG( A( I, I ) )
                 end if
                 B( I, J ) = TEMP
  130             CONTINUE
  140          CONTINUE
        ELSE
           DO 180, J = 1, N
              DO 170, I = M, 1, -1
                 TEMP = ALPHA*B( I, J )
                 if (  NOCONJ )THEN
                    DO 150, K = I + 1, M
                       TEMP = TEMP - A( K, I )*B( K, J )
  150                   CONTINUE
                    if (  NOUNIT ) &
                       TEMP = TEMP/A( I, I )
                 ELSE
                    DO 160, K = I + 1, M
                       TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  160                   CONTINUE
                    if (  NOUNIT ) &
                       TEMP = TEMP/CONJG( A( I, I ) )
                 end if
                 B( I, J ) = TEMP
  170             CONTINUE
  180          CONTINUE
        end if
     end if
  ELSE
     if (  LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
        if (  UPPER )THEN
           DO 230, J = 1, N
              if (  ALPHA /= ONE )THEN
                 DO 190, I = 1, M
                    B( I, J ) = ALPHA*B( I, J )
  190                CONTINUE
              end if
              DO 210, K = 1, J - 1
                 if (  A( K, J ) /= ZERO )THEN
                    DO 200, I = 1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  200                   CONTINUE
                 end if
  210             CONTINUE
              if (  NOUNIT )THEN
                 TEMP = ONE/A( J, J )
                 DO 220, I = 1, M
                    B( I, J ) = TEMP*B( I, J )
  220                CONTINUE
              end if
  230          CONTINUE
        ELSE
           DO 280, J = N, 1, -1
              if (  ALPHA /= ONE )THEN
                 DO 240, I = 1, M
                    B( I, J ) = ALPHA*B( I, J )
  240                CONTINUE
              end if
              DO 260, K = J + 1, N
                 if (  A( K, J ) /= ZERO )THEN
                    DO 250, I = 1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  250                   CONTINUE
                 end if
  260             CONTINUE
              if (  NOUNIT )THEN
                 TEMP = ONE/A( J, J )
                 DO 270, I = 1, M
                   B( I, J ) = TEMP*B( I, J )
  270                CONTINUE
              end if
  280          CONTINUE
        end if
     ELSE
!
!           Form  B := alpha*B*inv( A' )
!           or    B := alpha*B*inv( conjg( A' ) ).
!
        if (  UPPER )THEN
           DO 330, K = N, 1, -1
              if (  NOUNIT )THEN
                 if (  NOCONJ )THEN
                    TEMP = ONE/A( K, K )
                 ELSE
                    TEMP = ONE/CONJG( A( K, K ) )
                 end if
                 DO 290, I = 1, M
                    B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
              end if
              DO 310, J = 1, K - 1
                 if (  A( J, K ) /= ZERO )THEN
                    if (  NOCONJ )THEN
                       TEMP = A( J, K )
                    ELSE
                       TEMP = CONJG( A( J, K ) )
                    end if
                    DO 300, I = 1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
  300                   CONTINUE
                 end if
  310             CONTINUE
              if (  ALPHA /= ONE )THEN
                 DO 320, I = 1, M
                    B( I, K ) = ALPHA*B( I, K )
  320                CONTINUE
              end if
  330          CONTINUE
        ELSE
           DO 380, K = 1, N
              if (  NOUNIT )THEN
                 if (  NOCONJ )THEN
                    TEMP = ONE/A( K, K )
                 ELSE
                    TEMP = ONE/CONJG( A( K, K ) )
                 end if
                 DO 340, I = 1, M
                    B( I, K ) = TEMP*B( I, K )
  340                CONTINUE
              end if
              DO 360, J = K + 1, N
                 if (  A( J, K ) /= ZERO )THEN
                    if (  NOCONJ )THEN
                       TEMP = A( J, K )
                    ELSE
                       TEMP = CONJG( A( J, K ) )
                    end if
                    DO 350, I = 1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
  350                   CONTINUE
                 end if
  360             CONTINUE
              if (  ALPHA /= ONE )THEN
                 DO 370, I = 1, M
                    B( I, K ) = ALPHA*B( I, K )
  370                CONTINUE
              end if
  380          CONTINUE
        end if
     end if
  end if
!
  return
!
!     End of CTRSM .
!
end
