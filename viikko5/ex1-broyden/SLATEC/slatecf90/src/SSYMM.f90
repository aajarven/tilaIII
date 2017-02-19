subroutine SSYMM (SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, &
     C, LDC)
!
!! SSYMM multiplies a real general matrix by a real symmetric matrix.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B6
!***TYPE      SINGLE PRECISION (SSYMM-S, DSYMM-D, CSYMM-C)
!***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S. (NAG)
!***DESCRIPTION
!
!  SSYMM  performs one of the matrix-matrix operations
!
!     C := alpha*A*B + beta*C,
!
!  or
!
!     C := alpha*B*A + beta*C,
!
!  where alpha and beta are scalars,  A is a symmetric matrix and  B and
!  C are  m by n matrices.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry,  SIDE  specifies whether  the  symmetric matrix  A
!           appears on the  left or right  in the  operation as follows:
!
!              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
!
!              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!           triangular  part  of  the  symmetric  matrix   A  is  to  be
!           referenced as follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of the
!                                  symmetric matrix is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of the
!                                  symmetric matrix is to be referenced.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies the number of rows of the matrix  C.
!           M  must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix C.
!           N  must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
!           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
!           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
!           the array  A  must contain the  symmetric matrix,  such that
!           when  UPLO = 'U' or 'u', the leading m by m upper triangular
!           part of the array  A  must contain the upper triangular part
!           of the  symmetric matrix and the  strictly  lower triangular
!           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
!           the leading  m by m  lower triangular part  of the  array  A
!           must  contain  the  lower triangular part  of the  symmetric
!           matrix and the  strictly upper triangular part of  A  is not
!           referenced.
!           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
!           the array  A  must contain the  symmetric matrix,  such that
!           when  UPLO = 'U' or 'u', the leading n by n upper triangular
!           part of the array  A  must contain the upper triangular part
!           of the  symmetric matrix and the  strictly  lower triangular
!           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
!           the leading  n by n  lower triangular part  of the  array  A
!           must  contain  the  lower triangular part  of the  symmetric
!           matrix and the  strictly upper triangular part of  A  is not
!           referenced.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  B      - REAL             array of DIMENSION ( LDB, n ).
!           Before entry, the leading  m by n part of the array  B  must
!           contain the matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!  BETA   - REAL            .
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - REAL             array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n updated
!           matrix.
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
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
!***END PROLOGUE  SSYMM
!     .. Scalar Arguments ..
  CHARACTER*1        SIDE, UPLO
  INTEGER            M, N, LDA, LDB, LDC
  REAL               ALPHA, BETA
!     .. Array Arguments ..
  REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            I, INFO, J, K, NROWA
  REAL               TEMP1, TEMP2
!     .. Parameters ..
  REAL               ONE         , ZERO
  PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!***FIRST EXECUTABLE STATEMENT  SSYMM
!
!     Set NROWA as the number of rows of A.
!
  if (  LSAME( SIDE, 'L' ) )THEN
     NROWA = M
  ELSE
     NROWA = N
  end if
  UPPER = LSAME( UPLO, 'U' )
!
!     Test the input parameters.
!
  INFO = 0
  if (       ( .NOT.LSAME( SIDE, 'L' ) ).AND. &
           ( .NOT.LSAME( SIDE, 'R' ) )      )THEN
     INFO = 1
  ELSE if (  ( .NOT.UPPER              ).AND. &
           ( .NOT.LSAME( UPLO, 'L' ) )      )THEN
     INFO = 2
  ELSE if (  M   < 0               )THEN
     INFO = 3
  ELSE if (  N   < 0               )THEN
     INFO = 4
  ELSE if (  LDA < MAX( 1, NROWA ) )THEN
     INFO = 7
  ELSE if (  LDB < MAX( 1, M     ) )THEN
     INFO = 9
  ELSE if (  LDC < MAX( 1, M     ) )THEN
     INFO = 12
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'SSYMM ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if (  ( M == 0 ).OR.( N == 0 ).OR. &
      ( ( ALPHA == ZERO ).AND.( BETA == ONE ) ) ) &
     return
!
!     And when  alpha.eq.zero.
!
  if (  ALPHA == ZERO )THEN
     if (  BETA == ZERO )THEN
        DO 20, J = 1, N
           DO 10, I = 1, M
              C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
     ELSE
        DO 40, J = 1, N
           DO 30, I = 1, M
              C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
     end if
     return
  end if
!
!     Start the operations.
!
  if (  LSAME( SIDE, 'L' ) )THEN
!
!        Form  C := alpha*A*B + beta*C.
!
     if (  UPPER )THEN
        DO 70, J = 1, N
           DO 60, I = 1, M
              TEMP1 = ALPHA*B( I, J )
              TEMP2 = ZERO
              DO 50, K = 1, I - 1
                 C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                 TEMP2     = TEMP2     + B( K, J )*A( K, I )
   50             CONTINUE
              if (  BETA == ZERO )THEN
                 C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
              ELSE
                 C( I, J ) = BETA *C( I, J ) + &
                             TEMP1*A( I, I ) + ALPHA*TEMP2
              end if
   60          CONTINUE
   70       CONTINUE
     ELSE
        DO 100, J = 1, N
           DO 90, I = M, 1, -1
              TEMP1 = ALPHA*B( I, J )
              TEMP2 = ZERO
              DO 80, K = I + 1, M
                 C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                 TEMP2     = TEMP2     + B( K, J )*A( K, I )
   80             CONTINUE
              if (  BETA == ZERO )THEN
                 C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
              ELSE
                 C( I, J ) = BETA *C( I, J ) + &
                             TEMP1*A( I, I ) + ALPHA*TEMP2
              end if
   90          CONTINUE
  100       CONTINUE
     end if
  ELSE
!
!        Form  C := alpha*B*A + beta*C.
!
     DO 170, J = 1, N
        TEMP1 = ALPHA*A( J, J )
        if (  BETA == ZERO )THEN
           DO 110, I = 1, M
              C( I, J ) = TEMP1*B( I, J )
  110          CONTINUE
        ELSE
           DO 120, I = 1, M
              C( I, J ) = BETA*C( I, J ) + TEMP1*B( I, J )
  120          CONTINUE
        end if
        DO 140, K = 1, J - 1
           if (  UPPER )THEN
              TEMP1 = ALPHA*A( K, J )
           ELSE
              TEMP1 = ALPHA*A( J, K )
           end if
           DO 130, I = 1, M
              C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  130          CONTINUE
  140       CONTINUE
        DO 160, K = J + 1, N
           if (  UPPER )THEN
              TEMP1 = ALPHA*A( J, K )
           ELSE
              TEMP1 = ALPHA*A( K, J )
           end if
           DO 150, I = 1, M
              C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE
  end if
!
  return
!
!     End of SSYMM .
!
end
