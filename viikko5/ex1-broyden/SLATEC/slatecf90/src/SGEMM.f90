subroutine SGEMM (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
     BETA, C, LDC)
!
!! SGEMM multiplies a real general matrix by a real general matrix.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B6
!***TYPE      SINGLE PRECISION (SGEMM-S, DGEMM-D, CGEMM-C)
!***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S. (NAG)
!***DESCRIPTION
!
!  SGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
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
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
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
!***END PROLOGUE  SGEMM
!     .. Scalar Arguments ..
  CHARACTER*1        TRANSA, TRANSB
  INTEGER            M, N, K, LDA, LDB, LDC
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
  LOGICAL            NOTA, NOTB
  INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
  REAL               TEMP
!     .. Parameters ..
  REAL               ONE         , ZERO
  PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!***FIRST EXECUTABLE STATEMENT  SGEMM
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
  NOTA  = LSAME( TRANSA, 'N' )
  NOTB  = LSAME( TRANSB, 'N' )
  if (  NOTA )THEN
     NROWA = M
     NCOLA = K
  ELSE
     NROWA = K
     NCOLA = M
  end if
  if (  NOTB )THEN
     NROWB = K
  ELSE
     NROWB = N
  end if
!
!     Test the input parameters.
!
  INFO = 0
  if (       ( .NOT.NOTA                 ).AND. &
           ( .NOT.LSAME( TRANSA, 'C' ) ).AND. &
           ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
     INFO = 1
  ELSE if (  ( .NOT.NOTB                 ).AND. &
           ( .NOT.LSAME( TRANSB, 'C' ) ).AND. &
           ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
     INFO = 2
  ELSE if (  M   < 0               )THEN
     INFO = 3
  ELSE if (  N   < 0               )THEN
     INFO = 4
  ELSE if (  K   < 0               )THEN
     INFO = 5
  ELSE if (  LDA < MAX( 1, NROWA ) )THEN
     INFO = 8
  ELSE if (  LDB < MAX( 1, NROWB ) )THEN
     INFO = 10
  ELSE if (  LDC < MAX( 1, M     ) )THEN
     INFO = 13
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'SGEMM ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if (  ( M == 0 ).OR.( N == 0 ).OR. &
      ( ( ( ALPHA == ZERO ).OR.( K == 0 ) ).AND.( BETA == ONE ) ) ) &
     return
!
!     And if  alpha.eq.zero.
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
  if (  NOTB )THEN
     if (  NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
        DO 90, J = 1, N
           if (  BETA == ZERO )THEN
              DO 50, I = 1, M
                 C( I, J ) = ZERO
   50             CONTINUE
           ELSE if (  BETA /= ONE )THEN
              DO 60, I = 1, M
                 C( I, J ) = BETA*C( I, J )
   60             CONTINUE
           end if
           DO 80, L = 1, K
              if (  B( L, J ) /= ZERO )THEN
                 TEMP = ALPHA*B( L, J )
                 DO 70, I = 1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
              end if
   80          CONTINUE
   90       CONTINUE
     ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
        DO 120, J = 1, N
           DO 110, I = 1, M
              TEMP = ZERO
              DO 100, L = 1, K
                 TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
              if (  BETA == ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
  110          CONTINUE
  120       CONTINUE
     end if
  ELSE
     if (  NOTA )THEN
!
!           Form  C := alpha*A*B' + beta*C
!
        DO 170, J = 1, N
           if (  BETA == ZERO )THEN
              DO 130, I = 1, M
                 C( I, J ) = ZERO
  130             CONTINUE
           ELSE if (  BETA /= ONE )THEN
              DO 140, I = 1, M
                 C( I, J ) = BETA*C( I, J )
  140             CONTINUE
           end if
           DO 160, L = 1, K
              if (  B( J, L ) /= ZERO )THEN
                 TEMP = ALPHA*B( J, L )
                 DO 150, I = 1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
              end if
  160          CONTINUE
  170       CONTINUE
     ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
        DO 200, J = 1, N
           DO 190, I = 1, M
              TEMP = ZERO
              DO 180, L = 1, K
                 TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
              if (  BETA == ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
  190          CONTINUE
  200       CONTINUE
     end if
  end if
!
  return
!
!     End of SGEMM .
!
end
