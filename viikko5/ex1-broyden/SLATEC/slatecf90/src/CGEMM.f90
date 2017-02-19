subroutine CGEMM (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
     BETA, C, LDC)
!
!! CGEMM multiplies a complex general matrix by a complex general matrix.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B6
!***TYPE      COMPLEX (SGEMM-S, DGEMM-D, CGEMM-C)
!***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S. (NAG)
!***DESCRIPTION
!
!  CGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
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
!              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
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
!              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
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
!  ALPHA  - COMPLEX         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
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
!  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
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
!  BETA   - COMPLEX         .
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - COMPLEX          array of DIMENSION ( LDC, n ).
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
!***END PROLOGUE  CGEMM
!     .. Scalar Arguments ..
  CHARACTER*1        TRANSA, TRANSB
  INTEGER            M, N, K, LDA, LDB, LDC
  COMPLEX            ALPHA, BETA
!     .. Array Arguments ..
  COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          CONJG, MAX
!     .. Local Scalars ..
  LOGICAL            CONJA, CONJB, NOTA, NOTB
  INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
  COMPLEX            TEMP
!     .. Parameters ..
  COMPLEX            ONE
  PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
  COMPLEX            ZERO
  PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!***FIRST EXECUTABLE STATEMENT  CGEMM
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
!     B  respectively are to be  transposed but  not conjugated  and set
!     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
!     and the number of rows of  B  respectively.
!
  NOTA  = LSAME( TRANSA, 'N' )
  NOTB  = LSAME( TRANSB, 'N' )
  CONJA = LSAME( TRANSA, 'C' )
  CONJB = LSAME( TRANSB, 'C' )
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
           ( .NOT.CONJA                ).AND. &
           ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
     INFO = 1
  ELSE if (  ( .NOT.NOTB                 ).AND. &
           ( .NOT.CONJB                ).AND. &
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
     call XERBLA( 'CGEMM ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if (  ( M == 0 ).OR.( N == 0 ).OR. &
      ( ( ( ALPHA == ZERO ).OR.( K == 0 ) ).AND.( BETA == ONE ) ) ) &
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
     ELSE if (  CONJA )THEN
!
!           Form  C := alpha*conjg( A' )*B + beta*C.
!
        DO 120, J = 1, N
           DO 110, I = 1, M
              TEMP = ZERO
              DO 100, L = 1, K
                 TEMP = TEMP + CONJG( A( L, I ) )*B( L, J )
  100             CONTINUE
              if (  BETA == ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
  110          CONTINUE
  120       CONTINUE
     ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
        DO 150, J = 1, N
           DO 140, I = 1, M
              TEMP = ZERO
              DO 130, L = 1, K
                 TEMP = TEMP + A( L, I )*B( L, J )
  130             CONTINUE
              if (  BETA == ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
  140          CONTINUE
  150       CONTINUE
     end if
  ELSE if (  NOTA )THEN
     if (  CONJB )THEN
!
!           Form  C := alpha*A*conjg( B' ) + beta*C.
!
        DO 200, J = 1, N
           if (  BETA == ZERO )THEN
              DO 160, I = 1, M
                 C( I, J ) = ZERO
  160             CONTINUE
           ELSE if (  BETA /= ONE )THEN
              DO 170, I = 1, M
                 C( I, J ) = BETA*C( I, J )
  170             CONTINUE
           end if
           DO 190, L = 1, K
              if (  B( J, L ) /= ZERO )THEN
                 TEMP = ALPHA*CONJG( B( J, L ) )
                 DO 180, I = 1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
  180                CONTINUE
              end if
  190          CONTINUE
  200       CONTINUE
     ELSE
!
!           Form  C := alpha*A*B'          + beta*C
!
        DO 250, J = 1, N
           if (  BETA == ZERO )THEN
              DO 210, I = 1, M
                 C( I, J ) = ZERO
  210             CONTINUE
           ELSE if (  BETA /= ONE )THEN
              DO 220, I = 1, M
                 C( I, J ) = BETA*C( I, J )
  220             CONTINUE
           end if
           DO 240, L = 1, K
              if (  B( J, L ) /= ZERO )THEN
                 TEMP = ALPHA*B( J, L )
                 DO 230, I = 1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
  230                CONTINUE
              end if
  240          CONTINUE
  250       CONTINUE
     end if
  ELSE if (  CONJA )THEN
     if (  CONJB )THEN
!
!           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
!
        DO 280, J = 1, N
           DO 270, I = 1, M
              TEMP = ZERO
              DO 260, L = 1, K
                 TEMP = TEMP + CONJG( A( L, I ) )*CONJG( B( J, L ) )
  260             CONTINUE
              if (  BETA == ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
  270          CONTINUE
  280       CONTINUE
     ELSE
!
!           Form  C := alpha*conjg( A' )*B' + beta*C
!
        DO 310, J = 1, N
           DO 300, I = 1, M
              TEMP = ZERO
              DO 290, L = 1, K
                 TEMP = TEMP + CONJG( A( L, I ) )*B( J, L )
  290             CONTINUE
              if (  BETA == ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
  300          CONTINUE
  310       CONTINUE
     end if
  ELSE
     if (  CONJB )THEN
!
!           Form  C := alpha*A'*conjg( B' ) + beta*C
!
        DO 340, J = 1, N
           DO 330, I = 1, M
              TEMP = ZERO
              DO 320, L = 1, K
                 TEMP = TEMP + A( L, I )*CONJG( B( J, L ) )
  320             CONTINUE
              if (  BETA == ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
  330          CONTINUE
  340       CONTINUE
     ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
        DO 370, J = 1, N
           DO 360, I = 1, M
              TEMP = ZERO
              DO 350, L = 1, K
                 TEMP = TEMP + A( L, I )*B( J, L )
  350             CONTINUE
              if (  BETA == ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
  360          CONTINUE
  370       CONTINUE
     end if
  end if
!
  return
!
!     End of CGEMM .
!
end
