subroutine CHERK (UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC)
!
!! CHERK performs Hermitian rank k update of a complex Hermitian matrix.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B6
!***TYPE      COMPLEX (SHERK-S, DHERK-D, CHERK-C)
!***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S. (NAG)
!***DESCRIPTION
!
!  CHERK  performs one of the hermitian rank k operations
!
!     C := alpha*A*conjg( A' ) + beta*C,
!
!  or
!
!     C := alpha*conjg( A' )*A + beta*C,
!
!  where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
!  matrix and  A  is an  n by k  matrix in the  first case and a  k by n
!  matrix in the second case.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!           triangular  part  of the  array  C  is to be  referenced  as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry,  TRANS  specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   C := alpha*A*conjg( A' ) + beta*C.
!
!              TRANS = 'C' or 'c'   C := alpha*conjg( A' )*A + beta*C.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N specifies the order of the matrix C.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
!           of  columns   of  the   matrix   A,   and  on   entry   with
!           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
!           matrix A.  K must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by n  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!           then  LDA must be at least  max( 1, n ), otherwise  LDA must
!           be at least  max( 1, k ).
!           Unchanged on exit.
!
!  BETA   - REAL            .
!           On entry, BETA specifies the scalar beta.
!           Unchanged on exit.
!
!  C      - COMPLEX          array of DIMENSION ( LDC, n ).
!           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
!           upper triangular part of the array C must contain the upper
!           triangular part  of the  hermitian matrix  and the strictly
!           lower triangular part of C is not referenced.  On exit, the
!           upper triangular part of the array  C is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
!           lower triangular part of the array C must contain the lower
!           triangular part  of the  hermitian matrix  and the strictly
!           upper triangular part of C is not referenced.  On exit, the
!           lower triangular part of the array  C is overwritten by the
!           lower triangular part of the updated matrix.
!           Note that the imaginary parts of the diagonal elements need
!           not be set,  they are assumed to be zero,  and on exit they
!           are set to zero.
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, n ).
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
!***END PROLOGUE  CHERK
!     .. Scalar Arguments ..
  CHARACTER*1        UPLO, TRANS
  INTEGER            N, K, LDA, LDC
  REAL               ALPHA, BETA
!     .. Array Arguments ..
  COMPLEX            A( LDA, * ), C( LDC, * )
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          CMPLX, CONJG, MAX, REAL
!     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            I, INFO, J, L, NROWA
  REAL               RTEMP
  COMPLEX            TEMP
!     .. Parameters ..
  REAL               ONE ,         ZERO
  PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!***FIRST EXECUTABLE STATEMENT  CHERK
!
!     Test the input parameters.
!
  if (  LSAME( TRANS, 'N' ) )THEN
     NROWA = N
  ELSE
     NROWA = K
  end if
  UPPER = LSAME( UPLO, 'U' )
!
  INFO = 0
  if (       ( .NOT.UPPER               ).AND. &
           ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
     INFO = 1
  ELSE if (  ( .NOT.LSAME( TRANS, 'N' ) ).AND. &
           ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
     INFO = 2
  ELSE if (  N   < 0               )THEN
     INFO = 3
  ELSE if (  K   < 0               )THEN
     INFO = 4
  ELSE if (  LDA < MAX( 1, NROWA ) )THEN
     INFO = 7
  ELSE if (  LDC < MAX( 1, N     ) )THEN
     INFO = 10
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'CHERK ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if (  ( N == 0 ).OR. &
      ( ( ( ALPHA == ZERO ).OR.( K == 0 ) ).AND.( BETA == ONE ) ) ) &
     return
!
!     And when  alpha.eq.zero.
!
  if (  ALPHA == ZERO )THEN
     if (  UPPER )THEN
        if (  BETA == ZERO )THEN
           DO 20, J = 1, N
              DO 10, I = 1, J
                 C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
        ELSE
           DO 40, J = 1, N
              DO 30, I = 1, J - 1
                 C( I, J ) = BETA*C( I, J )
   30             CONTINUE
              C( J, J ) = BETA*REAL( C( J, J ) )
   40          CONTINUE
        end if
     ELSE
        if (  BETA == ZERO )THEN
           DO 60, J = 1, N
              DO 50, I = J, N
                 C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
        ELSE
           DO 80, J = 1, N
              C( J, J ) = BETA*REAL( C( J, J ) )
              DO 70, I = J + 1, N
                 C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
        end if
     end if
     return
  end if
!
!     Start the operations.
!
  if (  LSAME( TRANS, 'N' ) )THEN
!
!        Form  C := alpha*A*conjg( A' ) + beta*C.
!
     if (  UPPER )THEN
        DO 130, J = 1, N
           if (  BETA == ZERO )THEN
              DO 90, I = 1, J
                 C( I, J ) = ZERO
   90             CONTINUE
           ELSE if (  BETA /= ONE )THEN
              DO 100, I = 1, J - 1
                 C( I, J ) = BETA*C( I, J )
  100             CONTINUE
              C( J, J ) = BETA*REAL( C( J, J ) )
           end if
           DO 120, L = 1, K
              if (  A( J, L ) /= CMPLX( ZERO ) )THEN
                 TEMP = ALPHA*CONJG( A( J, L ) )
                 DO 110, I = 1, J - 1
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
  110                CONTINUE
                 C( J, J ) = REAL( C( J, J )      ) + &
                             REAL( TEMP*A( I, L ) )
              end if
  120          CONTINUE
  130       CONTINUE
     ELSE
        DO 180, J = 1, N
           if (  BETA == ZERO )THEN
              DO 140, I = J, N
                 C( I, J ) = ZERO
  140             CONTINUE
           ELSE if (  BETA /= ONE )THEN
              C( J, J ) = BETA*REAL( C( J, J ) )
              DO 150, I = J + 1, N
                 C( I, J ) = BETA*C( I, J )
  150             CONTINUE
           end if
           DO 170, L = 1, K
              if (  A( J, L ) /= CMPLX( ZERO ) )THEN
                 TEMP      = ALPHA*CONJG( A( J, L ) )
                 C( J, J ) = REAL( C( J, J )      )   + &
                             REAL( TEMP*A( J, L ) )
                 DO 160, I = J + 1, N
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
  160                CONTINUE
              end if
  170          CONTINUE
  180       CONTINUE
     end if
  ELSE
!
!        Form  C := alpha*conjg( A' )*A + beta*C.
!
     if (  UPPER )THEN
        DO 220, J = 1, N
           DO 200, I = 1, J - 1
              TEMP = ZERO
              DO 190, L = 1, K
                 TEMP = TEMP + CONJG( A( L, I ) )*A( L, J )
  190             CONTINUE
              if (  BETA == ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
  200          CONTINUE
           RTEMP = ZERO
           DO 210, L = 1, K
              RTEMP = RTEMP + CONJG( A( L, J ) )*A( L, J )
  210          CONTINUE
           if (  BETA == ZERO )THEN
              C( J, J ) = ALPHA*RTEMP
           ELSE
              C( J, J ) = ALPHA*RTEMP + BETA*REAL( C( J, J ) )
           end if
  220       CONTINUE
     ELSE
        DO 260, J = 1, N
           RTEMP = ZERO
           DO 230, L = 1, K
              RTEMP = RTEMP + CONJG( A( L, J ) )*A( L, J )
  230          CONTINUE
           if (  BETA == ZERO )THEN
              C( J, J ) = ALPHA*RTEMP
           ELSE
              C( J, J ) = ALPHA*RTEMP + BETA*REAL( C( J, J ) )
           end if
           DO 250, I = J + 1, N
              TEMP = ZERO
              DO 240, L = 1, K
                 TEMP = TEMP + CONJG( A( L, I ) )*A( L, J )
  240             CONTINUE
              if (  BETA == ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
  250          CONTINUE
  260       CONTINUE
     end if
  end if
!
  return
!
!     End of CHERK .
!
end
