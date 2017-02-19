subroutine CTPSV (UPLO, TRANS, DIAG, N, AP, X, INCX)
!
!! CTPSV solves a triangular system of linear equations.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B4
!***TYPE      COMPLEX (STPSV-S, DTPSV-D, CTPSV-C)
!***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J. J., (ANL)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!  CTPSV  solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,   or   conjg( A')*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular matrix, supplied in packed form.
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
!              TRANS = 'C' or 'c'   conjg( A' )*x = b.
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
!  AP     - COMPLEX          array of DIMENSION at least
!           ( ( n*( n + 1 ) )/2 ).
!           Before entry with  UPLO = 'U' or 'u', the array AP must
!           contain the upper triangular matrix packed sequentially,
!           column by column, so that AP( 1 ) contains a( 1, 1 ),
!           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
!           respectively, and so on.
!           Before entry with UPLO = 'L' or 'l', the array AP must
!           contain the lower triangular matrix packed sequentially,
!           column by column, so that AP( 1 ) contains a( 1, 1 ),
!           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
!           respectively, and so on.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  X      - COMPLEX          array of dimension at least
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
!***END PROLOGUE  CTPSV
!     .. Scalar Arguments ..
  INTEGER            INCX, N
  CHARACTER*1        DIAG, TRANS, UPLO
!     .. Array Arguments ..
  COMPLEX            AP( * ), X( * )
!     .. Parameters ..
  COMPLEX            ZERO
  PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     .. Local Scalars ..
  COMPLEX            TEMP
  INTEGER            I, INFO, IX, J, JX, K, KK, KX
  LOGICAL            NOCONJ, NOUNIT
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
  INTRINSIC          CONJG
!***FIRST EXECUTABLE STATEMENT  CTPSV
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
  ELSE if (  INCX == 0 )THEN
     INFO = 7
  end if
  if (  INFO /= 0 )THEN
     call XERBLA( 'CTPSV ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if (  N == 0 ) &
     return
!
  NOCONJ = LSAME( TRANS, 'T' )
  NOUNIT = LSAME( DIAG , 'N' )
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
!     Start the operations. In this version the elements of AP are
!     accessed sequentially with one pass through AP.
!
  if (  LSAME( TRANS, 'N' ) )THEN
!
!        Form  x := inv( A )*x.
!
     if (  LSAME( UPLO, 'U' ) )THEN
        KK = ( N*( N + 1 ) )/2
        if (  INCX == 1 )THEN
           DO 20, J = N, 1, -1
              if (  X( J ) /= ZERO )THEN
                 if (  NOUNIT ) &
                    X( J ) = X( J )/AP( KK )
                 TEMP = X( J )
                 K    = KK     - 1
                 DO 10, I = J - 1, 1, -1
                    X( I ) = X( I ) - TEMP*AP( K )
                    K      = K      - 1
   10                CONTINUE
              end if
              KK = KK - J
   20          CONTINUE
        ELSE
           JX = KX + ( N - 1 )*INCX
           DO 40, J = N, 1, -1
              if (  X( JX ) /= ZERO )THEN
                 if (  NOUNIT ) &
                    X( JX ) = X( JX )/AP( KK )
                 TEMP = X( JX )
                 IX   = JX
                 DO 30, K = KK - 1, KK - J + 1, -1
                    IX      = IX      - INCX
                    X( IX ) = X( IX ) - TEMP*AP( K )
   30                CONTINUE
              end if
              JX = JX - INCX
              KK = KK - J
   40          CONTINUE
        end if
     ELSE
        KK = 1
        if (  INCX == 1 )THEN
           DO 60, J = 1, N
              if (  X( J ) /= ZERO )THEN
                 if (  NOUNIT ) &
                    X( J ) = X( J )/AP( KK )
                 TEMP = X( J )
                 K    = KK     + 1
                 DO 50, I = J + 1, N
                    X( I ) = X( I ) - TEMP*AP( K )
                    K      = K      + 1
   50                CONTINUE
              end if
              KK = KK + ( N - J + 1 )
   60          CONTINUE
        ELSE
           JX = KX
           DO 80, J = 1, N
              if (  X( JX ) /= ZERO )THEN
                 if (  NOUNIT ) &
                    X( JX ) = X( JX )/AP( KK )
                 TEMP = X( JX )
                 IX   = JX
                 DO 70, K = KK + 1, KK + N - J
                    IX      = IX      + INCX
                    X( IX ) = X( IX ) - TEMP*AP( K )
   70                CONTINUE
              end if
              JX = JX + INCX
              KK = KK + ( N - J + 1 )
   80          CONTINUE
        end if
     end if
  ELSE
!
!        Form  x := inv( A' )*x  or  x := inv( conjg( A' ) )*x.
!
     if (  LSAME( UPLO, 'U' ) )THEN
        KK = 1
        if (  INCX == 1 )THEN
           DO 110, J = 1, N
              TEMP = X( J )
              K    = KK
              if (  NOCONJ )THEN
                 DO 90, I = 1, J - 1
                    TEMP = TEMP - AP( K )*X( I )
                    K    = K    + 1
   90                CONTINUE
                 if (  NOUNIT ) &
                    TEMP = TEMP/AP( KK + J - 1 )
              ELSE
                 DO 100, I = 1, J - 1
                    TEMP = TEMP - CONJG( AP( K ) )*X( I )
                    K    = K    + 1
  100                CONTINUE
                 if (  NOUNIT ) &
                    TEMP = TEMP/CONJG( AP( KK + J - 1 ) )
              end if
              X( J ) = TEMP
              KK     = KK   + J
  110          CONTINUE
        ELSE
           JX = KX
           DO 140, J = 1, N
              TEMP = X( JX )
              IX   = KX
              if (  NOCONJ )THEN
                 DO 120, K = KK, KK + J - 2
                    TEMP = TEMP - AP( K )*X( IX )
                    IX   = IX   + INCX
  120                CONTINUE
                 if (  NOUNIT ) &
                    TEMP = TEMP/AP( KK + J - 1 )
              ELSE
                 DO 130, K = KK, KK + J - 2
                    TEMP = TEMP - CONJG( AP( K ) )*X( IX )
                    IX   = IX   + INCX
  130                CONTINUE
                 if (  NOUNIT ) &
                    TEMP = TEMP/CONJG( AP( KK + J - 1 ) )
              end if
              X( JX ) = TEMP
              JX      = JX   + INCX
              KK      = KK   + J
  140          CONTINUE
        end if
     ELSE
        KK = ( N*( N + 1 ) )/2
        if (  INCX == 1 )THEN
           DO 170, J = N, 1, -1
              TEMP = X( J )
              K    = KK
              if (  NOCONJ )THEN
                 DO 150, I = N, J + 1, -1
                    TEMP = TEMP - AP( K )*X( I )
                    K    = K    - 1
  150                CONTINUE
                 if (  NOUNIT ) &
                    TEMP = TEMP/AP( KK - N + J )
              ELSE
                 DO 160, I = N, J + 1, -1
                    TEMP = TEMP - CONJG( AP( K ) )*X( I )
                    K    = K    - 1
  160                CONTINUE
                 if (  NOUNIT ) &
                    TEMP = TEMP/CONJG( AP( KK - N + J ) )
              end if
              X( J ) = TEMP
              KK     = KK   - ( N - J + 1 )
  170          CONTINUE
        ELSE
           KX = KX + ( N - 1 )*INCX
           JX = KX
           DO 200, J = N, 1, -1
              TEMP = X( JX )
              IX   = KX
              if (  NOCONJ )THEN
                 DO 180, K = KK, KK - ( N - ( J + 1 ) ), -1
                    TEMP = TEMP - AP( K )*X( IX )
                    IX   = IX   - INCX
  180                CONTINUE
                 if (  NOUNIT ) &
                    TEMP = TEMP/AP( KK - N + J )
              ELSE
                 DO 190, K = KK, KK - ( N - ( J + 1 ) ), -1
                    TEMP = TEMP - CONJG( AP( K ) )*X( IX )
                    IX   = IX   - INCX
  190                CONTINUE
                 if (  NOUNIT ) &
                    TEMP = TEMP/CONJG( AP( KK - N + J ) )
              end if
              X( JX ) = TEMP
              JX      = JX   - INCX
              KK      = KK   - ( N - J + 1 )
  200          CONTINUE
        end if
     end if
  end if
!
  return
!
!     End of CTPSV .
!
end
