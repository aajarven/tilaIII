subroutine DQFORM (M, N, Q, LDQ, WA)
!
!! DQFORM explicitly forms the Q matrix of an implicit QR factorization.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (QFORM-S, DQFORM-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine proceeds from the computed QR factorization of
!     an M by N matrix A to accumulate the M by M orthogonal matrix
!     Q from its factored form.
!
!     The subroutine statement is
!
!       SUBROUTINE DQFORM(M,N,Q,LDQ,WA)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of A and the order of Q.
!
!       N is a positive integer input variable set to the number
!         of columns of A.
!
!       Q is an M by M array. On input the full lower trapezoid in
!         the first MIN(M,N) columns of Q contains the factored form.
!         On output Q has been accumulated into a square matrix.
!
!       LDQ is a positive integer input variable not less than M
!         which specifies the leading dimension of the array Q.
!
!       WA is a work array of length M.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DQFORM
  INTEGER I, J, JM1, K, L, LDQ, M, MINMN, N, NP1
  DOUBLE PRECISION ONE, Q(LDQ,*), SUM, TEMP, WA(*), ZERO
  SAVE ONE, ZERO
  DATA ONE,ZERO /1.0D0,0.0D0/
!
!     ZERO OUT UPPER TRIANGLE OF Q IN THE FIRST MIN(M,N) COLUMNS.
!
!***FIRST EXECUTABLE STATEMENT  DQFORM
  MINMN = MIN(M,N)
  if (MINMN  <  2) go to 30
  DO 20 J = 2, MINMN
     JM1 = J - 1
     DO 10 I = 1, JM1
        Q(I,J) = ZERO
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
!
!     INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.
!
  NP1 = N + 1
  if (M  <  NP1) go to 60
  DO 50 J = NP1, M
     DO 40 I = 1, M
        Q(I,J) = ZERO
   40       CONTINUE
     Q(J,J) = ONE
   50    CONTINUE
   60 CONTINUE
!
!     ACCUMULATE Q FROM ITS FACTORED FORM.
!
  DO 120 L = 1, MINMN
     K = MINMN - L + 1
     DO 70 I = K, M
        WA(I) = Q(I,K)
        Q(I,K) = ZERO
   70       CONTINUE
     Q(K,K) = ONE
     if (WA(K)  ==  ZERO) go to 110
     DO 100 J = K, M
        SUM = ZERO
        DO 80 I = K, M
           SUM = SUM + Q(I,J)*WA(I)
   80          CONTINUE
        TEMP = SUM/WA(K)
        DO 90 I = K, M
           Q(I,J) = Q(I,J) - TEMP*WA(I)
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
  120    CONTINUE
  return
!
!     LAST CARD OF SUBROUTINE DQFORM.
!
end
