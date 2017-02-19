subroutine DOGLEG (N, R, LR, DIAG, QTB, DELTA, X, WA1, WA2)
!
!! DOGLEG is subsidiary to SNSQ and SNSQE.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (DOGLEG-S, DDOGLG-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Given an M by N matrix A, an N by N nonsingular DIAGONAL
!     matrix D, an M-vector B, and a positive number DELTA, the
!     problem is to determine the convex combination X of the
!     Gauss-Newton and scaled gradient directions that minimizes
!     (A*X - B) in the least squares sense, subject to the
!     restriction that the Euclidean norm of D*X be at most DELTA.
!
!     This subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     QR factorization of A. That is, if A = Q*R, where Q has
!     orthogonal columns and R is an upper triangular matrix,
!     then DOGLEG expects the full upper triangle of R and
!     the first N components of (Q TRANSPOSE)*B.
!
!     The subroutine statement is
!
!       SUBROUTINE DOGLEG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
!
!     where
!
!       N is a positive integer input variable set to the order of R.
!
!       R is an input array of length LR which must contain the upper
!         triangular matrix R stored by rows.
!
!       LR is a positive integer input variable not less than
!         (N*(N+1))/2.
!
!       DIAG is an input array of length N which must contain the
!         diagonal elements of the matrix D.
!
!       QTB is an input array of length N which must contain the first
!         N elements of the vector (Q TRANSPOSE)*B.
!
!       DELTA is a positive input variable which specifies an upper
!         bound on the Euclidean norm of D*X.
!
!       X is an output array of length N which contains the desired
!         convex combination of the Gauss-Newton direction and the
!         scaled gradient direction.
!
!       WA1 and WA2 are work arrays of length N.
!
!***SEE ALSO  SNSQ, SNSQE
!***ROUTINES CALLED  ENORM, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DOGLEG
  INTEGER N,LR
  REAL DELTA
  REAL R(LR),DIAG(*),QTB(*),X(*),WA1(*),WA2(*)
  INTEGER I,J,JJ,JP1,K,L
  REAL ALPHA,BNORM,EPSMCH,GNORM,ONE,QNORM,SGNORM,SUM,TEMP,ZERO
  REAL R1MACH,ENORM
  SAVE ONE, ZERO
  DATA ONE,ZERO /1.0E0,0.0E0/
!***FIRST EXECUTABLE STATEMENT  DOGLEG
  EPSMCH = R1MACH(4)
!
!     FIRST, CALCULATE THE GAUSS-NEWTON DIRECTION.
!
  JJ = (N*(N + 1))/2 + 1
  DO 50 K = 1, N
     J = N - K + 1
     JP1 = J + 1
     JJ = JJ - K
     L = JJ + 1
     SUM = ZERO
     if (N  <  JP1) go to 20
     DO 10 I = JP1, N
        SUM = SUM + R(L)*X(I)
        L = L + 1
   10       CONTINUE
   20    CONTINUE
     TEMP = R(JJ)
     if (TEMP  /=  ZERO) go to 40
     L = J
     DO 30 I = 1, J
        TEMP = MAX(TEMP,ABS(R(L)))
        L = L + N - I
   30       CONTINUE
     TEMP = EPSMCH*TEMP
     if (TEMP  ==  ZERO) TEMP = EPSMCH
   40    CONTINUE
     X(J) = (QTB(J) - SUM)/TEMP
   50    CONTINUE
!
!     TEST WHETHER THE GAUSS-NEWTON DIRECTION IS ACCEPTABLE.
!
  DO 60 J = 1, N
     WA1(J) = ZERO
     WA2(J) = DIAG(J)*X(J)
   60    CONTINUE
  QNORM = ENORM(N,WA2)
  if (QNORM  <=  DELTA) go to 140
!
!     THE GAUSS-NEWTON DIRECTION IS NOT ACCEPTABLE.
!     NEXT, CALCULATE THE SCALED GRADIENT DIRECTION.
!
  L = 1
  DO 80 J = 1, N
     TEMP = QTB(J)
     DO 70 I = J, N
        WA1(I) = WA1(I) + R(L)*TEMP
        L = L + 1
   70       CONTINUE
     WA1(J) = WA1(J)/DIAG(J)
   80    CONTINUE
!
!     CALCULATE THE NORM OF THE SCALED GRADIENT DIRECTION,
!     NORMALIZE, AND RESCALE THE GRADIENT.
!
  GNORM = ENORM(N,WA1)
  SGNORM = ZERO
  ALPHA = DELTA/QNORM
  if (GNORM  ==  ZERO) go to 120
  DO 90 J = 1, N
     WA1(J) = (WA1(J)/GNORM)/DIAG(J)
   90    CONTINUE
!
!     CALCULATE THE POINT ALONG THE SCALED GRADIENT
!     AT WHICH THE QUADRATIC IS MINIMIZED.
!
  L = 1
  DO 110 J = 1, N
     SUM = ZERO
     DO 100 I = J, N
        SUM = SUM + R(L)*WA1(I)
        L = L + 1
  100       CONTINUE
     WA2(J) = SUM
  110    CONTINUE
  TEMP = ENORM(N,WA2)
  SGNORM = (GNORM/TEMP)/TEMP
!
!     TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE.
!
  ALPHA = ZERO
  if (SGNORM  >=  DELTA) go to 120
!
!     THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE.
!     FINALLY, CALCULATE THE POINT ALONG THE DOGLEG
!     AT WHICH THE QUADRATIC IS MINIMIZED.
!
  BNORM = ENORM(N,QTB)
  TEMP = (BNORM/GNORM)*(BNORM/QNORM)*(SGNORM/DELTA)
  TEMP = TEMP - (DELTA/QNORM)*(SGNORM/DELTA)**2 &
         + SQRT((TEMP-(DELTA/QNORM))**2 &
                +(ONE-(DELTA/QNORM)**2)*(ONE-(SGNORM/DELTA)**2))
  ALPHA = ((DELTA/QNORM)*(ONE - (SGNORM/DELTA)**2))/TEMP
  120 CONTINUE
!
!     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON
!     DIRECTION AND THE SCALED GRADIENT DIRECTION.
!
  TEMP = (ONE - ALPHA)*MIN(SGNORM,DELTA)
  DO 130 J = 1, N
     X(J) = TEMP*WA1(J) + ALPHA*X(J)
  130    CONTINUE
  140 CONTINUE
  return
!
!     LAST CARD OF SUBROUTINE DOGLEG.
!
end
