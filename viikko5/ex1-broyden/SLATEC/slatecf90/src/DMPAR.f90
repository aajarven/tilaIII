subroutine DMPAR (N, R, LDR, IPVT, DIAG, QTB, DELTA, PAR, X, &
     SIGMA, WA1, WA2)
!
!! DMPAR is subsidiary to DNLS1 and DNLS1E.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNLS1 and DNLS1E
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (LMPAR-S, DMPAR-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   **** Double Precision version of LMPAR ****
!
!     Given an M by N matrix A, an N by N nonsingular DIAGONAL
!     matrix D, an M-vector B, and a positive number DELTA,
!     the problem is to determine a value for the parameter
!     PAR such that if X solves the system
!
!           A*X = B ,     SQRT(PAR)*D*X = 0 ,
!
!     in the least squares sense, and DXNORM is the Euclidean
!     norm of D*X, then either PAR is zero and
!
!           (DXNORM-DELTA)  <=  0.1*DELTA ,
!
!     or PAR is positive and
!
!           ABS(DXNORM-DELTA)  <=  0.1*DELTA .
!
!     This subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     QR factorization, with column pivoting, of A. That is, if
!     A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!     columns, and R is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then DMPAR expects
!     the full upper triangle of R, the permutation matrix P,
!     and the first N components of (Q TRANSPOSE)*B. On output
!     DMPAR also provides an upper triangular matrix S such that
!
!            T   T                   T
!           P *(A *A + PAR*D*D)*P = S *S .
!
!     S is employed within DMPAR and may be of separate interest.
!
!     Only a few iterations are generally needed for convergence
!     of the algorithm. If, however, the limit of 10 iterations
!     is reached, then the output PAR will contain the best
!     value obtained so far.
!
!     The subroutine statement is
!
!       SUBROUTINE DMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SIGMA,
!                        WA1,WA2)
!
!     where
!
!       N is a positive integer input variable set to the order of R.
!
!       R is an N by N array. On input the full upper triangle
!         must contain the full upper triangle of the matrix R.
!         On output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix S.
!
!       LDR is a positive integer input variable not less than N
!         which specifies the leading dimension of the array R.
!
!       IPVT is an integer input array of length N which defines the
!         permutation matrix P such that A*P = Q*R. Column J of P
!         is column IPVT(J) of the identity matrix.
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
!       PAR is a nonnegative variable. On input PAR contains an
!         initial estimate of the Levenberg-Marquardt parameter.
!         On output PAR contains the final estimate.
!
!       X is an output array of length N which contains the least
!         squares solution of the system A*X = B, SQRT(PAR)*D*X = 0,
!         for the output PAR.
!
!       SIGMA is an output array of length N which contains the
!         diagonal elements of the upper triangular matrix S.
!
!       WA1 and WA2 are work arrays of length N.
!
!***SEE ALSO  DNLS1, DNLS1E
!***ROUTINES CALLED  D1MACH, DENORM, DQRSLV
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DMPAR
  INTEGER N,LDR
  INTEGER IPVT(*)
  DOUBLE PRECISION DELTA,PAR
  DOUBLE PRECISION R(LDR,*),DIAG(*),QTB(*),X(*),SIGMA(*),WA1(*), &
   WA2(*)
  INTEGER I,ITER,J,JM1,JP1,K,L,NSING
  DOUBLE PRECISION DXNORM,DWARF,FP,GNORM,PARC,PARL,PARU,P1,P001, &
   SUM,TEMP,ZERO
  DOUBLE PRECISION D1MACH,DENORM
  SAVE P1, P001, ZERO
  DATA P1,P001,ZERO /1.0D-1,1.0D-3,0.0D0/
!***FIRST EXECUTABLE STATEMENT  DMPAR
  DWARF = D1MACH(1)
!
!     COMPUTE AND STORE IN X THE GAUSS-NEWTON DIRECTION. if THE
!     JACOBIAN IS RANK-DEFICIENT, OBTAIN A LEAST SQUARES SOLUTION.
!
  NSING = N
  DO 10 J = 1, N
     WA1(J) = QTB(J)
     if (R(J,J)  ==  ZERO .AND. NSING  ==  N) NSING = J - 1
     if (NSING  <  N) WA1(J) = ZERO
   10    CONTINUE
  if (NSING  <  1) go to 50
  DO 40 K = 1, NSING
     J = NSING - K + 1
     WA1(J) = WA1(J)/R(J,J)
     TEMP = WA1(J)
     JM1 = J - 1
     if (JM1  <  1) go to 30
     DO 20 I = 1, JM1
        WA1(I) = WA1(I) - R(I,J)*TEMP
   20       CONTINUE
   30    CONTINUE
   40    CONTINUE
   50 CONTINUE
  DO 60 J = 1, N
     L = IPVT(J)
     X(L) = WA1(J)
   60    CONTINUE
!
!     INITIALIZE THE ITERATION COUNTER.
!     EVALUATE THE FUNCTION AT THE ORIGIN, AND TEST
!     FOR ACCEPTANCE OF THE GAUSS-NEWTON DIRECTION.
!
  ITER = 0
  DO 70 J = 1, N
     WA2(J) = DIAG(J)*X(J)
   70    CONTINUE
  DXNORM = DENORM(N,WA2)
  FP = DXNORM - DELTA
  if (FP  <=  P1*DELTA) go to 220
!
!     if THE JACOBIAN IS NOT RANK DEFICIENT, THE NEWTON
!     STEP PROVIDES A LOWER BOUND, PARL, FOR THE ZERO OF
!     THE FUNCTION. OTHERWISE SET THIS BOUND TO ZERO.
!
  PARL = ZERO
  if (NSING  <  N) go to 120
  DO 80 J = 1, N
     L = IPVT(J)
     WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
   80    CONTINUE
  DO 110 J = 1, N
     SUM = ZERO
     JM1 = J - 1
     if (JM1  <  1) go to 100
     DO 90 I = 1, JM1
        SUM = SUM + R(I,J)*WA1(I)
   90       CONTINUE
  100    CONTINUE
     WA1(J) = (WA1(J) - SUM)/R(J,J)
  110    CONTINUE
  TEMP = DENORM(N,WA1)
  PARL = ((FP/DELTA)/TEMP)/TEMP
  120 CONTINUE
!
!     CALCULATE AN UPPER BOUND, PARU, FOR THE ZERO OF THE FUNCTION.
!
  DO 140 J = 1, N
     SUM = ZERO
     DO 130 I = 1, J
        SUM = SUM + R(I,J)*QTB(I)
  130       CONTINUE
     L = IPVT(J)
     WA1(J) = SUM/DIAG(L)
  140    CONTINUE
  GNORM = DENORM(N,WA1)
  PARU = GNORM/DELTA
  if (PARU  ==  ZERO) PARU = DWARF/MIN(DELTA,P1)
!
!     if THE INPUT PAR LIES OUTSIDE OF THE INTERVAL (PARL,PARU),
!     SET PAR TO THE CLOSER ENDPOINT.
!
  PAR = MAX(PAR,PARL)
  PAR = MIN(PAR,PARU)
  if (PAR  ==  ZERO) PAR = GNORM/DXNORM
!
!     BEGINNING OF AN ITERATION.
!
  150 CONTINUE
     ITER = ITER + 1
!
!        EVALUATE THE FUNCTION AT THE CURRENT VALUE OF PAR.
!
     if (PAR  ==  ZERO) PAR = MAX(DWARF,P001*PARU)
     TEMP = SQRT(PAR)
     DO 160 J = 1, N
        WA1(J) = TEMP*DIAG(J)
  160       CONTINUE
     call DQRSLV(N,R,LDR,IPVT,WA1,QTB,X,SIGMA,WA2)
     DO 170 J = 1, N
        WA2(J) = DIAG(J)*X(J)
  170       CONTINUE
     DXNORM = DENORM(N,WA2)
     TEMP = FP
     FP = DXNORM - DELTA
!
!        if THE FUNCTION IS SMALL ENOUGH, ACCEPT THE CURRENT VALUE
!        OF PAR. ALSO TEST FOR THE EXCEPTIONAL CASES WHERE PARL
!        IS ZERO OR THE NUMBER OF ITERATIONS HAS REACHED 10.
!
     if (ABS(FP)  <=  P1*DELTA &
         .OR. PARL  ==  ZERO .AND. FP  <=  TEMP &
              .AND. TEMP  <  ZERO .OR. ITER  ==  10) go to 220
!
!        COMPUTE THE NEWTON CORRECTION.
!
     DO 180 J = 1, N
        L = IPVT(J)
        WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
  180       CONTINUE
     DO 210 J = 1, N
        WA1(J) = WA1(J)/SIGMA(J)
        TEMP = WA1(J)
        JP1 = J + 1
        if (N  <  JP1) go to 200
        DO 190 I = JP1, N
           WA1(I) = WA1(I) - R(I,J)*TEMP
  190          CONTINUE
  200       CONTINUE
  210       CONTINUE
     TEMP = DENORM(N,WA1)
     PARC = ((FP/DELTA)/TEMP)/TEMP
!
!        DEPENDING ON THE SIGN OF THE FUNCTION, UPDATE PARL OR PARU.
!
     if (FP  >  ZERO) PARL = MAX(PARL,PAR)
     if (FP  <  ZERO) PARU = MIN(PARU,PAR)
!
!        COMPUTE AN IMPROVED ESTIMATE FOR PAR.
!
     PAR = MAX(PARL,PAR+PARC)
!
!        END OF AN ITERATION.
!
     go to 150
  220 CONTINUE
!
!     TERMINATION.
!
  if (ITER  ==  0) PAR = ZERO
  return
!
!     LAST CARD OF SUBROUTINE DMPAR.
!
end
