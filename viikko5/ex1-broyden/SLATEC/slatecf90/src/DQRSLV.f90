subroutine DQRSLV (N, R, LDR, IPVT, DIAG, QTB, X, SIGMA, WA)
!
!! DQRSLV solves a least squares problem.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNLS1 and DNLS1E
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (QRSOLV-S, DQRSLV-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!  **** Double Precision version of QRSOLV ****
!
!     Given an M by N matrix A, an N by N diagonal matrix D,
!     and an M-vector B, the problem is to determine an X which
!     solves the system
!
!           A*X = B ,     D*X = 0 ,
!
!     in the least squares sense.
!
!     This subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     QR factorization, with column pivoting, of A. That is, if
!     A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!     columns, and R is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then DQRSLV expects
!     the full upper triangle of R, the permutation matrix P,
!     and the first N components of (Q TRANSPOSE)*B. The system
!     A*X = B, D*X = 0, is then equivalent to
!
!                  T       T
!           R*Z = Q *B ,  P *D*P*Z = 0 ,
!
!     where X = P*Z. If this system does not have full rank,
!     then a least squares solution is obtained. On output DQRSLV
!     also provides an upper triangular matrix S such that
!
!            T   T               T
!           P *(A *A + D*D)*P = S *S .
!
!     S is computed within DQRSLV and may be of separate interest.
!
!     The subroutine statement is
!
!       SUBROUTINE DQRSLV(N,R,LDR,IPVT,DIAG,QTB,X,SIGMA,WA)
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
!       X is an output array of length N which contains the least
!         squares solution of the system A*X = B, D*X = 0.
!
!       SIGMA is an output array of length N which contains the
!         diagonal elements of the upper triangular matrix S.
!
!       WA is a work array of length N.
!
!***SEE ALSO  DNLS1, DNLS1E
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DQRSLV
  INTEGER N,LDR
  INTEGER IPVT(*)
  DOUBLE PRECISION R(LDR,*),DIAG(*),QTB(*),X(*),SIGMA(*),WA(*)
  INTEGER I,J,JP1,K,KP1,L,NSING
  DOUBLE PRECISION COS,COTAN,P5,P25,QTBPJ,SIN,SUM,TAN,TEMP,ZERO
  SAVE P5, P25, ZERO
  DATA P5,P25,ZERO /5.0D-1,2.5D-1,0.0D0/
!***FIRST EXECUTABLE STATEMENT  DQRSLV
  DO 20 J = 1, N
     DO 10 I = J, N
        R(I,J) = R(J,I)
   10       CONTINUE
     X(J) = R(J,J)
     WA(J) = QTB(J)
   20    CONTINUE
!
!     ELIMINATE THE DIAGONAL MATRIX D USING A GIVENS ROTATION.
!
  DO 100 J = 1, N
!
!        PREPARE THE ROW OF D TO BE ELIMINATED, LOCATING THE
!        DIAGONAL ELEMENT USING P FROM THE QR FACTORIZATION.
!
     L = IPVT(J)
     if (DIAG(L)  ==  ZERO) go to 90
     DO 30 K = J, N
        SIGMA(K) = ZERO
   30       CONTINUE
     SIGMA(J) = DIAG(L)
!
!        THE TRANSFORMATIONS TO ELIMINATE THE ROW OF D
!        MODIFY ONLY A SINGLE ELEMENT OF (Q TRANSPOSE)*B
!        BEYOND THE FIRST N, WHICH IS INITIALLY ZERO.
!
     QTBPJ = ZERO
     DO 80 K = J, N
!
!           DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
!           APPROPRIATE ELEMENT IN THE CURRENT ROW OF D.
!
        if (SIGMA(K)  ==  ZERO) go to 70
        if (ABS(R(K,K))  >=  ABS(SIGMA(K))) go to 40
           COTAN = R(K,K)/SIGMA(K)
           SIN = P5/SQRT(P25+P25*COTAN**2)
           COS = SIN*COTAN
           go to 50
   40       CONTINUE
           TAN = SIGMA(K)/R(K,K)
           COS = P5/SQRT(P25+P25*TAN**2)
           SIN = COS*TAN
   50       CONTINUE
!
!           COMPUTE THE MODIFIED DIAGONAL ELEMENT OF R AND
!           THE MODIFIED ELEMENT OF ((Q TRANSPOSE)*B,0).
!
        R(K,K) = COS*R(K,K) + SIN*SIGMA(K)
        TEMP = COS*WA(K) + SIN*QTBPJ
        QTBPJ = -SIN*WA(K) + COS*QTBPJ
        WA(K) = TEMP
!
!           ACCUMULATE THE TRANSFORMATION IN THE ROW OF S.
!
        KP1 = K + 1
        if (N  <  KP1) go to 70
        DO 60 I = KP1, N
           TEMP = COS*R(I,K) + SIN*SIGMA(I)
           SIGMA(I) = -SIN*R(I,K) + COS*SIGMA(I)
           R(I,K) = TEMP
   60          CONTINUE
   70       CONTINUE
   80       CONTINUE
   90    CONTINUE
!
!        STORE THE DIAGONAL ELEMENT OF S AND RESTORE
!        THE CORRESPONDING DIAGONAL ELEMENT OF R.
!
     SIGMA(J) = R(J,J)
     R(J,J) = X(J)
  100    CONTINUE
!
!     SOLVE THE TRIANGULAR SYSTEM FOR Z. if THE SYSTEM IS
!     SINGULAR, THEN OBTAIN A LEAST SQUARES SOLUTION.
!
  NSING = N
  DO 110 J = 1, N
     if (SIGMA(J)  ==  ZERO .AND. NSING  ==  N) NSING = J - 1
     if (NSING  <  N) WA(J) = ZERO
  110    CONTINUE
  if (NSING  <  1) go to 150
  DO 140 K = 1, NSING
     J = NSING - K + 1
     SUM = ZERO
     JP1 = J + 1
     if (NSING  <  JP1) go to 130
     DO 120 I = JP1, NSING
        SUM = SUM + R(I,J)*WA(I)
  120       CONTINUE
  130    CONTINUE
     WA(J) = (WA(J) - SUM)/SIGMA(J)
  140    CONTINUE
  150 CONTINUE
!
!     PERMUTE THE COMPONENTS OF Z BACK TO COMPONENTS OF X.
!
  DO 160 J = 1, N
     L = IPVT(J)
     X(L) = WA(J)
  160    CONTINUE
  return
!
!     LAST CARD OF SUBROUTINE DQRSLV.
!
end
