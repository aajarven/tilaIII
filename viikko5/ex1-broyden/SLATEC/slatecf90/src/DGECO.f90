subroutine DGECO (A, LDA, N, IPVT, RCOND, Z)
!
!! DGECO factors a matrix using Gaussian elimination and estimate ...
!            the condition number of the matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2A1
!***TYPE      DOUBLE PRECISION (SGECO-S, DGECO-D, CGECO-C)
!***KEYWORDS  CONDITION NUMBER, GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGECO factors a double precision matrix by Gaussian elimination
!     and estimates the condition of the matrix.
!
!     If  RCOND  is not needed, DGEFA is slightly faster.
!     To solve  A*X = B , follow DGECO by DGESL.
!     To compute  INVERSE(A)*C , follow DGECO by DGESL.
!     To compute  DETERMINANT(A) , follow DGECO by DGEDI.
!     To compute  INVERSE(A) , follow DGECO by DGEDI.
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the matrix to be factored.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an INTEGER vector of pivot indices.
!
!        RCOND   DOUBLE PRECISION
!                an estimate of the reciprocal condition of  A .
!                For the system  A*X = B , relative perturbations
!                in  A  and  B  of size  EPSILON  may cause
!                relative perturbations in  X  of size  EPSILON/RCOND .
!                If  RCOND  is so small that the logical expression
!                           1.0 + RCOND  ==  1.0
!                is true, then  A  may be singular to working
!                precision.  In particular,  RCOND  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        Z       DOUBLE PRECISION(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DASUM, DAXPY, DDOT, DGEFA, DSCAL
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGECO
  INTEGER LDA,N,IPVT(*)
  DOUBLE PRECISION A(LDA,*),Z(*)
  DOUBLE PRECISION RCOND
!
  DOUBLE PRECISION DDOT,EK,T,WK,WKM
  DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
  INTEGER INFO,J,K,KB,KP1,L
!
!     COMPUTE 1-NORM OF A
!
!***FIRST EXECUTABLE STATEMENT  DGECO
  ANORM = 0.0D0
  DO 10 J = 1, N
     ANORM = MAX(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
!
!     FACTOR
!
  call DGEFA(A,LDA,N,IPVT,INFO)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.
!
!     SOLVE TRANS(U)*W = E
!
  EK = 1.0D0
  DO 20 J = 1, N
     Z(J) = 0.0D0
   20 CONTINUE
  DO 100 K = 1, N
     if (Z(K)  /=  0.0D0) EK = SIGN(EK,-Z(K))
     if (ABS(EK-Z(K))  <=  ABS(A(K,K))) go to 30
        S = ABS(A(K,K))/ABS(EK-Z(K))
        call DSCAL(N,S,Z,1)
        EK = S*EK
   30    CONTINUE
     WK = EK - Z(K)
     WKM = -EK - Z(K)
     S = ABS(WK)
     SM = ABS(WKM)
     if (A(K,K)  ==  0.0D0) go to 40
        WK = WK/A(K,K)
        WKM = WKM/A(K,K)
     go to 50
   40    CONTINUE
        WK = 1.0D0
        WKM = 1.0D0
   50    CONTINUE
     KP1 = K + 1
     if (KP1  >  N) go to 90
        DO 60 J = KP1, N
           SM = SM + ABS(Z(J)+WKM*A(K,J))
           Z(J) = Z(J) + WK*A(K,J)
           S = S + ABS(Z(J))
   60       CONTINUE
        if (S  >=  SM) go to 80
           T = WKM - WK
           WK = WKM
           DO 70 J = KP1, N
              Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
     Z(K) = WK
  100 CONTINUE
  S = 1.0D0/DASUM(N,Z,1)
  call DSCAL(N,S,Z,1)
!
!     SOLVE TRANS(L)*Y = W
!
  DO 120 KB = 1, N
     K = N + 1 - KB
     if (K  <  N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
     if (ABS(Z(K))  <=  1.0D0) go to 110
        S = 1.0D0/ABS(Z(K))
        call DSCAL(N,S,Z,1)
  110    CONTINUE
     L = IPVT(K)
     T = Z(L)
     Z(L) = Z(K)
     Z(K) = T
  120 CONTINUE
  S = 1.0D0/DASUM(N,Z,1)
  call DSCAL(N,S,Z,1)
!
  YNORM = 1.0D0
!
!     SOLVE L*V = Y
!
  DO 140 K = 1, N
     L = IPVT(K)
     T = Z(L)
     Z(L) = Z(K)
     Z(K) = T
     if (K  <  N) call DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
     if (ABS(Z(K))  <=  1.0D0) go to 130
        S = 1.0D0/ABS(Z(K))
        call DSCAL(N,S,Z,1)
        YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
  S = 1.0D0/DASUM(N,Z,1)
  call DSCAL(N,S,Z,1)
  YNORM = S*YNORM
!
!     SOLVE  U*Z = V
!
  DO 160 KB = 1, N
     K = N + 1 - KB
     if (ABS(Z(K))  <=  ABS(A(K,K))) go to 150
        S = ABS(A(K,K))/ABS(Z(K))
        call DSCAL(N,S,Z,1)
        YNORM = S*YNORM
  150    CONTINUE
     if (A(K,K)  /=  0.0D0) Z(K) = Z(K)/A(K,K)
     if (A(K,K)  ==  0.0D0) Z(K) = 1.0D0
     T = -Z(K)
     call DAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
!     MAKE ZNORM = 1.0
  S = 1.0D0/DASUM(N,Z,1)
  call DSCAL(N,S,Z,1)
  YNORM = S*YNORM
!
  if (ANORM  /=  0.0D0) RCOND = YNORM/ANORM
  if (ANORM  ==  0.0D0) RCOND = 0.0D0
  return
end
