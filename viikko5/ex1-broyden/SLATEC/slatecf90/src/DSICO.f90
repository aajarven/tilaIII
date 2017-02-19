subroutine DSICO (A, LDA, N, KPVT, RCOND, Z)
!
!! DSICO factors a symmetric matrix by elimination with symmetric ...
!            pivoting and estimate the condition number of the matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1A
!***TYPE      DOUBLE PRECISION (SSICO-S, DSICO-D, CHICO-C, CSICO-C)
!***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION, SYMMETRIC
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DSICO factors a double precision symmetric matrix by elimination
!     with symmetric pivoting and estimates the condition of the
!     matrix.
!
!     If  RCOND  is not needed, DSIFA is slightly faster.
!     To solve  A*X = B , follow DSICO by DSISL.
!     To compute  INVERSE(A)*C , follow DSICO by DSISL.
!     To compute  INVERSE(A) , follow DSICO by DSIDI.
!     To compute  DETERMINANT(A) , follow DSICO by DSIDI.
!     To compute  INERTIA(A), follow DSICO by DSIDI.
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the symmetric matrix to be factored.
!                Only the diagonal and upper triangle are used.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     Output
!
!        A       a block diagonal matrix and the multipliers which
!                were used to obtain it.
!                The factorization can be written  A = U*D*TRANS(U)
!                where  U  is a product of permutation and unit
!                upper triangular matrices, TRANS(U) is the
!                transpose of  U , and  D  is block diagonal
!                with 1 by 1 and 2 by 2 blocks.
!
!        KPVT    INTEGER(N)
!                an integer vector of pivot indices.
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
!***ROUTINES CALLED  DASUM, DAXPY, DDOT, DSCAL, DSIFA
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891107  Modified routine equivalence list.  (WRB)
!   891107  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DSICO
  INTEGER LDA,N,KPVT(*)
  DOUBLE PRECISION A(LDA,*),Z(*)
  DOUBLE PRECISION RCOND
!
  DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,EK,T
  DOUBLE PRECISION ANORM,S,DASUM,YNORM
  INTEGER I,INFO,J,JM1,K,KP,KPS,KS
!
!     FIND NORM OF A USING ONLY UPPER HALF
!
!***FIRST EXECUTABLE STATEMENT  DSICO
  DO J = 1, N
     Z(J) = DASUM(J,A(1,J),1)
     JM1 = J - 1
     DO I = 1, JM1
        Z(I) = Z(I) + ABS(A(I,J))
     end do
  end do

  ANORM = 0.0D0
  DO 40 J = 1, N
     ANORM = MAX(ANORM,Z(J))
   40 CONTINUE
!
!     FACTOR
!
  call DSIFA(A,LDA,N,KPVT,INFO)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!     SOLVE U*D*W = E
!
  EK = 1.0D0
  DO 50 J = 1, N
     Z(J) = 0.0D0
   50 CONTINUE
  K = N
   60 if (K  ==  0) go to 120
     KS = 1
     if (KPVT(K)  <  0) KS = 2
     KP = ABS(KPVT(K))
     KPS = K + 1 - KS
     if (KP  ==  KPS) go to 70
        T = Z(KPS)
        Z(KPS) = Z(KP)
        Z(KP) = T
   70    CONTINUE
     if (Z(K)  /=  0.0D0) EK = SIGN(EK,Z(K))
     Z(K) = Z(K) + EK
     call DAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
     if (KS  ==  1) go to 80
        if (Z(K-1)  /=  0.0D0) EK = SIGN(EK,Z(K-1))
        Z(K-1) = Z(K-1) + EK
        call DAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
   80    CONTINUE
     if (KS  ==  2) go to 100
        if (ABS(Z(K))  <=  ABS(A(K,K))) go to 90
           S = ABS(A(K,K))/ABS(Z(K))
           call DSCAL(N,S,Z,1)
           EK = S*EK
   90       CONTINUE
        if (A(K,K)  /=  0.0D0) Z(K) = Z(K)/A(K,K)
        if (A(K,K)  ==  0.0D0) Z(K) = 1.0D0
     go to 110
  100    CONTINUE
        AK = A(K,K)/A(K-1,K)
        AKM1 = A(K-1,K-1)/A(K-1,K)
        BK = Z(K)/A(K-1,K)
        BKM1 = Z(K-1)/A(K-1,K)
        DENOM = AK*AKM1 - 1.0D0
        Z(K) = (AKM1*BK - BKM1)/DENOM
        Z(K-1) = (AK*BKM1 - BK)/DENOM
  110    CONTINUE
     K = K - KS
  go to 60
  120 CONTINUE
  S = 1.0D0/DASUM(N,Z,1)
  call DSCAL(N,S,Z,1)
!
!     SOLVE TRANS(U)*Y = W
!
  K = 1
  130 if (K  >  N) go to 160
     KS = 1
     if (KPVT(K)  <  0) KS = 2
     if (K  ==  1) go to 150
        Z(K) = Z(K) + DDOT(K-1,A(1,K),1,Z(1),1)
        if (KS  ==  2) &
           Z(K+1) = Z(K+1) + DDOT(K-1,A(1,K+1),1,Z(1),1)
        KP = ABS(KPVT(K))
        if (KP  ==  K) go to 140
           T = Z(K)
           Z(K) = Z(KP)
           Z(KP) = T
  140       CONTINUE
  150    CONTINUE
     K = K + KS
  go to 130
  160 CONTINUE
  S = 1.0D0/DASUM(N,Z,1)
  call DSCAL(N,S,Z,1)
!
  YNORM = 1.0D0
!
!     SOLVE U*D*V = Y
!
  K = N
  170 if (K  ==  0) go to 230
     KS = 1
     if (KPVT(K)  <  0) KS = 2
     if (K  ==  KS) go to 190
        KP = ABS(KPVT(K))
        KPS = K + 1 - KS
        if (KP  ==  KPS) go to 180
           T = Z(KPS)
           Z(KPS) = Z(KP)
           Z(KP) = T
  180       CONTINUE
        call DAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
        if (KS  ==  2) call DAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
  190    CONTINUE
     if (KS  ==  2) go to 210
        if (ABS(Z(K))  <=  ABS(A(K,K))) go to 200
           S = ABS(A(K,K))/ABS(Z(K))
           call DSCAL(N,S,Z,1)
           YNORM = S*YNORM
  200       CONTINUE
        if (A(K,K)  /=  0.0D0) Z(K) = Z(K)/A(K,K)
        if (A(K,K)  ==  0.0D0) Z(K) = 1.0D0
     go to 220
  210    CONTINUE
        AK = A(K,K)/A(K-1,K)
        AKM1 = A(K-1,K-1)/A(K-1,K)
        BK = Z(K)/A(K-1,K)
        BKM1 = Z(K-1)/A(K-1,K)
        DENOM = AK*AKM1 - 1.0D0
        Z(K) = (AKM1*BK - BKM1)/DENOM
        Z(K-1) = (AK*BKM1 - BK)/DENOM
  220    CONTINUE
     K = K - KS
  go to 170
  230 CONTINUE
  S = 1.0D0/DASUM(N,Z,1)
  call DSCAL(N,S,Z,1)
  YNORM = S*YNORM
!
!     SOLVE TRANS(U)*Z = V
!
  K = 1
  240 if (K  >  N) go to 270
     KS = 1
     if (KPVT(K)  <  0) KS = 2
     if (K  ==  1) go to 260
        Z(K) = Z(K) + DDOT(K-1,A(1,K),1,Z(1),1)
        if (KS  ==  2) &
           Z(K+1) = Z(K+1) + DDOT(K-1,A(1,K+1),1,Z(1),1)
        KP = ABS(KPVT(K))
        if (KP  ==  K) go to 250
           T = Z(K)
           Z(K) = Z(KP)
           Z(KP) = T
  250       CONTINUE
  260    CONTINUE
     K = K + KS
  go to 240
  270 CONTINUE
!     MAKE ZNORM = 1.0
  S = 1.0D0/DASUM(N,Z,1)
  call DSCAL(N,S,Z,1)
  YNORM = S*YNORM
!
  if (ANORM  /=  0.0D0) RCOND = YNORM/ANORM
  if (ANORM  ==  0.0D0) RCOND = 0.0D0
  return
end
