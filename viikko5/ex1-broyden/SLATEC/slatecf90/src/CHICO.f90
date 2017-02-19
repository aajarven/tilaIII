subroutine CHICO (A, LDA, N, KPVT, RCOND, Z)
!
!! CHICO factors a complex Hermitian matrix by elimination with symmetric ...
!  pivoting and estimate the condition of the matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2D1A
!***TYPE      COMPLEX (SSICO-S, DSICO-D, CHICO-C, CSICO-C)
!***KEYWORDS  CONDITION NUMBER, HERMITIAN, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CHICO factors a complex Hermitian matrix by elimination with
!     symmetric pivoting and estimates the condition of the matrix.
!
!     If  RCOND  is not needed, CHIFA is slightly faster.
!     To solve  A*X = B , follow CHICO by CHISL.
!     To compute  INVERSE(A)*C , follow CHICO by CHISL.
!     To compute  INVERSE(A) , follow CHICO by CHIDI.
!     To compute  DETERMINANT(A) , follow CHICO by CHIDI.
!     To compute  INERTIA(A), follow CHICO by CHIDI.
!
!     On Entry
!
!        A       COMPLEX(LDA, N)
!                the Hermitian matrix to be factored.
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
!                The factorization can be written  A = U*D*CTRANS(U)
!                where  U  is a product of permutation and unit
!                upper triangular matrices , CTRANS(U) is the
!                conjugate transpose of  U , and  D  is block diagonal
!                with 1 by 1 and 2 by 2 blocks.
!
!        KVPT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        RCOND   REAL
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
!        Z       COMPLEX(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CDOTC, CHIFA, CSSCAL, SCASUM
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
!***END PROLOGUE  CHICO
  INTEGER LDA,N,KPVT(*)
  COMPLEX A(LDA,*),Z(*)
  REAL RCOND
!
  COMPLEX AK,AKM1,BK,BKM1,CDOTC,DENOM,EK,T
  REAL ANORM,S,SCASUM,YNORM
  INTEGER I,INFO,J,JM1,K,KP,KPS,KS
  COMPLEX ZDUM,ZDUM2,CSIGN1
  REAL CABS1
  CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
  CSIGN1(ZDUM,ZDUM2) = CABS1(ZDUM)*(ZDUM2/CABS1(ZDUM2))
!
!     FIND NORM OF A USING ONLY UPPER HALF
!
!***FIRST EXECUTABLE STATEMENT  CHICO
  DO 30 J = 1, N
     Z(J) = CMPLX(SCASUM(J,A(1,J),1),0.0E0)
     JM1 = J - 1
     if (JM1  <  1) go to 20
     DO 10 I = 1, JM1
        Z(I) = CMPLX(REAL(Z(I))+CABS1(A(I,J)),0.0E0)
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
  ANORM = 0.0E0
  DO 40 J = 1, N
     ANORM = MAX(ANORM,REAL(Z(J)))
   40 CONTINUE
!
!     FACTOR
!
  call CHIFA(A,LDA,N,KPVT,INFO)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!     SOLVE U*D*W = E
!
  EK = (1.0E0,0.0E0)
  DO 50 J = 1, N
     Z(J) = (0.0E0,0.0E0)
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
     if (CABS1(Z(K))  /=  0.0E0) EK = CSIGN1(EK,Z(K))
     Z(K) = Z(K) + EK
     call CAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
     if (KS  ==  1) go to 80
        if (CABS1(Z(K-1))  /=  0.0E0) EK = CSIGN1(EK,Z(K-1))
        Z(K-1) = Z(K-1) + EK
        call CAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
   80    CONTINUE
     if (KS  ==  2) go to 100
        if (CABS1(Z(K))  <=  CABS1(A(K,K))) go to 90
           S = CABS1(A(K,K))/CABS1(Z(K))
           call CSSCAL(N,S,Z,1)
           EK = CMPLX(S,0.0E0)*EK
   90       CONTINUE
        if (CABS1(A(K,K))  /=  0.0E0) Z(K) = Z(K)/A(K,K)
        if (CABS1(A(K,K))  ==  0.0E0) Z(K) = (1.0E0,0.0E0)
     go to 110
  100    CONTINUE
        AK = A(K,K)/CONJG(A(K-1,K))
        AKM1 = A(K-1,K-1)/A(K-1,K)
        BK = Z(K)/CONJG(A(K-1,K))
        BKM1 = Z(K-1)/A(K-1,K)
        DENOM = AK*AKM1 - 1.0E0
        Z(K) = (AKM1*BK - BKM1)/DENOM
        Z(K-1) = (AK*BKM1 - BK)/DENOM
  110    CONTINUE
     K = K - KS
  go to 60
  120 CONTINUE
  S = 1.0E0/SCASUM(N,Z,1)
  call CSSCAL(N,S,Z,1)
!
!     SOLVE CTRANS(U)*Y = W
!
  K = 1
  130 if (K  >  N) go to 160
     KS = 1
     if (KPVT(K)  <  0) KS = 2
     if (K  ==  1) go to 150
        Z(K) = Z(K) + CDOTC(K-1,A(1,K),1,Z(1),1)
        if (KS  ==  2) &
           Z(K+1) = Z(K+1) + CDOTC(K-1,A(1,K+1),1,Z(1),1)
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
  S = 1.0E0/SCASUM(N,Z,1)
  call CSSCAL(N,S,Z,1)
!
  YNORM = 1.0E0
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
        call CAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
        if (KS  ==  2) call CAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
  190    CONTINUE
     if (KS  ==  2) go to 210
        if (CABS1(Z(K))  <=  CABS1(A(K,K))) go to 200
           S = CABS1(A(K,K))/CABS1(Z(K))
           call CSSCAL(N,S,Z,1)
           YNORM = S*YNORM
  200       CONTINUE
        if (CABS1(A(K,K))  /=  0.0E0) Z(K) = Z(K)/A(K,K)
        if (CABS1(A(K,K))  ==  0.0E0) Z(K) = (1.0E0,0.0E0)
     go to 220
  210    CONTINUE
        AK = A(K,K)/CONJG(A(K-1,K))
        AKM1 = A(K-1,K-1)/A(K-1,K)
        BK = Z(K)/CONJG(A(K-1,K))
        BKM1 = Z(K-1)/A(K-1,K)
        DENOM = AK*AKM1 - 1.0E0
        Z(K) = (AKM1*BK - BKM1)/DENOM
        Z(K-1) = (AK*BKM1 - BK)/DENOM
  220    CONTINUE
     K = K - KS
  go to 170
  230 CONTINUE
  S = 1.0E0/SCASUM(N,Z,1)
  call CSSCAL(N,S,Z,1)
  YNORM = S*YNORM
!
!     SOLVE CTRANS(U)*Z = V
!
  K = 1
  240 if (K  >  N) go to 270
     KS = 1
     if (KPVT(K)  <  0) KS = 2
     if (K  ==  1) go to 260
        Z(K) = Z(K) + CDOTC(K-1,A(1,K),1,Z(1),1)
        if (KS  ==  2) &
           Z(K+1) = Z(K+1) + CDOTC(K-1,A(1,K+1),1,Z(1),1)
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
  S = 1.0E0/SCASUM(N,Z,1)
  call CSSCAL(N,S,Z,1)
  YNORM = S*YNORM
!
  if (ANORM  /=  0.0E0) RCOND = YNORM/ANORM
  if (ANORM  ==  0.0E0) RCOND = 0.0E0
  return
end
