subroutine CPOCO (A, LDA, N, RCOND, Z, INFO)
!
!! CPOCO factors a complex Hermitian positive definite matrix ...
!  and estimates the condition number of the matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2D1B
!***TYPE      COMPLEX (SPOCO-S, DPOCO-D, CPOCO-C)
!***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION, POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CPOCO factors a complex Hermitian positive definite matrix
!     and estimates the condition of the matrix.
!
!     If  RCOND  is not needed, CPOFA is slightly faster.
!     To solve  A*X = B , follow CPOCO by CPOSL.
!     To compute  INVERSE(A)*C , follow CPOCO by CPOSL.
!     To compute  DETERMINANT(A) , follow CPOCO by CPODI.
!     To compute  INVERSE(A) , follow CPOCO by CPODI.
!
!     On Entry
!
!        A       COMPLEX(LDA, N)
!                the Hermitian matrix to be factored.  Only the
!                diagonal and upper triangle are used.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       an upper triangular matrix  R  so that  A =
!                CTRANS(R)*R where  CTRANS(R)  is the conjugate
!                transpose.  The strict lower triangle is unaltered.
!                If  INFO  /=  0 , the factorization is not complete.
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
!                underflows.  If INFO  /=  0 , RCOND is unchanged.
!
!        Z       COMPLEX(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!                If  INFO  /=  0 , Z  is unchanged.
!
!        INFO    INTEGER
!                = 0  for normal return.
!                = K  signals an error condition.  The leading minor
!                     of order  K  is not positive definite.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CDOTC, CPOFA, CSSCAL, SCASUM
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CPOCO
  INTEGER LDA,N,INFO
  COMPLEX A(LDA,*),Z(*)
  REAL RCOND
!
  COMPLEX CDOTC,EK,T,WK,WKM
  REAL ANORM,S,SCASUM,SM,YNORM
  INTEGER I,J,JM1,K,KB,KP1
  COMPLEX ZDUM,ZDUM2,CSIGN1
  REAL CABS1
  CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
  CSIGN1(ZDUM,ZDUM2) = CABS1(ZDUM)*(ZDUM2/CABS1(ZDUM2))
!
!     FIND NORM OF A USING ONLY UPPER HALF
!
!***FIRST EXECUTABLE STATEMENT  CPOCO
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
  call CPOFA(A,LDA,N,INFO)
  if (INFO  /=  0) go to 180
!
!        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
!        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!        GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(R)*W = E .
!        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!        SOLVE CTRANS(R)*W = E
!
     EK = (1.0E0,0.0E0)
     DO 50 J = 1, N
        Z(J) = (0.0E0,0.0E0)
   50    CONTINUE
     DO 110 K = 1, N
        if (CABS1(Z(K))  /=  0.0E0) EK = CSIGN1(EK,-Z(K))
        if (CABS1(EK-Z(K))  <=  REAL(A(K,K))) go to 60
           S = REAL(A(K,K))/CABS1(EK-Z(K))
           call CSSCAL(N,S,Z,1)
           EK = CMPLX(S,0.0E0)*EK
   60       CONTINUE
        WK = EK - Z(K)
        WKM = -EK - Z(K)
        S = CABS1(WK)
        SM = CABS1(WKM)
        WK = WK/A(K,K)
        WKM = WKM/A(K,K)
        KP1 = K + 1
        if (KP1  >  N) go to 100
           DO 70 J = KP1, N
              SM = SM + CABS1(Z(J)+WKM*CONJG(A(K,J)))
              Z(J) = Z(J) + WK*CONJG(A(K,J))
              S = S + CABS1(Z(J))
   70          CONTINUE
           if (S  >=  SM) go to 90
              T = WKM - WK
              WK = WKM
              DO 80 J = KP1, N
                 Z(J) = Z(J) + T*CONJG(A(K,J))
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
        Z(K) = WK
  110    CONTINUE
     S = 1.0E0/SCASUM(N,Z,1)
     call CSSCAL(N,S,Z,1)
!
!        SOLVE R*Y = W
!
     DO 130 KB = 1, N
        K = N + 1 - KB
        if (CABS1(Z(K))  <=  REAL(A(K,K))) go to 120
           S = REAL(A(K,K))/CABS1(Z(K))
           call CSSCAL(N,S,Z,1)
  120       CONTINUE
        Z(K) = Z(K)/A(K,K)
        T = -Z(K)
        call CAXPY(K-1,T,A(1,K),1,Z(1),1)
  130    CONTINUE
     S = 1.0E0/SCASUM(N,Z,1)
     call CSSCAL(N,S,Z,1)
!
     YNORM = 1.0E0
!
!        SOLVE CTRANS(R)*V = Y
!
     DO 150 K = 1, N
        Z(K) = Z(K) - CDOTC(K-1,A(1,K),1,Z(1),1)
        if (CABS1(Z(K))  <=  REAL(A(K,K))) go to 140
           S = REAL(A(K,K))/CABS1(Z(K))
           call CSSCAL(N,S,Z,1)
           YNORM = S*YNORM
  140       CONTINUE
        Z(K) = Z(K)/A(K,K)
  150    CONTINUE
     S = 1.0E0/SCASUM(N,Z,1)
     call CSSCAL(N,S,Z,1)
     YNORM = S*YNORM
!
!        SOLVE R*Z = V
!
     DO 170 KB = 1, N
        K = N + 1 - KB
        if (CABS1(Z(K))  <=  REAL(A(K,K))) go to 160
           S = REAL(A(K,K))/CABS1(Z(K))
           call CSSCAL(N,S,Z,1)
           YNORM = S*YNORM
  160       CONTINUE
        Z(K) = Z(K)/A(K,K)
        T = -Z(K)
        call CAXPY(K-1,T,A(1,K),1,Z(1),1)
  170    CONTINUE
!        MAKE ZNORM = 1.0
     S = 1.0E0/SCASUM(N,Z,1)
     call CSSCAL(N,S,Z,1)
     YNORM = S*YNORM
!
     if (ANORM  /=  0.0E0) RCOND = YNORM/ANORM
     if (ANORM  ==  0.0E0) RCOND = 0.0E0
  180 CONTINUE
  return
end
