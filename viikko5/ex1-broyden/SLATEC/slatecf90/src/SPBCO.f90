subroutine SPBCO (ABD, LDA, N, M, RCOND, Z, INFO)
!
!! SPBCO factors a real symmetric positive definite matrix stored in ...
!            band form and estimate the condition number of the matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B2
!***TYPE      SINGLE PRECISION (SPBCO-S, DPBCO-D, CPBCO-C)
!***KEYWORDS  BANDED, CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION, POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     SPBCO factors a real symmetric positive definite matrix
!     stored in band form and estimates the condition of the matrix.
!
!     If  RCOND  is not needed, SPBFA is slightly faster.
!     To solve  A*X = B , follow SPBCO by SPBSL.
!     To compute  INVERSE(A)*C , follow SPBCO by SPBSL.
!     To compute  DETERMINANT(A) , follow SPBCO by SPBDI.
!
!     On Entry
!
!        ABD     REAL(LDA, N)
!                the matrix to be factored.  The columns of the upper
!                triangle are stored in the columns of ABD and the
!                diagonals of the upper triangle are stored in the
!                rows of ABD .  See the comments below for details.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!                LDA must be  >=  M + 1 .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        M       INTEGER
!                the number of diagonals above the main diagonal.
!                0  <=  M  <  N .
!
!     On Return
!
!        ABD     an upper triangular matrix  R , stored in band
!                form, so that  A = TRANS(R)*R .
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
!        Z       REAL(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is singular to working precision, then  Z  is
!                an approximate null vector in the sense that
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!                If  INFO  /=  0 , Z  is unchanged.
!
!        INFO    INTEGER
!                = 0  for normal return.
!                = K  signals an error condition.  The leading minor
!                     of order  K  is not positive definite.
!
!     Band Storage
!
!           If  A  is a symmetric positive definite band matrix,
!           the following program segment will set up the input.
!
!                   M = (band width above diagonal)
!                   DO 20 J = 1, N
!                      I1 = MAX(1, J-M)
!                      DO 10 I = I1, J
!                         K = I-J+M+1
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           This uses  M + 1  rows of  A , except for the  M by M
!           upper left triangle, which is ignored.
!
!     Example:  If the original matrix is
!
!           111213  0  0  0
!           12222324  0  0
!           1323333435  0
!            02434444546
!            0  035455556
!            0  0  0465666
!
!     then  N = 6 , M = 2  and  ABD  should contain
!
!            *  * 13243546
!            * 1223344556
!           112233445566
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  SASUM, SAXPY, SDOT, SPBFA, SSCAL
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SPBCO
  INTEGER LDA,N,M,INFO
  REAL ABD(LDA,*),Z(*)
  REAL RCOND
!
  REAL SDOT,EK,T,WK,WKM
  REAL ANORM,S,SASUM,SM,YNORM
  INTEGER I,J,J2,K,KB,KP1,L,LA,LB,LM,MU
!
!     FIND NORM OF A
!
!***FIRST EXECUTABLE STATEMENT  SPBCO
  DO 30 J = 1, N
     L = MIN(J,M+1)
     MU = MAX(M+2-J,1)
     Z(J) = SASUM(L,ABD(MU,J),1)
     K = J - L
     if (M  <  MU) go to 20
     DO 10 I = MU, M
        K = K + 1
        Z(K) = Z(K) + ABS(ABD(I,J))
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
  ANORM = 0.0E0
  DO 40 J = 1, N
     ANORM = MAX(ANORM,Z(J))
   40 CONTINUE
!
!     FACTOR
!
  call SPBFA(ABD,LDA,N,M,INFO)
  if (INFO  /=  0) go to 180
!
!        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
!        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .
!        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!        SOLVE TRANS(R)*W = E
!
     EK = 1.0E0
     DO 50 J = 1, N
        Z(J) = 0.0E0
   50    CONTINUE
     DO 110 K = 1, N
        if (Z(K)  /=  0.0E0) EK = SIGN(EK,-Z(K))
        if (ABS(EK-Z(K))  <=  ABD(M+1,K)) go to 60
           S = ABD(M+1,K)/ABS(EK-Z(K))
           call SSCAL(N,S,Z,1)
           EK = S*EK
   60       CONTINUE
        WK = EK - Z(K)
        WKM = -EK - Z(K)
        S = ABS(WK)
        SM = ABS(WKM)
        WK = WK/ABD(M+1,K)
        WKM = WKM/ABD(M+1,K)
        KP1 = K + 1
        J2 = MIN(K+M,N)
        I = M + 1
        if (KP1  >  J2) go to 100
           DO 70 J = KP1, J2
              I = I - 1
              SM = SM + ABS(Z(J)+WKM*ABD(I,J))
              Z(J) = Z(J) + WK*ABD(I,J)
              S = S + ABS(Z(J))
   70          CONTINUE
           if (S  >=  SM) go to 90
              T = WKM - WK
              WK = WKM
              I = M + 1
              DO 80 J = KP1, J2
                 I = I - 1
                 Z(J) = Z(J) + T*ABD(I,J)
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
        Z(K) = WK
  110    CONTINUE
     S = 1.0E0/SASUM(N,Z,1)
     call SSCAL(N,S,Z,1)
!
!        SOLVE  R*Y = W
!
     DO 130 KB = 1, N
        K = N + 1 - KB
        if (ABS(Z(K))  <=  ABD(M+1,K)) go to 120
           S = ABD(M+1,K)/ABS(Z(K))
           call SSCAL(N,S,Z,1)
  120       CONTINUE
        Z(K) = Z(K)/ABD(M+1,K)
        LM = MIN(K-1,M)
        LA = M + 1 - LM
        LB = K - LM
        T = -Z(K)
        call SAXPY(LM,T,ABD(LA,K),1,Z(LB),1)
  130    CONTINUE
     S = 1.0E0/SASUM(N,Z,1)
     call SSCAL(N,S,Z,1)
!
     YNORM = 1.0E0
!
!        SOLVE TRANS(R)*V = Y
!
     DO 150 K = 1, N
        LM = MIN(K-1,M)
        LA = M + 1 - LM
        LB = K - LM
        Z(K) = Z(K) - SDOT(LM,ABD(LA,K),1,Z(LB),1)
        if (ABS(Z(K))  <=  ABD(M+1,K)) go to 140
           S = ABD(M+1,K)/ABS(Z(K))
           call SSCAL(N,S,Z,1)
           YNORM = S*YNORM
  140       CONTINUE
        Z(K) = Z(K)/ABD(M+1,K)
  150    CONTINUE
     S = 1.0E0/SASUM(N,Z,1)
     call SSCAL(N,S,Z,1)
     YNORM = S*YNORM
!
!        SOLVE  R*Z = W
!
     DO 170 KB = 1, N
        K = N + 1 - KB
        if (ABS(Z(K))  <=  ABD(M+1,K)) go to 160
           S = ABD(M+1,K)/ABS(Z(K))
           call SSCAL(N,S,Z,1)
           YNORM = S*YNORM
  160       CONTINUE
        Z(K) = Z(K)/ABD(M+1,K)
        LM = MIN(K-1,M)
        LA = M + 1 - LM
        LB = K - LM
        T = -Z(K)
        call SAXPY(LM,T,ABD(LA,K),1,Z(LB),1)
  170    CONTINUE
!        MAKE ZNORM = 1.0
     S = 1.0E0/SASUM(N,Z,1)
     call SSCAL(N,S,Z,1)
     YNORM = S*YNORM
!
     if (ANORM  /=  0.0E0) RCOND = YNORM/ANORM
     if (ANORM  ==  0.0E0) RCOND = 0.0E0
  180 CONTINUE
  return
end
