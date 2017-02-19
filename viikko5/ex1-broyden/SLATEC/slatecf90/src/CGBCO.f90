subroutine CGBCO (ABD, LDA, N, ML, MU, IPVT, RCOND, Z)
!
!! CGBCO factors a band matrix by Gaussian elimination and ...
!            estimate the condition number of the matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2C2
!***TYPE      COMPLEX (SGBCO-S, DGBCO-D, CGBCO-C)
!***KEYWORDS  BANDED, CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CGBCO factors a complex band matrix by Gaussian
!     elimination and estimates the condition of the matrix.
!
!     If  RCOND  is not needed, CGBFA is slightly faster.
!     To solve  A*X = B , follow CGBCO by CGBSL.
!     To compute  INVERSE(A)*C , follow CGBCO by CGBSL.
!     To compute  DETERMINANT(A) , follow CGBCO by CGBDI.
!
!     On Entry
!
!        ABD     COMPLEX(LDA, N)
!                contains the matrix in band storage.  The columns
!                of the matrix are stored in the columns of  ABD  and
!                the diagonals of the matrix are stored in rows
!                ML+1 through 2*ML+MU+1 of  ABD .
!                See the comments below for details.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!                LDA must be  >=  2*ML + MU + 1 .
!
!        N       INTEGER
!                the order of the original matrix.
!
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!                0  <=  ML  <  N .
!
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!                0  <=  MU  <  N .
!                More efficient if  ML  <=  MU .
!
!     On Return
!
!        ABD     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        RCOND   REAL
!                an estimate of the reciprocal condition of  A .
!                For the system  A*X = B , relative perturbations
!                in  A  And  B  of size  EPSILON  may cause
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
!     Band Storage
!
!           if  A  is a band matrix, the following program segment
!           will set up the input.
!
!                   ML = (band width below the diagonal)
!                   MU = (band width above the diagonal)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX(1, J-MU)
!                      I2 = MIN(N, J+Ml)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
!           In addition, the first  ML  rows in  ABD  are used for
!           elements generated during the triangularization.
!           The total number of rows needed in  ABD  is  2*ML+MU+1 .
!           The  ML+MU by ML+MU  upper left triangle and the
!           ML by ML  lower right triangle are not referenced.
!
!     Example:  If the original matrix is
!
!           111213  0  0  0
!           21222324  0  0
!            032333435  0
!            0  043444546
!            0  0  0545556
!            0  0  0  06566
!
!      then  N = 6, ML = 1, MU = 2, LDA  >=  5  and ABD should contain
!
!            *  *  *  +  +  +  , * = not used
!            *  * 13243546  , + = used for pivoting
!            * 1223344556
!           112233445566
!           2132435465  *
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CDOTC, CGBFA, CSSCAL, SCASUM
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CGBCO
  INTEGER LDA,N,ML,MU,IPVT(*)
  COMPLEX ABD(LDA,*),Z(*)
  REAL RCOND
!
  COMPLEX CDOTC,EK,T,WK,WKM
  REAL ANORM,S,SCASUM,SM,YNORM
  INTEGER IS,INFO,J,JU,K,KB,KP1,L,LA,LM,LZ,M,MM
  COMPLEX ZDUM,ZDUM1,ZDUM2,CSIGN1
  REAL CABS1
  CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
  CSIGN1(ZDUM1,ZDUM2) = CABS1(ZDUM1)*(ZDUM2/CABS1(ZDUM2))
!
!     COMPUTE 1-NORM OF A
!
!***FIRST EXECUTABLE STATEMENT  CGBCO
  ANORM = 0.0E0
  L = ML + 1
  IS = L + MU
  DO 10 J = 1, N
     ANORM = MAX(ANORM,SCASUM(L,ABD(IS,J),1))
     if (IS  >  ML + 1) IS = IS - 1
     if (J  <=  MU) L = L + 1
     if (J  >=  N - ML) L = L - 1
   10 CONTINUE
!
!     FACTOR
!
  call CGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E .
!     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A .
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E .
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!     SOLVE CTRANS(U)*W = E
!
  EK = (1.0E0,0.0E0)
  DO 20 J = 1, N
     Z(J) = (0.0E0,0.0E0)
   20 CONTINUE
  M = ML + MU + 1
  JU = 0
  DO 100 K = 1, N
     if (CABS1(Z(K))  /=  0.0E0) EK = CSIGN1(EK,-Z(K))
     if (CABS1(EK-Z(K))  <=  CABS1(ABD(M,K))) go to 30
        S = CABS1(ABD(M,K))/CABS1(EK-Z(K))
        call CSSCAL(N,S,Z,1)
        EK = CMPLX(S,0.0E0)*EK
   30    CONTINUE
     WK = EK - Z(K)
     WKM = -EK - Z(K)
     S = CABS1(WK)
     SM = CABS1(WKM)
     if (CABS1(ABD(M,K))  ==  0.0E0) go to 40
        WK = WK/CONJG(ABD(M,K))
        WKM = WKM/CONJG(ABD(M,K))
     go to 50
   40    CONTINUE
        WK = (1.0E0,0.0E0)
        WKM = (1.0E0,0.0E0)
   50    CONTINUE
     KP1 = K + 1
     JU = MIN(MAX(JU,MU+IPVT(K)),N)
     MM = M
     if (KP1  >  JU) go to 90
        DO 60 J = KP1, JU
           MM = MM - 1
           SM = SM + CABS1(Z(J)+WKM*CONJG(ABD(MM,J)))
           Z(J) = Z(J) + WK*CONJG(ABD(MM,J))
           S = S + CABS1(Z(J))
   60       CONTINUE
        if (S  >=  SM) go to 80
           T = WKM - WK
           WK = WKM
           MM = M
           DO 70 J = KP1, JU
              MM = MM - 1
              Z(J) = Z(J) + T*CONJG(ABD(MM,J))
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
     Z(K) = WK
  100 CONTINUE
  S = 1.0E0/SCASUM(N,Z,1)
  call CSSCAL(N,S,Z,1)
!
!     SOLVE CTRANS(L)*Y = W
!
  DO 120 KB = 1, N
     K = N + 1 - KB
     LM = MIN(ML,N-K)
     if (K  <  N) Z(K) = Z(K) + CDOTC(LM,ABD(M+1,K),1,Z(K+1),1)
     if (CABS1(Z(K))  <=  1.0E0) go to 110
        S = 1.0E0/CABS1(Z(K))
        call CSSCAL(N,S,Z,1)
  110    CONTINUE
     L = IPVT(K)
     T = Z(L)
     Z(L) = Z(K)
     Z(K) = T
  120 CONTINUE
  S = 1.0E0/SCASUM(N,Z,1)
  call CSSCAL(N,S,Z,1)
!
  YNORM = 1.0E0
!
!     SOLVE L*V = Y
!
  DO 140 K = 1, N
     L = IPVT(K)
     T = Z(L)
     Z(L) = Z(K)
     Z(K) = T
     LM = MIN(ML,N-K)
     if (K  <  N) call CAXPY(LM,T,ABD(M+1,K),1,Z(K+1),1)
     if (CABS1(Z(K))  <=  1.0E0) go to 130
        S = 1.0E0/CABS1(Z(K))
        call CSSCAL(N,S,Z,1)
        YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
  S = 1.0E0/SCASUM(N,Z,1)
  call CSSCAL(N,S,Z,1)
  YNORM = S*YNORM
!
!     SOLVE  U*Z = W
!
  DO 160 KB = 1, N
     K = N + 1 - KB
     if (CABS1(Z(K))  <=  CABS1(ABD(M,K))) go to 150
        S = CABS1(ABD(M,K))/CABS1(Z(K))
        call CSSCAL(N,S,Z,1)
        YNORM = S*YNORM
  150    CONTINUE
     if (CABS1(ABD(M,K))  /=  0.0E0) Z(K) = Z(K)/ABD(M,K)
     if (CABS1(ABD(M,K))  ==  0.0E0) Z(K) = (1.0E0,0.0E0)
     LM = MIN(K,M) - 1
     LA = M - LM
     LZ = K - LM
     T = -Z(K)
     call CAXPY(LM,T,ABD(LA,K),1,Z(LZ),1)
  160 CONTINUE
!     MAKE ZNORM = 1.0
  S = 1.0E0/SCASUM(N,Z,1)
  call CSSCAL(N,S,Z,1)
  YNORM = S*YNORM
!
  if (ANORM  /=  0.0E0) RCOND = YNORM/ANORM
  if (ANORM  ==  0.0E0) RCOND = 0.0E0
  return
end
