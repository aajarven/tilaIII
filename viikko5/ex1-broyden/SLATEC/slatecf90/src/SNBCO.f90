subroutine SNBCO (ABE, LDA, N, ML, MU, IPVT, RCOND, Z)
!
!! SNBCO factors a band matrix using Gaussian elimination and estimates ...
!  the condition number.
!
!***LIBRARY   SLATEC
!***CATEGORY  D2A2
!***TYPE      SINGLE PRECISION (SNBCO-S, DNBCO-D, CNBCO-C)
!***KEYWORDS  BANDED, LINEAR EQUATIONS, MATRIX FACTORIZATION,
!             NONSYMMETRIC
!***AUTHOR  Voorhees, E. A., (LANL)
!***DESCRIPTION
!
!     SNBCO factors a real band matrix by Gaussian
!     elimination and estimates the condition of the matrix.
!
!     If RCOND is not needed, SNBFA is slightly faster.
!     To solve  A*X = B , follow SNBCO by SNBSL.
!     To compute  INVERSE(A)*C , follow SNBCO by SNBSL.
!     To compute  DETERMINANT(A) , follow SNBCO by SNBDI.
!
!     On Entry
!
!        ABE     REAL(LDA, NC)
!                contains the matrix in band storage.  The rows
!                of the original matrix are stored in the rows
!                of ABE and the diagonals of the original matrix
!                are stored in columns 1 through ML+MU+1 of ABE.
!                NC must be  >=  2*ML+MU+1 .
!                See the comments below for details.
!
!        LDA     INTEGER
!                the leading dimension of the array ABE.
!                LDA must be  >=  N .
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
!                More efficient if ML  <=  MU .
!
!     On Return
!
!        ABE     an upper triangular matrix in band storage
!                and the multipliers which were used to obtain it.
!                The factorization can be written  A = L*U , where
!                L is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        RCOND   REAL
!                an estimate of the reciprocal condition of  A .
!                For the system  A*X = B , relative perturbations
!                in  A  and  B  of size  EPSILON  may cause
!                relative perturbations in  X  of size  EPSILON/RCOND .
!                If  RCOND  is so small that the logical expression
!                         1.0 + RCOND  ==  1.0
!                is true, then  A  may be singular to working
!                precision.  In particular,  RCOND  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        Z       REAL(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!     Band Storage
!
!           If  A  is a band matrix, the following program segment
!           will set up the input.
!
!                   ML = (band width below the diagonal)
!                   MU = (band width above the diagonal)
!                   DO 20 I = 1, N
!                      J1 = MAX(1, I-ML)
!                      J2 = MIN(N, I+MU)
!                      DO 10 J = J1, J2
!                         K = J - I + ML + 1
!                         ABE(I,K) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           This uses columns  1  through  ML+MU+1  of ABE .
!           Furthermore,  ML  additional columns are needed in
!           ABE  starting with column  ML+MU+2  for elements
!           generated during the triangularization.  The total
!           number of columns needed in  ABE  is  2*ML+MU+1 .
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
!      then  N = 6, ML = 1, MU = 2, LDA  >=  5  and ABE should contain
!
!            * 111213  +     , * = not used
!           21222324  +     , + = used for pivoting
!           32333435  +
!           43444546  +
!           545556  *  +
!           6566  *  *  +
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  SASUM, SAXPY, SDOT, SNBFA, SSCAL
!***REVISION HISTORY  (YYMMDD)
!   800723  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SNBCO
  INTEGER LDA,N,ML,MU,IPVT(*)
  REAL ABE(LDA,*),Z(*)
  REAL RCOND
!
  REAL SDOT,EK,T,WK,WKM
  REAL ANORM,S,SASUM,SM,YNORM
  INTEGER I,INFO,J,JU,K,KB,KP1,L,LDB,LM,LZ,M,ML1,MM,NL,NU
!***FIRST EXECUTABLE STATEMENT  SNBCO
  ML1=ML+1
  LDB = LDA - 1
  ANORM = 0.0E0
  DO 10 J = 1, N
    NU = MIN(MU,J-1)
    NL = MIN(ML,N-J)
    L = 1 + NU + NL
    ANORM = MAX(ANORM,SASUM(L,ABE(J+NL,ML1-NL),LDB))
   10 CONTINUE
!
!     FACTOR
!
  call SNBFA(ABE,LDA,N,ML,MU,IPVT,INFO)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF  W WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.
!
!     SOLVE TRANS(U)*W = E
!
  EK = 1.0E0
  DO 20 J = 1, N
    Z(J) = 0.0E0
   20 CONTINUE
  M = ML + MU + 1
  JU = 0
  DO 100 K = 1, N
    if (Z(K)  /=  0.0E0) EK = SIGN(EK,-Z(K))
    if (ABS(EK-Z(K))  <=  ABS(ABE(K,ML1))) go to 30
      S = ABS(ABE(K,ML1))/ABS(EK-Z(K))
      call SSCAL(N,S,Z,1)
      EK = S*EK
   30   CONTINUE
    WK = EK - Z(K)
    WKM = -EK - Z(K)
    S = ABS(WK)
    SM = ABS(WKM)
    if (ABE(K,ML1)  ==  0.0E0) go to 40
      WK = WK/ABE(K,ML1)
      WKM = WKM/ABE(K,ML1)
    go to 50
   40   CONTINUE
      WK = 1.0E0
      WKM = 1.0E0
   50   CONTINUE
    KP1 = K + 1
    JU = MIN(MAX(JU,MU+IPVT(K)),N)
    MM = ML1
    if (KP1  >  JU) go to 90
      DO 60 I = KP1, JU
        MM = MM + 1
        SM = SM + ABS(Z(I)+WKM*ABE(K,MM))
        Z(I) = Z(I) + WK*ABE(K,MM)
        S = S + ABS(Z(I))
   60     CONTINUE
      if (S  >=  SM) go to 80
        T = WKM -WK
        WK = WKM
        MM = ML1
        DO 70 I = KP1, JU
          MM = MM + 1
          Z(I) = Z(I) + T*ABE(K,MM)
   70       CONTINUE
   80     CONTINUE
   90   CONTINUE
  Z(K) = WK
  100 CONTINUE
  S = 1.0E0/SASUM(N,Z,1)
  call SSCAL(N,S,Z,1)
!
!     SOLVE TRANS(L)*Y = W
!
  DO 120 KB = 1, N
    K = N + 1 - KB
    NL = MIN(ML,N-K)
    if (K  <  N) Z(K) = Z(K) + SDOT(NL,ABE(K+NL,ML1-NL),-LDB,Z(K+1) &
    ,1)
    if (ABS(Z(K))  <=  1.0E0) go to 110
      S = 1.0E0/ABS(Z(K))
      call SSCAL(N,S,Z,1)
  110   CONTINUE
    L = IPVT(K)
    T = Z(L)
    Z(L) = Z(K)
    Z(K) = T
  120 CONTINUE
  S = 1.0E0/SASUM(N,Z,1)
  call SSCAL(N,S,Z,1)
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
    NL = MIN(ML,N-K)
    if (K  <  N) call SAXPY(NL,T,ABE(K+NL,ML1-NL),-LDB,Z(K+1),1)
    if (ABS(Z(K))  <=  1.0E0) go to 130
      S = 1.0E0/ABS(Z(K))
      call SSCAL(N,S,Z,1)
      YNORM = S*YNORM
  130   CONTINUE
  140 CONTINUE
  S = 1.0E0/SASUM(N,Z,1)
  call SSCAL(N,S,Z,1)
  YNORM = S*YNORM
!
!     SOLVE  U*Z = V
!
  DO 160 KB = 1, N
    K = N + 1 - KB
    if (ABS(Z(K))  <=  ABS(ABE(K,ML1))) go to 150
      S = ABS(ABE(K,ML1))/ABS(Z(K))
      call SSCAL(N,S,Z,1)
      YNORM = S*YNORM
  150   CONTINUE
    if (ABE(K,ML1)  /=  0.0E0) Z(K) = Z(K)/ABE(K,ML1)
    if (ABE(K,ML1)  ==  0.0E0) Z(K) = 1.0E0
    LM = MIN(K,M) - 1
    LZ = K - LM
    T = -Z(K)
    call SAXPY(LM,T,ABE(K-1,ML+2),-LDB,Z(LZ),1)
  160 CONTINUE
!     MAKE ZNORM = 1.0E0
  S = 1.0E0/SASUM(N,Z,1)
  call SSCAL(N,S,Z,1)
  YNORM = S*YNORM
!
  if (ANORM  /=  0.0E0) RCOND = YNORM/ANORM
  if (ANORM  ==  0.0E0) RCOND = 0.0E0
  return
end
