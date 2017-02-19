subroutine CGBFA (ABD, LDA, N, ML, MU, IPVT, INFO)
!
!! CGBFA factors a band matrix using Gaussian elimination.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2C2
!***TYPE      COMPLEX (SGBFA-S, DGBFA-D, CGBFA-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CGBFA factors a complex band matrix by elimination.
!
!     CGBFA is usually called by CGBCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
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
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K)  ==  0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that CGBSL will divide by zero if
!                     called.  Use  RCOND  in CGBCO for a reliable
!                     indication of singularity.
!
!     Band Storage
!
!           If  A  is a band matrix, the following program segment
!           will set up the input.
!
!                   ML = (band width below the diagonal)
!                   MU = (band width above the diagonal)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX(1, J-MU)
!                      I2 = MIN(N, J+ML)
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
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CSCAL, ICAMAX
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CGBFA
  INTEGER LDA,N,ML,MU,IPVT(*),INFO
  COMPLEX ABD(LDA,*)
!
  COMPLEX T
  INTEGER I,ICAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
  COMPLEX ZDUM
  REAL CABS1
  CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
!
!***FIRST EXECUTABLE STATEMENT  CGBFA
  M = ML + MU + 1
  INFO = 0
!
!     ZERO INITIAL FILL-IN COLUMNS
!
  J0 = MU + 2
  J1 = MIN(N,M) - 1
  if (J1  <  J0) go to 30
  DO 20 JZ = J0, J1
     I0 = M + 1 - JZ
     DO 10 I = I0, ML
        ABD(I,JZ) = (0.0E0,0.0E0)
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
  JZ = J1
  JU = 0
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
  NM1 = N - 1
  if (NM1  <  1) go to 130
  DO 120 K = 1, NM1
     KP1 = K + 1
!
!        ZERO NEXT FILL-IN COLUMN
!
     JZ = JZ + 1
     if (JZ  >  N) go to 50
     if (ML  <  1) go to 50
        DO 40 I = 1, ML
           ABD(I,JZ) = (0.0E0,0.0E0)
   40       CONTINUE
   50    CONTINUE
!
!        FIND L = PIVOT INDEX
!
     LM = MIN(ML,N-K)
     L = ICAMAX(LM+1,ABD(M,K),1) + M - 1
     IPVT(K) = L + K - M
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
     if (CABS1(ABD(L,K))  ==  0.0E0) go to 100
!
!           INTERCHANGE if NECESSARY
!
        if (L  ==  M) go to 60
           T = ABD(L,K)
           ABD(L,K) = ABD(M,K)
           ABD(M,K) = T
   60       CONTINUE
!
!           COMPUTE MULTIPLIERS
!
        T = -(1.0E0,0.0E0)/ABD(M,K)
        call CSCAL(LM,T,ABD(M+1,K),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
        JU = MIN(MAX(JU,MU+IPVT(K)),N)
        MM = M
        if (JU  <  KP1) go to 90
        DO 80 J = KP1, JU
           L = L - 1
           MM = MM - 1
           T = ABD(L,J)
           if (L  ==  MM) go to 70
              ABD(L,J) = ABD(MM,J)
              ABD(MM,J) = T
   70          CONTINUE
           call CAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
     go to 110
  100    CONTINUE
        INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
  IPVT(N) = N
  if (CABS1(ABD(M,N))  ==  0.0E0) INFO = N
  return
end
