subroutine CNBFA (ABE, LDA, N, ML, MU, IPVT, INFO)
!
!! CNBFA factors a band matrix by elimination.
!
!***LIBRARY   SLATEC
!***CATEGORY  D2C2
!***TYPE      COMPLEX (SNBFA-S, DNBFA-D, CNBFA-C)
!***KEYWORDS  BANDED, LINEAR EQUATIONS, MATRIX FACTORIZATION,
!             NONSYMMETRIC
!***AUTHOR  Voorhees, E. A., (LANL)
!***DESCRIPTION
!
!    CNBFA factors a complex band matrix by elimination.
!
!     CNBFA is usually called by CNBCO, but it can be called
!     directly with a saving in time if RCOND is not needed.
!
!     On Entry
!
!        ABE     COMPLEX(LDA, NC)
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
!                the factorization can be written  A = L*U  where
!                L is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                =0  normal value
!                =K  if  U(K,K)  ==  0.0 .  This is not an error
!                condition for this subroutine, but it does
!                indicate that CNBSL will divide by zero if
!                called.  Use RCOND in CNBCO for a reliable
!                indication of singularity.
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
!***ROUTINES CALLED  CAXPY, CSCAL, CSWAP, ICAMAX
!***REVISION HISTORY  (YYMMDD)
!   800730  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CNBFA
  INTEGER LDA,N,ML,MU,IPVT(*),INFO
  COMPLEX ABE(LDA,*)
!
  INTEGER ML1,MB,M,N1,LDB,I,J,K,L,LM,LM1,LM2,MP,ICAMAX
  COMPLEX T
  COMPLEX ZDUM
  REAL CABS1
  CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
!
!***FIRST EXECUTABLE STATEMENT  CNBFA
  ML1=ML+1
  MB=ML+MU
  M=ML+MU+1
  N1=N-1
  LDB=LDA-1
  INFO=0
!
!     SET FILL-IN COLUMNS TO ZERO
!
  if ( N <= 1)go to 50
  if ( ML <= 0)go to 7
  DO 6 J=1,ML
    DO 5 I=1,N
      ABE(I,M+J)=(0.0E0,0.0E0)
    5   CONTINUE
    6 CONTINUE
    7 CONTINUE
!
!     GAUSSIAN ELIMINATION WITH PARTIAL ELIMINATION
!
  DO 40 K=1,N1
    LM=MIN(N-K,ML)
    LM1=LM+1
    LM2=ML1-LM
!
!     SEARCH FOR PIVOT INDEX
!
    L=-ICAMAX(LM1,ABE(LM+K,LM2),LDB)+LM1+K
    IPVT(K)=L
    MP=MIN(MB,N-K)
!
!     SWAP ROWS if NECESSARY
!
    if ( L /= K)CALL CSWAP(MP+1,ABE(K,ML1),LDA,ABE(L,ML1+K-L),LDA)
!
!     SKIP COLUMN REDUCTION if PIVOT IS ZERO
!
    if ( CABS1(ABE(K,ML1)) == 0.0E0) go to 20
!
!     COMPUTE MULTIPLIERS
!
    T=-(1.0E0,0.0E0)/ABE(K,ML1)
    call CSCAL(LM,T,ABE(LM+K,LM2),LDB)
!
!     ROW ELIMINATION WITH COLUMN INDEXING
!
    DO 10 J=1,MP
      call CAXPY (LM,ABE(K,ML1+J),ABE(LM+K,LM2),LDB,ABE(LM+K,LM2+J), &
                  LDB)
   10   CONTINUE
    go to 30
   20   CONTINUE
    INFO=K
   30   CONTINUE
   40 CONTINUE
   50 CONTINUE
  IPVT(N)=N
  if ( CABS1(ABE(N,ML1)) == 0.0E0) INFO=N
  return
end
