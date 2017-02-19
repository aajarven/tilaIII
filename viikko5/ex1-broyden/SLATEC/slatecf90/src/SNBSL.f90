subroutine SNBSL (ABE, LDA, N, ML, MU, IPVT, B, JOB)
!
!! SNBSL solves a real band system using the factors computed by SNBCO or SNBFA.
!
!***LIBRARY   SLATEC
!***CATEGORY  D2A2
!***TYPE      SINGLE PRECISION (SNBSL-S, DNBSL-D, CNBSL-C)
!***KEYWORDS  BANDED, LINEAR EQUATIONS, NONSYMMETRIC, SOLVE
!***AUTHOR  Voorhees, E. A., (LANL)
!***DESCRIPTION
!
!     SNBSL solves the real band system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by SNBCO or SNBFA.
!
!     On Entry
!
!        ABE     REAL(LDA, NC)
!                the output from SNBCO or SNBFA.
!                NC must be  >=  2*ML+MU+1 .
!
!        LDA     INTEGER
!                the leading dimension of the array  ABE .
!
!        N       INTEGER
!                the order of the original matrix.
!
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!
!        IPVT    INTEGER(N)
!                the pivot vector from SNBCO or SNBFA.
!
!        B       REAL(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B .
!                = nonzero   to solve  TRANS(A)*X = B , where
!                            TRANS(A)  is the transpose.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically, this indicates singularity,
!        but it is often caused by improper arguments or improper
!        setting of LDA.  It will not occur if the subroutines are
!        called correctly and if SNBCO has set RCOND  >  0.0
!        or SNBFA has set INFO  ==  0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           call SNBCO(ABE,LDA,N,ML,MU,IPVT,RCOND,Z)
!           if (RCOND is too small) go to ...
!           DO 10 J = 1, P
!             call SNBSL(ABE,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  SAXPY, SDOT
!***REVISION HISTORY  (YYMMDD)
!   800717  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SNBSL
  INTEGER LDA,N,ML,MU,IPVT(*),JOB
  REAL ABE(LDA,*),B(*)
!
  REAL SDOT,T
  INTEGER K,KB,L,LB,LDB,LM,M,MLM,NM1
!***FIRST EXECUTABLE STATEMENT  SNBSL
  M=MU+ML+1
  NM1=N-1
  LDB=1-LDA
  if ( JOB /= 0)go to 50
!
!       JOB = 0 , SOLVE  A * X = B
!       FIRST SOLVE L*Y = B
!
    if ( ML == 0)go to 30
    if ( NM1 < 1)go to 30
      DO 20 K=1,NM1
        LM=MIN(ML,N-K)
        L=IPVT(K)
        T=B(L)
        if ( L == K)go to 10
          B(L)=B(K)
          B(K)=T
   10       CONTINUE
        MLM=ML-(LM-1)
        call SAXPY(LM,T,ABE(K+LM,MLM),LDB,B(K+1),1)
   20     CONTINUE
   30   CONTINUE
!
!       NOW SOLVE  U*X = Y
!
    DO 40 KB=1,N
      K=N+1-KB
      B(K)=B(K)/ABE(K,ML+1)
      LM=MIN(K,M)-1
      LB=K-LM
      T=-B(K)
      call SAXPY(LM,T,ABE(K-1,ML+2),LDB,B(LB),1)
   40   CONTINUE
  go to 100
   50 CONTINUE
!
!       JOB = NONZERO, SOLVE TRANS(A) * X = B
!       FIRST SOLVE  TRANS(U)*Y = B
!
    DO 60 K = 1, N
      LM = MIN(K,M) - 1
      LB = K - LM
      T = SDOT(LM,ABE(K-1,ML+2),LDB,B(LB),1)
      B(K) = (B(K) - T)/ABE(K,ML+1)
   60   CONTINUE
!
!       NOW SOLVE TRANS(L)*X = Y
!
    if (ML  ==  0) go to 90
    if (NM1  <  1) go to 90
      DO 80 KB = 1, NM1
        K = N - KB
        LM = MIN(ML,N-K)
        MLM = ML - (LM - 1)
        B(K) = B(K) + SDOT(LM,ABE(K+LM,MLM),LDB,B(K+1),1)
        L = IPVT(K)
        if (L  ==  K) go to 70
          T = B(L)
          B(L) = B(K)
          B(K) = T
   70       CONTINUE
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
  return
end
