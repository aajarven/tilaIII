subroutine CNBSL (ABE, LDA, N, ML, MU, IPVT, B, JOB)
!
!! CNBSL solves a complex band system with factors computed by CNBCO or CNBFA.
!
!***LIBRARY   SLATEC
!***CATEGORY  D2C2
!***TYPE      COMPLEX (SNBSL-S, DNBSL-D, CNBSL-C)
!***KEYWORDS  BANDED, LINEAR EQUATIONS, NONSYMMETRIC, SOLVE
!***AUTHOR  Voorhees, E. A., (LANL)
!***DESCRIPTION
!
!     CNBSL solves the complex band system
!     A * X = B  or  CTRANS(A) * X = B
!     using the factors computed by CNBCO or CNBFA.
!
!     On Entry
!
!        ABE     COMPLEX(LDA, NC)
!                the output from CNBCO or CNBFA.
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
!                the pivot vector from CNBCO or CNBFA.
!
!        B       COMPLEX(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B .
!                = nonzero   to solve  CTRANS(A)*X = B , where
!                            CTRANS(A)  is the conjugate transpose.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA.  It will not occur if the subroutines are
!        called correctly and if CNBCO has set RCOND  >  0.0
!        or CNBFA has set INFO  ==  0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           call CNBCO(ABE,LDA,N,ML,MU,IPVT,RCOND,Z)
!           if (RCOND is too small) go to ...
!           DO 10 J = 1, P
!             call CNBSL(ABE,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CDOTC
!***REVISION HISTORY  (YYMMDD)
!   800730  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CNBSL
  INTEGER LDA,N,ML,MU,IPVT(*),JOB
  COMPLEX ABE(LDA,*),B(*)
!
  COMPLEX CDOTC,T
  INTEGER K,KB,L,LB,LDB,LM,M,MLM,NM1
!***FIRST EXECUTABLE STATEMENT  CNBSL
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
        call CAXPY(LM,T,ABE(K+LM,MLM),LDB,B(K+1),1)
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
      call CAXPY(LM,T,ABE(K-1,ML+2),LDB,B(LB),1)
   40   CONTINUE
  go to 100
   50 CONTINUE
!
!       JOB = NONZERO, SOLVE CTRANS(A) * X = B
!       FIRST SOLVE  CTRANS(U)*Y = B
!
    DO 60 K = 1, N
      LM = MIN(K,M) - 1
      LB = K - LM
      T = CDOTC(LM,ABE(K-1,ML+2),LDB,B(LB),1)
      B(K) = (B(K) - T)/CONJG(ABE(K,ML+1))
   60   CONTINUE
!
!       NOW SOLVE CTRANS(L)*X = Y
!
    if (ML  ==  0) go to 90
    if (NM1  <  1) go to 90
      DO 80 KB = 1, NM1
        K = N - KB
        LM = MIN(ML,N-K)
        MLM = ML - (LM - 1)
        B(K) = B(K) + CDOTC(LM,ABE(K+LM,MLM),LDB,B(K+1),1)
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
