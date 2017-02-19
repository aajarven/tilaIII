subroutine CGBSL (ABD, LDA, N, ML, MU, IPVT, B, JOB)
!
!! CGBSL solves the complex band system A*X=B or CTRANS(A)*X=B using ...
!            the factors computed by CGBCO or CGBFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2C2
!***TYPE      COMPLEX (SGBSL-S, DGBSL-D, CGBSL-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CGBSL solves the complex band system
!     A * X = B  or  CTRANS(A) * X = B
!     using the factors computed by CGBCO or CGBFA.
!
!     On Entry
!
!        ABD     COMPLEX(LDA, N)
!                the output from CGBCO or CGBFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
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
!                the pivot vector from CGBCO or CGBFA.
!
!        B       COMPLEX(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
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
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if CGBCO has set RCOND  >  0.0
!        or CGBFA has set INFO  ==  0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           call CGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!           if (RCOND is too small) go to ...
!           DO 10 J = 1, P
!              call CGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CDOTC
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CGBSL
  INTEGER LDA,N,ML,MU,IPVT(*),JOB
  COMPLEX ABD(LDA,*),B(*)
!
  COMPLEX CDOTC,T
  INTEGER K,KB,L,LA,LB,LM,M,NM1
!***FIRST EXECUTABLE STATEMENT  CGBSL
  M = MU + ML + 1
  NM1 = N - 1
  if (JOB  /=  0) go to 50
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE L*Y = B
!
     if (ML  ==  0) go to 30
     if (NM1  <  1) go to 30
        DO 20 K = 1, NM1
           LM = MIN(ML,N-K)
           L = IPVT(K)
           T = B(L)
           if (L  ==  K) go to 10
              B(L) = B(K)
              B(K) = T
   10          CONTINUE
           call CAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
!
!        NOW SOLVE  U*X = Y
!
     DO 40 KB = 1, N
        K = N + 1 - KB
        B(K) = B(K)/ABD(M,K)
        LM = MIN(K,M) - 1
        LA = M - LM
        LB = K - LM
        T = -B(K)
        call CAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
  go to 100
   50 CONTINUE
!
!        JOB = NONZERO, SOLVE  CTRANS(A) * X = B
!        FIRST SOLVE  CTRANS(U)*Y = B
!
     DO 60 K = 1, N
        LM = MIN(K,M) - 1
        LA = M - LM
        LB = K - LM
        T = CDOTC(LM,ABD(LA,K),1,B(LB),1)
        B(K) = (B(K) - T)/CONJG(ABD(M,K))
   60    CONTINUE
!
!        NOW SOLVE CTRANS(L)*X = Y
!
     if (ML  ==  0) go to 90
     if (NM1  <  1) go to 90
        DO 80 KB = 1, NM1
           K = N - KB
           LM = MIN(ML,N-K)
           B(K) = B(K) + CDOTC(LM,ABD(M+1,K),1,B(K+1),1)
           L = IPVT(K)
           if (L  ==  K) go to 70
              T = B(L)
              B(L) = B(K)
              B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
  return
end
