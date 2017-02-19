subroutine DGESL (A, LDA, N, IPVT, B, JOB)
!
!! DGESL solves the real system A*X=B or TRANS(A)*X=B using the ...
!            factors computed by DGECO or DGEFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2A1
!***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGESL solves the double precision system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DGECO or DGEFA.
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the output from DGECO or DGEFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        IPVT    INTEGER(N)
!                the pivot vector from DGECO or DGEFA.
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B  where
!                            TRANS(A)  is the transpose.
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
!        called correctly and if DGECO has set RCOND  >  0.0
!        or DGEFA has set INFO  ==  0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           call DGECO(A,LDA,N,IPVT,RCOND,Z)
!           if (RCOND is too small) go to ...
!           DO 10 J = 1, P
!              call DGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGESL
  INTEGER LDA,N,IPVT(*),JOB
  DOUBLE PRECISION A(LDA,*),B(*)
!
  DOUBLE PRECISION DDOT,T
  INTEGER K,KB,L,NM1
!***FIRST EXECUTABLE STATEMENT  DGESL
  NM1 = N - 1
  if (JOB  /=  0) go to 50
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE  L*Y = B
!
     if (NM1  <  1) go to 30
     DO 20 K = 1, NM1
        L = IPVT(K)
        T = B(L)
        if (L  ==  K) go to 10
           B(L) = B(K)
           B(K) = T
   10       CONTINUE
        call DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
!
!        NOW SOLVE  U*X = Y
!
     DO 40 KB = 1, N
        K = N + 1 - KB
        B(K) = B(K)/A(K,K)
        T = -B(K)
        call DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
  go to 100
   50 CONTINUE
!
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B
!
     DO 60 K = 1, N
        T = DDOT(K-1,A(1,K),1,B(1),1)
        B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
!
!        NOW SOLVE TRANS(L)*X = Y
!
     if (NM1  <  1) go to 90
     DO 80 KB = 1, NM1
        K = N - KB
        B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
        L = IPVT(K)
        if (L  ==  K) go to 70
           T = B(L)
           B(L) = B(K)
           B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
  return
end
