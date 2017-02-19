subroutine SPBSL (ABD, LDA, N, M, B)
!
!! SPBSL solves a real symmetric positive definite band system using the ...
!  factors computed by SPBCO or SPBFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B2
!***TYPE      SINGLE PRECISION (SPBSL-S, DPBSL-D, CPBSL-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX,
!             POSITIVE DEFINITE, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     SPBSL solves the real symmetric positive definite band
!     system  A*X = B
!     using the factors computed by SPBCO or SPBFA.
!
!     On Entry
!
!        ABD     REAL(LDA, N)
!                the output from SPBCO or SPBFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        M       INTEGER
!                the number of diagonals above the main diagonal.
!
!        B       REAL(N)
!                the right hand side vector.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains
!        a zero on the diagonal.  Technically, this indicates
!        singularity, but it is usually caused by improper subroutine
!        arguments.  It will not occur if the subroutines are called
!        correctly and  INFO  ==  0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           call SPBCO(ABD,LDA,N,RCOND,Z,INFO)
!           if (RCOND is too small .OR. INFO  /=  0) go to ...
!           DO 10 J = 1, P
!              call SPBSL(ABD,LDA,N,C(1,J))
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  SAXPY, SDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SPBSL
  INTEGER LDA,N,M
  REAL ABD(LDA,*),B(*)
!
  REAL SDOT,T
  INTEGER K,KB,LA,LB,LM
!
!     SOLVE TRANS(R)*Y = B
!
!***FIRST EXECUTABLE STATEMENT  SPBSL
  DO 10 K = 1, N
     LM = MIN(K-1,M)
     LA = M + 1 - LM
     LB = K - LM
     T = SDOT(LM,ABD(LA,K),1,B(LB),1)
     B(K) = (B(K) - T)/ABD(M+1,K)
   10 CONTINUE
!
!     SOLVE R*X = Y
!
  DO 20 KB = 1, N
     K = N + 1 - KB
     LM = MIN(K-1,M)
     LA = M + 1 - LM
     LB = K - LM
     B(K) = B(K)/ABD(M+1,K)
     T = -B(K)
     call SAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   20 CONTINUE
  return
end
