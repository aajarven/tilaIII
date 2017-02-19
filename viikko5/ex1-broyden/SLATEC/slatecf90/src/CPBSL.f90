subroutine CPBSL (ABD, LDA, N, M, B)
!
!! CPBSL solves the complex Hermitian positive definite band system ...
!  using the factors computed by CPBCO or CPBFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2D2
!***TYPE      COMPLEX (SPBSL-S, DPBSL-D, CPBSL-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX,
!             POSITIVE DEFINITE, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CPBSL solves the complex Hermitian positive definite band
!     system  A*X = B
!     using the factors computed by CPBCO or CPBFA.
!
!     On Entry
!
!        ABD     COMPLEX(LDA, N)
!                the output from CPBCO or CPBFA.
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
!        B       COMPLEX(N)
!                the right hand side vector.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains
!        a zero on the diagonal.  Technically this indicates
!        singularity but it is usually caused by improper subroutine
!        arguments.  It will not occur if the subroutines are called
!        correctly and  INFO  ==  0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           call CPBCO(ABD,LDA,N,RCOND,Z,INFO)
!           if (RCOND is too small .OR. INFO  /=  0) go to ...
!           DO 10 J = 1, P
!              call CPBSL(ABD,LDA,N,C(1,J))
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
!***END PROLOGUE  CPBSL
  INTEGER LDA,N,M
  COMPLEX ABD(LDA,*),B(*)
!
  COMPLEX CDOTC,T
  INTEGER K,KB,LA,LB,LM
!
!     SOLVE CTRANS(R)*Y = B
!
!***FIRST EXECUTABLE STATEMENT  CPBSL
  DO 10 K = 1, N
     LM = MIN(K-1,M)
     LA = M + 1 - LM
     LB = K - LM
     T = CDOTC(LM,ABD(LA,K),1,B(LB),1)
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
     call CAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   20 CONTINUE
  return
end
