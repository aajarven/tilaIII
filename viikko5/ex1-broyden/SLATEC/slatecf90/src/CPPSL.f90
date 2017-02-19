subroutine CPPSL (AP, N, B)
!
!! CPPSL solves the complex Hermitian positive definite system using ...
!            the factors computed by CPPCO or CPPFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2D1B
!***TYPE      COMPLEX (SPPSL-S, DPPSL-D, CPPSL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, PACKED,
!             POSITIVE DEFINITE, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CPPSL solves the complex Hermitian positive definite system
!     A * X = B
!     using the factors computed by CPPCO or CPPFA.
!
!     On Entry
!
!        AP      COMPLEX (N*(N+1)/2)
!                the output from CPPCO or CPPFA.
!
!        N       INTEGER
!                the order of the matrix  A .
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
!           call CPPCO(AP,N,RCOND,Z,INFO)
!           if (RCOND is too small .OR. INFO  /=  0) go to ...
!           DO 10 J = 1, P
!              call CPPSL(AP,N,C(1,J))
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CDOTC
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CPPSL
  INTEGER N
  COMPLEX AP(*),B(*)
!
  COMPLEX CDOTC,T
  INTEGER K,KB,KK
!***FIRST EXECUTABLE STATEMENT  CPPSL
  KK = 0
  DO 10 K = 1, N
     T = CDOTC(K-1,AP(KK+1),1,B(1),1)
     KK = KK + K
     B(K) = (B(K) - T)/AP(KK)
   10 CONTINUE
  DO 20 KB = 1, N
     K = N + 1 - KB
     B(K) = B(K)/AP(KK)
     KK = KK - K
     T = -B(K)
     call CAXPY(K-1,T,AP(KK+1),1,B(1),1)
   20 CONTINUE
  return
end
