subroutine CPOFA (A, LDA, N, INFO)
!
!! CPOFA factors a complex Hermitian positive definite matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2D1B
!***TYPE      COMPLEX (SPOFA-S, DPOFA-D, CPOFA-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION,
!             POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CPOFA factors a complex Hermitian positive definite matrix.
!
!     CPOFA is usually called by CPOCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (Time for CPOCO) = (1 + 18/N)*(Time for CPOFA) .
!
!     On Entry
!
!        A       COMPLEX(LDA, N)
!                the Hermitian matrix to be factored.  Only the
!                diagonal and upper triangle are used.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       an upper triangular matrix  R  so that  A =
!                CTRANS(R)*R where  CTRANS(R)  is the conjugate
!                transpose.  The strict lower triangle is unaltered.
!                If  INFO  /=  0 , the factorization is not complete.
!
!        INFO    INTEGER
!                = 0  for normal return.
!                = K  signals an error condition.  The leading minor
!                     of order  K  is not positive definite.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CDOTC
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CPOFA
  INTEGER LDA,N,INFO
  COMPLEX A(LDA,*)
!
  COMPLEX CDOTC,T
  REAL S
  INTEGER J,JM1,K
!***FIRST EXECUTABLE STATEMENT  CPOFA
     DO 30 J = 1, N
        INFO = J
        S = 0.0E0
        JM1 = J - 1
        if (JM1  <  1) go to 20
        DO 10 K = 1, JM1
           T = A(K,J) - CDOTC(K-1,A(1,K),1,A(1,J),1)
           T = T/A(K,K)
           A(K,J) = T
           S = S + REAL(T*CONJG(T))
   10       CONTINUE
   20       CONTINUE
        S = REAL(A(J,J)) - S
        if (S  <=  0.0E0 .OR. AIMAG(A(J,J))  /=  0.0E0) go to 40
        A(J,J) = CMPLX(SQRT(S),0.0E0)
   30    CONTINUE
     INFO = 0
   40 CONTINUE
  return
end
