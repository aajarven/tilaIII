subroutine DPOFA (A, LDA, N, INFO)
!
!! DPOFA factors a real symmetric positive definite matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1B
!***TYPE      DOUBLE PRECISION (SPOFA-S, DPOFA-D, CPOFA-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION,
!             POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DPOFA factors a double precision symmetric positive definite
!     matrix.
!
!     DPOFA is usually called by DPOCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (time for DPOCO) = (1 + 18/N)*(time for DPOFA) .
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the symmetric matrix to be factored.  Only the
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
!        A       an upper triangular matrix  R  so that  A = TRANS(R)*R
!                where  TRANS(R)  is the transpose.
!                The strict lower triangle is unaltered.
!                If  INFO  /=  0 , the factorization is not complete.
!
!        INFO    INTEGER
!                = 0  for normal return.
!                = K  signals an error condition.  The leading minor
!                     of order  K  is not positive definite.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DPOFA
  INTEGER LDA,N,INFO
  DOUBLE PRECISION A(LDA,*)
!
  DOUBLE PRECISION DDOT,T
  DOUBLE PRECISION S
  INTEGER J,JM1,K
!***FIRST EXECUTABLE STATEMENT  DPOFA
     DO 30 J = 1, N
        INFO = J
        S = 0.0D0
        JM1 = J - 1
        if (JM1  <  1) go to 20
        DO 10 K = 1, JM1
           T = A(K,J) - DDOT(K-1,A(1,K),1,A(1,J),1)
           T = T/A(K,K)
           A(K,J) = T
           S = S + T*T
   10       CONTINUE
   20       CONTINUE
        S = A(J,J) - S
        if (S  <=  0.0D0) go to 40
        A(J,J) = SQRT(S)
   30    CONTINUE
     INFO = 0
   40 CONTINUE
  return
end
