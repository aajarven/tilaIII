subroutine SPPFA (AP, N, INFO)
!
!! SPPFA factors a real symmetric positive definite matrix in packed form.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1B
!***TYPE      SINGLE PRECISION (SPPFA-S, DPPFA-D, CPPFA-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, PACKED,
!             POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     SPPFA factors a real symmetric positive definite matrix
!     stored in packed form.
!
!     SPPFA is usually called by SPPCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (Time for SPPCO) = (1 + 18/N)*(Time for SPPFA) .
!
!     On Entry
!
!        AP      REAL (N*(N+1)/2)
!                the packed form of a symmetric matrix  A .  The
!                columns of the upper triangle are stored sequentially
!                in a one-dimensional array of length  N*(N+1)/2 .
!                See comments below for details.
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        AP      an upper triangular matrix  R , stored in packed
!                form, so that  A = TRANS(R)*R .
!
!        INFO    INTEGER
!                = 0  for normal return.
!                = K  if the leading minor of order  K  is not
!                     positive definite.
!
!
!     Packed Storage
!
!          The following program segment will pack the upper
!          triangle of a symmetric matrix.
!
!                K = 0
!                DO 20 J = 1, N
!                   DO 10 I = 1, J
!                      K = K + 1
!                      AP(K) = A(I,J)
!             10    CONTINUE
!             20 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  SDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SPPFA
  INTEGER N,INFO
  REAL AP(*)
!
  REAL SDOT,T
  REAL S
  INTEGER J,JJ,JM1,K,KJ,KK
!***FIRST EXECUTABLE STATEMENT  SPPFA
     JJ = 0
     DO 30 J = 1, N
        INFO = J
        S = 0.0E0
        JM1 = J - 1
        KJ = JJ
        KK = 0
        if (JM1  <  1) go to 20
        DO 10 K = 1, JM1
           KJ = KJ + 1
           T = AP(KJ) - SDOT(K-1,AP(KK+1),1,AP(JJ+1),1)
           KK = KK + K
           T = T/AP(KK)
           AP(KJ) = T
           S = S + T*T
   10       CONTINUE
   20       CONTINUE
        JJ = JJ + J
        S = AP(JJ) - S
        if (S  <=  0.0E0) go to 40
        AP(JJ) = SQRT(S)
   30    CONTINUE
     INFO = 0
   40 CONTINUE
  return
end
