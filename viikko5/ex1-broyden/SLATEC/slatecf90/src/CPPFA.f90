subroutine CPPFA (AP, N, INFO)
!
!! CPPFA factors a complex Hermitian positive definite matrix in packed form.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2D1B
!***TYPE      COMPLEX (SPPFA-S, DPPFA-D, CPPFA-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, PACKED,
!             POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CPPFA factors a complex Hermitian positive definite matrix
!     stored in packed form.
!
!     CPPFA is usually called by CPPCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (Time for CPPCO) = (1 + 18/N)*(Time for CPPFA) .
!
!     On Entry
!
!        AP      COMPLEX (N*(N+1)/2)
!                the packed form of a Hermitian matrix  A .  The
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
!                form, so that  A = CTRANS(R)*R .
!
!        INFO    INTEGER
!                = 0  for normal return.
!                = K  If the leading minor of order  K  is not
!                     positive definite.
!
!
!     Packed Storage
!
!          The following program segment will pack the upper
!          triangle of a Hermitian matrix.
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
!***ROUTINES CALLED  CDOTC
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CPPFA
  INTEGER N,INFO
  COMPLEX AP(*)
!
  COMPLEX CDOTC,T
  REAL S
  INTEGER J,JJ,JM1,K,KJ,KK
!***FIRST EXECUTABLE STATEMENT  CPPFA
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
           T = AP(KJ) - CDOTC(K-1,AP(KK+1),1,AP(JJ+1),1)
           KK = KK + K
           T = T/AP(KK)
           AP(KJ) = T
           S = S + REAL(T*CONJG(T))
   10       CONTINUE
   20       CONTINUE
        JJ = JJ + J
        S = REAL(AP(JJ)) - S
        if (S  <=  0.0E0 .OR. AIMAG(AP(JJ))  /=  0.0E0) go to 40
        AP(JJ) = CMPLX(SQRT(S),0.0E0)
   30    CONTINUE
     INFO = 0
   40 CONTINUE
  return
end
