subroutine CPBFA (ABD, LDA, N, M, INFO)
!
!! CPBFA factors a complex Hermitian positive definite band matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2D2
!***TYPE      COMPLEX (SPBFA-S, DPBFA-D, CPBFA-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION,
!             POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CPBFA factors a complex Hermitian positive definite matrix
!     stored in band form.
!
!     CPBFA is usually called by CPBCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!
!     On Entry
!
!        ABD     COMPLEX(LDA, N)
!                the matrix to be factored.  The columns of the upper
!                triangle are stored in the columns of ABD and the
!                diagonals of the upper triangle are stored in the
!                rows of ABD .  See the comments below for details.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!                LDA must be  >=  M + 1 .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        M       INTEGER
!                the number of diagonals above the main diagonal.
!                0  <=  M  <  N .
!
!     On Return
!
!        ABD     an upper triangular matrix  R , stored in band
!                form, so that  A = CTRANS(R)*R .
!
!        INFO    INTEGER
!                = 0  for normal return.
!                = K  if the leading minor of order  K  is not
!                     positive definite.
!
!     Band Storage
!
!           If  A  is a Hermitian positive definite band matrix,
!           the following program segment will set up the input.
!
!                   M = (band width above diagonal)
!                   DO 20 J = 1, N
!                      I1 = MAX(1, J-M)
!                      DO 10 I = I1, J
!                         K = I-J+M+1
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CDOTC
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CPBFA
  INTEGER LDA,N,M,INFO
  COMPLEX ABD(LDA,*)
!
  COMPLEX CDOTC,T
  REAL S
  INTEGER IK,J,JK,K,MU
!***FIRST EXECUTABLE STATEMENT  CPBFA
     DO 30 J = 1, N
        INFO = J
        S = 0.0E0
        IK = M + 1
        JK = MAX(J-M,1)
        MU = MAX(M+2-J,1)
        if (M  <  MU) go to 20
        DO 10 K = MU, M
           T = ABD(K,J) - CDOTC(K-MU,ABD(IK,JK),1,ABD(MU,J),1)
           T = T/ABD(M+1,JK)
           ABD(K,J) = T
           S = S + REAL(T*CONJG(T))
           IK = IK - 1
           JK = JK + 1
   10       CONTINUE
   20       CONTINUE
        S = REAL(ABD(M+1,J)) - S
        if (S  <=  0.0E0 .OR. AIMAG(ABD(M+1,J))  /=  0.0E0) &
           go to 40
        ABD(M+1,J) = CMPLX(SQRT(S),0.0E0)
   30    CONTINUE
     INFO = 0
   40 CONTINUE
  return
end
