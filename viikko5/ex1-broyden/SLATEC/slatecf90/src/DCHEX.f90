subroutine DCHEX (R, LDR, P, K, L, Z, LDZ, NZ, C, S, JOB)
!
!! DCHEX updates the Cholesky factorization  A=TRANS(R)*R  of a ...
!            positive definite matrix A of order P under diagonal ...
!            permutations of the form  TRANS(E)*A*E, where E is a ...
!            permutation matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D7B
!***TYPE      DOUBLE PRECISION (SCHEX-S, DCHEX-D, CCHEX-C)
!***KEYWORDS  CHOLESKY DECOMPOSITION, EXCHANGE, LINEAR ALGEBRA, LINPACK,
!             MATRIX, POSITIVE DEFINITE
!***AUTHOR  Stewart, G. W., (U. of Maryland)
!***DESCRIPTION
!
!     DCHEX updates the Cholesky factorization
!
!                   A = TRANS(R)*R
!
!     of a positive definite matrix A of order P under diagonal
!     permutations of the form
!
!                   TRANS(E)*A*E
!
!     where E is a permutation matrix.  Specifically, given
!     an upper triangular matrix R and a permutation matrix
!     E (which is specified by K, L, and JOB), DCHEX determines
!     an orthogonal matrix U such that
!
!                           U*R*E = RR,
!
!     where RR is upper triangular.  At the users option, the
!     transformation U will be multiplied into the array Z.
!     If A = TRANS(X)*X, so that R is the triangular part of the
!     QR factorization of X, then RR is the triangular part of the
!     QR factorization of X*E, i.e. X with its columns permuted.
!     For a less terse description of what DCHEX does and how
!     it may be applied, see the LINPACK guide.
!
!     The matrix Q is determined as the product U(L-K)*...*U(1)
!     of plane rotations of the form
!
!                           (    C(I)       S(I) )
!                           (                    ) ,
!                           (    -S(I)      C(I) )
!
!     where C(I) is double precision.  The rows these rotations operate
!     on are described below.
!
!     There are two types of permutations, which are determined
!     by the value of JOB.
!
!     1. Right circular shift (JOB = 1).
!
!         The columns are rearranged in the following order.
!
!                1,...,K-1,L,K,K+1,...,L-1,L+1,...,P.
!
!         U is the product of L-K rotations U(I), where U(I)
!         acts in the (L-I,L-I+1)-plane.
!
!     2. Left circular shift (JOB = 2).
!         The columns are rearranged in the following order
!
!                1,...,K-1,K+1,K+2,...,L,K,L+1,...,P.
!
!         U is the product of L-K rotations U(I), where U(I)
!         acts in the (K+I-1,K+I)-plane.
!
!     On Entry
!
!         R      DOUBLE PRECISION(LDR,P), where LDR  >=  P.
!                R contains the upper triangular factor
!                that is to be updated.  Elements of R
!                below the diagonal are not referenced.
!
!         LDR    INTEGER.
!                LDR is the leading dimension of the array R.
!
!         P      INTEGER.
!                P is the order of the matrix R.
!
!         K      INTEGER.
!                K is the first column to be permuted.
!
!         L      INTEGER.
!                L is the last column to be permuted.
!                L must be strictly greater than K.
!
!         Z      DOUBLE PRECISION(LDZ,N)Z), where LDZ  >=  P.
!                Z is an array of NZ P-vectors into which the
!                transformation U is multiplied.  Z is
!                not referenced if NZ = 0.
!
!         LDZ    INTEGER.
!                LDZ is the leading dimension of the array Z.
!
!         NZ     INTEGER.
!                NZ is the number of columns of the matrix Z.
!
!         JOB    INTEGER.
!                JOB determines the type of permutation.
!                       JOB = 1  right circular shift.
!                       JOB = 2  left circular shift.
!
!     On Return
!
!         R      contains the updated factor.
!
!         Z      contains the updated matrix Z.
!
!         C      DOUBLE PRECISION(P).
!                C contains the cosines of the transforming rotations.
!
!         S      DOUBLE PRECISION(P).
!                S contains the sines of the transforming rotations.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DROTG
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DCHEX
  INTEGER LDR,P,K,L,LDZ,NZ,JOB
  DOUBLE PRECISION R(LDR,*),Z(LDZ,*),S(*)
  DOUBLE PRECISION C(*)
!
  INTEGER I,II,IL,IU,J,JJ,KM1,KP1,LMK,LM1
  DOUBLE PRECISION T
!
!     INITIALIZE
!
!***FIRST EXECUTABLE STATEMENT  DCHEX
  KM1 = K - 1
  KP1 = K + 1
  LMK = L - K
  LM1 = L - 1
!
!     PERFORM THE APPROPRIATE TASK.
!
  go to (10,130), JOB
!
!     RIGHT CIRCULAR SHIFT.
!
   10 CONTINUE
!
!        REORDER THE COLUMNS.
!
     DO 20 I = 1, L
        II = L - I + 1
        S(I) = R(II,L)
   20    CONTINUE
     DO 40 JJ = K, LM1
        J = LM1 - JJ + K
        DO 30 I = 1, J
           R(I,J+1) = R(I,J)
   30       CONTINUE
        R(J+1,J+1) = 0.0D0
   40    CONTINUE
     if (K  ==  1) go to 60
        DO 50 I = 1, KM1
           II = L - I + 1
           R(I,K) = S(II)
   50       CONTINUE
   60    CONTINUE
!
!        CALCULATE THE ROTATIONS.
!
     T = S(1)
     DO 70 I = 1, LMK
        call DROTG(S(I+1),T,C(I),S(I))
        T = S(I+1)
   70    CONTINUE
     R(K,K) = T
     DO 90 J = KP1, P
        IL = MAX(1,L-J+1)
        DO 80 II = IL, LMK
           I = L - II
           T = C(II)*R(I,J) + S(II)*R(I+1,J)
           R(I+1,J) = C(II)*R(I+1,J) - S(II)*R(I,J)
           R(I,J) = T
   80       CONTINUE
   90    CONTINUE
!
!        if REQUIRED, APPLY THE TRANSFORMATIONS TO Z.
!
     if (NZ  <  1) go to 120
     DO 110 J = 1, NZ
        DO 100 II = 1, LMK
           I = L - II
           T = C(II)*Z(I,J) + S(II)*Z(I+1,J)
           Z(I+1,J) = C(II)*Z(I+1,J) - S(II)*Z(I,J)
           Z(I,J) = T
  100       CONTINUE
  110    CONTINUE
  120    CONTINUE
  go to 260
!
!     LEFT CIRCULAR SHIFT
!
  130 CONTINUE
!
!        REORDER THE COLUMNS
!
     DO 140 I = 1, K
        II = LMK + I
        S(II) = R(I,K)
  140    CONTINUE
     DO 160 J = K, LM1
        DO 150 I = 1, J
           R(I,J) = R(I,J+1)
  150       CONTINUE
        JJ = J - KM1
        S(JJ) = R(J+1,J+1)
  160    CONTINUE
     DO 170 I = 1, K
        II = LMK + I
        R(I,L) = S(II)
  170    CONTINUE
     DO 180 I = KP1, L
        R(I,L) = 0.0D0
  180    CONTINUE
!
!        REDUCTION LOOP.
!
     DO 220 J = K, P
        if (J  ==  K) go to 200
!
!              APPLY THE ROTATIONS.
!
           IU = MIN(J-1,L-1)
           DO 190 I = K, IU
              II = I - K + 1
              T = C(II)*R(I,J) + S(II)*R(I+1,J)
              R(I+1,J) = C(II)*R(I+1,J) - S(II)*R(I,J)
              R(I,J) = T
  190          CONTINUE
  200       CONTINUE
        if (J  >=  L) go to 210
           JJ = J - K + 1
           T = S(JJ)
           call DROTG(R(J,J),T,C(JJ),S(JJ))
  210       CONTINUE
  220    CONTINUE
!
!        APPLY THE ROTATIONS TO Z.
!
     if (NZ  <  1) go to 250
     DO 240 J = 1, NZ
        DO 230 I = K, LM1
           II = I - KM1
           T = C(II)*Z(I,J) + S(II)*Z(I+1,J)
           Z(I+1,J) = C(II)*Z(I+1,J) - S(II)*Z(I,J)
           Z(I,J) = T
  230       CONTINUE
  240    CONTINUE
  250    CONTINUE
  260 CONTINUE
  return
end
