subroutine BQR (NM, N, MB, A, T, R, IERR, NV, RV)
!
!! BQR computes some of the eigenvalues of a real symmetric ...
!  matrix using the QR method with shifts of origin.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A6
!***TYPE      SINGLE PRECISION (BQR-S)
!***KEYWORDS  EIGENVALUES, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure BQR,
!     NUM. MATH. 16, 85-92(1970) by Martin, Reinsch, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 266-272(1971).
!
!     This subroutine finds the eigenvalue of smallest (usually)
!     magnitude of a REAL SYMMETRIC BAND matrix using the
!     QR algorithm with shifts of origin.  Consecutive calls
!     can be made to find further eigenvalues.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, A, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        MB is the (half) band width of the matrix, defined as the
!          number of adjacent diagonals, including the principal
!          diagonal, required to specify the non-zero portion of the
!          lower triangle of the matrix.  MB is an INTEGER variable.
!          MB must be less than or equal to N on first call.
!
!        A contains the lower triangle of the symmetric band input
!          matrix stored as an N by MB array.  Its lowest subdiagonal
!          is stored in the last N+1-MB positions of the first column,
!          its next subdiagonal in the last N+2-MB positions of the
!          second column, further subdiagonals similarly, and finally
!          its principal diagonal in the N positions of the last column.
!          Contents of storages not part of the matrix are arbitrary.
!          On a subsequent call, its output contents from the previous
!          call should be passed.  A is a two-dimensional REAL array,
!          dimensioned A(NM,MB).
!
!        T specifies the shift (of eigenvalues) applied to the diagonal
!          of A in forming the input matrix. What is actually determined
!          is the eigenvalue of A+TI (I is the identity matrix) nearest
!          to T.  On a subsequent call, the output value of T from the
!          previous call should be passed if the next nearest eigenvalue
!          is sought.  T is a REAL variable.
!
!        R should be specified as zero on the first call, and as its
!          output value from the previous call on a subsequent call.
!          It is used to determine when the last row and column of
!          the transformed band matrix can be regarded as negligible.
!          R is a REAL variable.
!
!        NV must be set to the dimension of the array parameter RV
!          as declared in the calling program dimension statement.
!          NV is an INTEGER variable.
!
!     On OUTPUT
!
!        A contains the transformed band matrix.  The matrix A+TI
!          derived from the output parameters is similar to the
!          input A+TI to within rounding errors.  Its last row and
!          column are null (if IERR is zero).
!
!        T contains the computed eigenvalue of A+TI (if IERR is zero),
!          where I is the identity matrix.
!
!        R contains the maximum of its input value and the norm of the
!          last column of the input matrix A.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after a total of 30 iterations.
!
!        RV is a one-dimensional REAL array of dimension NV which is
!          at least (2*MB**2+4*MB-3), used for temporary storage.  The
!          first (3*MB-2) locations correspond to the ALGOL array B,
!          the next (2*MB-1) locations correspond to the ALGOL array H,
!          and the final (2*MB**2-MB) locations correspond to the MB
!          by (2*MB-1) ALGOL array U.
!
!     NOTE. For a subsequent call, N should be replaced by N-1, but
!     MB should not be altered even when it exceeds the current N.
!
!     Calls PYTHAG(A,B) for SQRT(A**2 + B**2).
!
!     Questions and comments should be directed to B. S. Garbow,
!     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  PYTHAG
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BQR
!
  INTEGER I,J,K,L,M,N,II,IK,JK,JM,KJ,KK,KM,LL,MB,MK,MN,MZ
  INTEGER M1,M2,M3,M4,NI,NM,NV,ITS,KJ1,M21,M31,IERR,IMULT
  REAL A(NM,*),RV(*)
  REAL F,G,Q,R,S,T,SCALE
  REAL PYTHAG
!
!***FIRST EXECUTABLE STATEMENT  BQR
  IERR = 0
  M1 = MIN(MB,N)
  M = M1 - 1
  M2 = M + M
  M21 = M2 + 1
  M3 = M21 + M
  M31 = M3 + 1
  M4 = M31 + M2
  MN = M + N
  MZ = MB - M1
  ITS = 0
!     .......... TEST FOR CONVERGENCE ..........
   40 G = A(N,MB)
  if (M  ==  0) go to 360
  F = 0.0E0
!
  DO 50 K = 1, M
     MK = K + MZ
     F = F + ABS(A(N,MK))
   50 CONTINUE
!
  if (ITS  ==  0 .AND. F  >  R) R = F
  if (R + F  <=  R) go to 360
  if (ITS  ==  30) go to 1000
  ITS = ITS + 1
!     .......... FORM SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
  if (F  >  0.25E0 * R .AND. ITS  <  5) go to 90
  F = A(N,MB-1)
  if (F  ==  0.0E0) go to 70
  Q = (A(N-1,MB) - G) / (2.0E0 * F)
  S = PYTHAG(Q,1.0E0)
  G = G - F / (Q + SIGN(S,Q))
   70 T = T + G
!
  DO 80 I = 1, N
   80 A(I,MB) = A(I,MB) - G
!
   90 DO 100 K = M31, M4
  100 RV(K) = 0.0E0
!
  DO 350 II = 1, MN
     I = II - M
     NI = N - II
     if (NI  <  0) go to 230
!     .......... FORM COLUMN OF SHIFTED MATRIX A-G*I ..........
     L = MAX(1,2-I)
!
     DO 110 K = 1, M3
  110    RV(K) = 0.0E0
!
     DO 120 K = L, M1
        KM = K + M
        MK = K + MZ
        RV(KM) = A(II,MK)
  120    CONTINUE
!
     LL = MIN(M,NI)
     if (LL  ==  0) go to 135
!
     DO 130 K = 1, LL
        KM = K + M21
        IK = II + K
        MK = MB - K
        RV(KM) = A(IK,MK)
  130    CONTINUE
!     .......... PRE-MULTIPLY WITH HOUSEHOLDER REFLECTIONS ..........
  135    LL = M2
     IMULT = 0
!     .......... MULTIPLICATION PROCEDURE ..........
  140    KJ = M4 - M1
!
     DO 170 J = 1, LL
        KJ = KJ + M1
        JM = J + M3
        if (RV(JM)  ==  0.0E0) go to 170
        F = 0.0E0
!
        DO 150 K = 1, M1
           KJ = KJ + 1
           JK = J + K - 1
           F = F + RV(KJ) * RV(JK)
  150       CONTINUE
!
        F = F / RV(JM)
        KJ = KJ - M1
!
        DO 160 K = 1, M1
           KJ = KJ + 1
           JK = J + K - 1
           RV(JK) = RV(JK) - RV(KJ) * F
  160       CONTINUE
!
        KJ = KJ - M1
  170    CONTINUE
!
     if (IMULT  /=  0) go to 280
!     .......... HOUSEHOLDER REFLECTION ..........
     F = RV(M21)
     S = 0.0E0
     RV(M4) = 0.0E0
     SCALE = 0.0E0
!
     DO 180 K = M21, M3
  180    SCALE = SCALE + ABS(RV(K))
!
     if (SCALE  ==  0.0E0) go to 210
!
     DO 190 K = M21, M3
  190    S = S + (RV(K)/SCALE)**2
!
     S = SCALE * SCALE * S
     G = -SIGN(SQRT(S),F)
     RV(M21) = G
     RV(M4) = S - F * G
     KJ = M4 + M2 * M1 + 1
     RV(KJ) = F - G
!
     DO 200 K = 2, M1
        KJ = KJ + 1
        KM = K + M2
        RV(KJ) = RV(KM)
  200    CONTINUE
!     .......... SAVE COLUMN OF TRIANGULAR FACTOR R ..........
  210    DO 220 K = L, M1
        KM = K + M
        MK = K + MZ
        A(II,MK) = RV(KM)
  220    CONTINUE
!
  230    L = MAX(1,M1+1-I)
     if (I  <=  0) go to 300
!     .......... PERFORM ADDITIONAL STEPS ..........
     DO 240 K = 1, M21
  240    RV(K) = 0.0E0
!
     LL = MIN(M1,NI+M1)
!     .......... GET ROW OF TRIANGULAR FACTOR R ..........
     DO 250 KK = 1, LL
        K = KK - 1
        KM = K + M1
        IK = I + K
        MK = MB - K
        RV(KM) = A(IK,MK)
  250    CONTINUE
!     .......... POST-MULTIPLY WITH HOUSEHOLDER REFLECTIONS ..........
     LL = M1
     IMULT = 1
     go to 140
!     .......... STORE COLUMN OF NEW A MATRIX ..........
  280    DO 290 K = L, M1
        MK = K + MZ
        A(I,MK) = RV(K)
  290    CONTINUE
!     .......... UPDATE HOUSEHOLDER REFLECTIONS ..........
  300    if (L  >  1) L = L - 1
     KJ1 = M4 + L * M1
!
     DO 320 J = L, M2
        JM = J + M3
        RV(JM) = RV(JM+1)
!
        DO 320 K = 1, M1
           KJ1 = KJ1 + 1
           KJ = KJ1 - M1
           RV(KJ) = RV(KJ1)
  320    CONTINUE
!
  350 CONTINUE
!
  go to 40
!     .......... CONVERGENCE ..........
  360 T = T + G
!
  DO 380 I = 1, N
  380 A(I,MB) = A(I,MB) - G
!
  DO 400 K = 1, M1
     MK = K + MZ
     A(N,MK) = 0.0E0
  400 CONTINUE
!
  go to 1001
!     .......... SET ERROR -- NO CONVERGENCE TO
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = N
 1001 RETURN
end
