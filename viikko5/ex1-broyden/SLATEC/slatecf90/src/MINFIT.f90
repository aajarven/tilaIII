subroutine MINFIT (NM, M, N, A, W, IP, B, IERR, RV1)
!
!! MINFIT computes the singular value decomposition of a rectangular ...
!            matrix and solve the related linear least squares problem.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D9
!***TYPE      SINGLE PRECISION (MINFIT-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure MINFIT,
!     NUM. MATH. 14, 403-420(1970) by Golub and Reinsch.
!     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
!
!     This subroutine determines, towards the solution of the linear
!                                                        T
!     system AX=B, the singular value decomposition A=USV  of a real
!                                         T
!     M by N rectangular matrix, forming U B rather than U.  Householder
!     bidiagonalization and a variant of the QR algorithm are used.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A and B, as declared in the calling
!          program dimension statement.  Note that NM must be at least
!          as large as the maximum of M and N.  NM is an INTEGER
!          variable.
!
!        M is the number of rows of A and B.  M is an INTEGER variable.
!
!        N is the number of columns of A and the order of V.  N is an
!          INTEGER variable.
!
!        A contains the rectangular coefficient matrix of the system.
!          A is a two-dimensional REAL array, dimensioned A(NM,N).
!
!        IP is the number of columns of B.  IP can be zero.
!
!        B contains the constant column matrix of the system if IP is
!          not zero.  Otherwise, B is not referenced.  B is a two-
!          dimensional REAL array, dimensioned B(NM,IP).
!
!     On OUTPUT
!
!        A has been overwritten by the matrix V (orthogonal) of the
!          decomposition in its first N rows and columns.  If an
!          error exit is made, the columns of V corresponding to
!          indices of correct singular values should be correct.
!
!        W contains the N (non-negative) singular values of A (the
!          diagonal elements of S).  They are unordered.  If an
!          error exit is made, the singular values should be correct
!          for indices IERR+1, IERR+2, ..., N.  W is a one-dimensional
!          REAL array, dimensioned W(N).
!
!                                   T
!        B has been overwritten by U B.  If an error exit is made,
!                       T
!          the rows of U B corresponding to indices of correct singular
!          values should be correct.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          K          if the K-th singular value has not been
!                     determined after 30 iterations.
!                     The singular values should be correct for
!                     indices IERR+1, IERR+2, ..., N.
!
!        RV1 is a one-dimensional REAL array used for temporary storage,
!          dimensioned RV1(N).
!
!     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
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
!***END PROLOGUE  MINFIT
!
  INTEGER I,J,K,L,M,N,II,IP,I1,KK,K1,LL,L1,M1,NM,ITS,IERR
  REAL A(NM,*),W(*),B(NM,IP),RV1(*)
  REAL C,F,G,H,S,X,Y,Z,SCALE,S1
  REAL PYTHAG
!
!***FIRST EXECUTABLE STATEMENT  MINFIT
  IERR = 0
!     .......... HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM ..........
  G = 0.0E0
  SCALE = 0.0E0
  S1 = 0.0E0
!
  DO 300 I = 1, N
     L = I + 1
     RV1(I) = SCALE * G
     G = 0.0E0
     S = 0.0E0
     SCALE = 0.0E0
     if (I  >  M) go to 210
!
     DO 120 K = I, M
  120    SCALE = SCALE + ABS(A(K,I))
!
     if (SCALE  ==  0.0E0) go to 210
!
     DO 130 K = I, M
        A(K,I) = A(K,I) / SCALE
        S = S + A(K,I)**2
  130    CONTINUE
!
     F = A(I,I)
     G = -SIGN(SQRT(S),F)
     H = F * G - S
     A(I,I) = F - G
     if (I  ==  N) go to 160
!
     DO 150 J = L, N
        S = 0.0E0
!
        DO 140 K = I, M
  140       S = S + A(K,I) * A(K,J)
!
        F = S / H
!
        DO 150 K = I, M
           A(K,J) = A(K,J) + F * A(K,I)
  150    CONTINUE
!
  160    if (IP  ==  0) go to 190
!
     DO 180 J = 1, IP
        S = 0.0E0
!
        DO 170 K = I, M
  170       S = S + A(K,I) * B(K,J)
!
        F = S / H
!
        DO 180 K = I, M
           B(K,J) = B(K,J) + F * A(K,I)
  180    CONTINUE
!
  190    DO 200 K = I, M
  200    A(K,I) = SCALE * A(K,I)
!
  210    W(I) = SCALE * G
     G = 0.0E0
     S = 0.0E0
     SCALE = 0.0E0
     if (I  >  M .OR. I  ==  N) go to 290
!
     DO 220 K = L, N
  220    SCALE = SCALE + ABS(A(I,K))
!
     if (SCALE  ==  0.0E0) go to 290
!
     DO 230 K = L, N
        A(I,K) = A(I,K) / SCALE
        S = S + A(I,K)**2
  230    CONTINUE
!
     F = A(I,L)
     G = -SIGN(SQRT(S),F)
     H = F * G - S
     A(I,L) = F - G
!
     DO 240 K = L, N
  240    RV1(K) = A(I,K) / H
!
     if (I  ==  M) go to 270
!
     DO 260 J = L, M
        S = 0.0E0
!
        DO 250 K = L, N
  250       S = S + A(J,K) * A(I,K)
!
        DO 260 K = L, N
           A(J,K) = A(J,K) + S * RV1(K)
  260    CONTINUE
!
  270    DO 280 K = L, N
  280    A(I,K) = SCALE * A(I,K)
!
  290    S1 = MAX(S1,ABS(W(I))+ABS(RV1(I)))
  300 CONTINUE
!     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS.
!                FOR I=N STEP -1 UNTIL 1 DO -- ..........
  DO 400 II = 1, N
     I = N + 1 - II
     if (I  ==  N) go to 390
     if (G  ==  0.0E0) go to 360
!
     DO 320 J = L, N
!     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
  320    A(J,I) = (A(I,J) / A(I,L)) / G
!
     DO 350 J = L, N
        S = 0.0E0
!
        DO 340 K = L, N
  340       S = S + A(I,K) * A(K,J)
!
        DO 350 K = L, N
           A(K,J) = A(K,J) + S * A(K,I)
  350    CONTINUE
!
  360    DO 380 J = L, N
        A(I,J) = 0.0E0
        A(J,I) = 0.0E0
  380    CONTINUE
!
  390    A(I,I) = 1.0E0
     G = RV1(I)
     L = I
  400 CONTINUE
!
  if (M  >=  N .OR. IP  ==  0) go to 510
  M1 = M + 1
!
  DO 500 I = M1, N
!
     DO 500 J = 1, IP
        B(I,J) = 0.0E0
  500 CONTINUE
!     .......... DIAGONALIZATION OF THE BIDIAGONAL FORM ..........
  510 CONTINUE
!     .......... FOR K=N STEP -1 UNTIL 1 DO -- ..........
  DO 700 KK = 1, N
     K1 = N - KK
     K = K1 + 1
     ITS = 0
!     .......... TEST FOR SPLITTING.
!                FOR L=K STEP -1 UNTIL 1 DO -- ..........
  520    DO 530 LL = 1, K
        L1 = K - LL
        L = L1 + 1
        if (S1 + ABS(RV1(L))  ==  S1) go to 565
!     .......... RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
        if (S1 + ABS(W(L1))  ==  S1) go to 540
  530    CONTINUE
!     .......... CANCELLATION OF RV1(L) if L GREATER THAN 1 ..........
  540    C = 0.0E0
     S = 1.0E0
!
     DO 560 I = L, K
        F = S * RV1(I)
        RV1(I) = C * RV1(I)
        if (S1 + ABS(F)  ==  S1) go to 565
        G = W(I)
        H = PYTHAG(F,G)
        W(I) = H
        C = G / H
        S = -F / H
        if (IP  ==  0) go to 560
!
        DO 550 J = 1, IP
           Y = B(L1,J)
           Z = B(I,J)
           B(L1,J) = Y * C + Z * S
           B(I,J) = -Y * S + Z * C
  550       CONTINUE
!
  560    CONTINUE
!     .......... TEST FOR CONVERGENCE ..........
  565    Z = W(K)
     if (L  ==  K) go to 650
!     .......... SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
     if (ITS  ==  30) go to 1000
     ITS = ITS + 1
     X = W(L)
     Y = W(K1)
     G = RV1(K1)
     H = RV1(K)
     F = 0.5E0 * (((G + Z) / H) * ((G - Z) / Y) + Y / H - H / Y)
     G = PYTHAG(F,1.0E0)
     F = X - (Z / X) * Z + (H / X) * (Y / (F + SIGN(G,F)) - H)
!     .......... NEXT QR TRANSFORMATION ..........
     C = 1.0E0
     S = 1.0E0
!
     DO 600 I1 = L, K1
        I = I1 + 1
        G = RV1(I)
        Y = W(I)
        H = S * G
        G = C * G
        Z = PYTHAG(F,H)
        RV1(I1) = Z
        C = F / Z
        S = H / Z
        F = X * C + G * S
        G = -X * S + G * C
        H = Y * S
        Y = Y * C
!
        DO 570 J = 1, N
           X = A(J,I1)
           Z = A(J,I)
           A(J,I1) = X * C + Z * S
           A(J,I) = -X * S + Z * C
  570       CONTINUE
!
        Z = PYTHAG(F,H)
        W(I1) = Z
!     .......... ROTATION CAN BE ARBITRARY if Z IS ZERO ..........
        if (Z  ==  0.0E0) go to 580
        C = F / Z
        S = H / Z
  580       F = C * G + S * Y
        X = -S * G + C * Y
        if (IP  ==  0) go to 600
!
        DO 590 J = 1, IP
           Y = B(I1,J)
           Z = B(I,J)
           B(I1,J) = Y * C + Z * S
           B(I,J) = -Y * S + Z * C
  590       CONTINUE
!
  600    CONTINUE
!
     RV1(L) = 0.0E0
     RV1(K) = F
     W(K) = X
     go to 520
!     .......... CONVERGENCE ..........
  650    if (Z  >=  0.0E0) go to 700
!     .......... W(K) IS MADE NON-NEGATIVE ..........
     W(K) = -Z
!
     DO 690 J = 1, N
  690    A(J,K) = -A(J,K)
!
  700 CONTINUE
!
  go to 1001
!     .......... SET ERROR -- NO CONVERGENCE TO A
!                SINGULAR VALUE AFTER 30 ITERATIONS ..........
 1000 IERR = K
 1001 RETURN
end
