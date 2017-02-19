subroutine SVD (NM, M, N, A, W, MATU, U, MATV, V, IERR, RV1)
!
!! SVD performs the singular value decomposition of a rectangular matrix.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (SVD-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure SVD,
!     NUM. MATH. 14, 403-420(1970) by Golub and Reinsch.
!     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
!
!     This subroutine determines the singular value decomposition
!          T
!     A=USV  of a REAL M by N rectangular matrix.  Householder
!     bidiagonalization and a variant of the QR algorithm are used.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A, U and V, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!          Note that NM must be at least as large as the maximum
!          of M and N.
!
!        M is the number of rows of A and U.
!
!        N is the number of columns of A and U and the order of V.
!
!        A contains the rectangular input matrix to be decomposed.  A is
!          a two-dimensional REAL array, dimensioned A(NM,N).
!
!        MATU should be set to .TRUE. if the U matrix in the
!          decomposition is desired, and to .FALSE. otherwise.
!          MATU is a LOGICAL variable.
!
!        MATV should be set to .TRUE. if the V matrix in the
!          decomposition is desired, and to .FALSE. otherwise.
!          MATV is a LOGICAL variable.
!
!     On Output
!
!        A is unaltered (unless overwritten by U or V).
!
!        W contains the N (non-negative) singular values of A (the
!          diagonal elements of S).  They are unordered.  If an
!          error exit is made, the singular values should be correct
!          for indices IERR+1, IERR+2, ..., N.  W is a one-dimensional
!          REAL array, dimensioned W(N).
!
!        U contains the matrix U (orthogonal column vectors) of the
!          decomposition if MATU has been set to .TRUE.  Otherwise,
!          U is used as a temporary array.  U may coincide with A.
!          If an error exit is made, the columns of U corresponding
!          to indices of correct singular values should be correct.
!          U is a two-dimensional REAL array, dimensioned U(NM,N).
!
!        V contains the matrix V (orthogonal) of the decomposition if
!          MATV has been set to .TRUE.  Otherwise, V is not referenced.
!          V may also coincide with A if U does not.  If an error
!          exit is made, the columns of V corresponding to indices of
!          correct singular values should be correct.  V is a two-
!          dimensional REAL array, dimensioned V(NM,N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          K          if the K-th singular value has not been
!                     determined after 30 iterations.
!
!        RV1 is a one-dimensional REAL array used for temporary storage,
!          dimensioned RV1(N).
!
!     CALLS PYTHAG(A,B) for sqrt(A**2 + B**2).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***SEE ALSO  EISDOC
!***ROUTINES CALLED  PYTHAG
!***REVISION HISTORY  (YYMMDD)
!   811101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  SVD
!
  INTEGER I,J,K,L,M,N,II,I1,KK,K1,LL,L1,MN,NM,ITS,IERR
  REAL A(NM,*),W(*),U(NM,*),V(NM,*),RV1(*)
  REAL C,F,G,H,S,X,Y,Z,SCALE,S1
  REAL PYTHAG
  LOGICAL MATU,MATV
!
!***FIRST EXECUTABLE STATEMENT  SVD
  IERR = 0
!
  DO 100 I = 1, M
!
     DO 100 J = 1, N
        U(I,J) = A(I,J)
  100 CONTINUE
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
  120    SCALE = SCALE + ABS(U(K,I))
!
     if (SCALE  ==  0.0E0) go to 210
!
     DO 130 K = I, M
        U(K,I) = U(K,I) / SCALE
        S = S + U(K,I)**2
  130    CONTINUE
!
     F = U(I,I)
     G = -SIGN(SQRT(S),F)
     H = F * G - S
     U(I,I) = F - G
     if (I  ==  N) go to 190
!
     DO 150 J = L, N
        S = 0.0E0
!
        DO 140 K = I, M
  140       S = S + U(K,I) * U(K,J)
!
        F = S / H
!
        DO 150 K = I, M
           U(K,J) = U(K,J) + F * U(K,I)
  150    CONTINUE
!
  190    DO 200 K = I, M
  200    U(K,I) = SCALE * U(K,I)
!
  210    W(I) = SCALE * G
     G = 0.0E0
     S = 0.0E0
     SCALE = 0.0E0
     if (I  >  M .OR. I  ==  N) go to 290
!
     DO 220 K = L, N
  220    SCALE = SCALE + ABS(U(I,K))
!
     if (SCALE  ==  0.0E0) go to 290
!
     DO 230 K = L, N
        U(I,K) = U(I,K) / SCALE
        S = S + U(I,K)**2
  230    CONTINUE
!
     F = U(I,L)
     G = -SIGN(SQRT(S),F)
     H = F * G - S
     U(I,L) = F - G
!
     DO 240 K = L, N
  240    RV1(K) = U(I,K) / H
!
     if (I  ==  M) go to 270
!
     DO 260 J = L, M
        S = 0.0E0
!
        DO 250 K = L, N
  250       S = S + U(J,K) * U(I,K)
!
        DO 260 K = L, N
           U(J,K) = U(J,K) + S * RV1(K)
  260    CONTINUE
!
  270    DO 280 K = L, N
  280    U(I,K) = SCALE * U(I,K)
!
  290    S1 = MAX(S1,ABS(W(I))+ABS(RV1(I)))
  300 CONTINUE
!     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS ..........
  if (.NOT. MATV) go to 410
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
  DO 400 II = 1, N
     I = N + 1 - II
     if (I  ==  N) go to 390
     if (G  ==  0.0E0) go to 360
!
     DO 320 J = L, N
!     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
  320    V(J,I) = (U(I,J) / U(I,L)) / G
!
     DO 350 J = L, N
        S = 0.0E0
!
        DO 340 K = L, N
  340       S = S + U(I,K) * V(K,J)
!
        DO 350 K = L, N
           V(K,J) = V(K,J) + S * V(K,I)
  350    CONTINUE
!
  360    DO 380 J = L, N
        V(I,J) = 0.0E0
        V(J,I) = 0.0E0
  380    CONTINUE
!
  390    V(I,I) = 1.0E0
     G = RV1(I)
     L = I
  400 CONTINUE
!     .......... ACCUMULATION OF LEFT-HAND TRANSFORMATIONS ..........
  410 if (.NOT. MATU) go to 510
!     ..........FOR I=MIN(M,N) STEP -1 UNTIL 1 DO -- ..........
  MN = N
  if (M  <  N) MN = M
!
  DO 500 II = 1, MN
     I = MN + 1 - II
     L = I + 1
     G = W(I)
     if (I  ==  N) go to 430
!
     DO 420 J = L, N
  420    U(I,J) = 0.0E0
!
  430    if (G  ==  0.0E0) go to 475
     if (I  ==  MN) go to 460
!
     DO 450 J = L, N
        S = 0.0E0
!
        DO 440 K = L, M
  440       S = S + U(K,I) * U(K,J)
!     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
        F = (S / U(I,I)) / G
!
        DO 450 K = I, M
           U(K,J) = U(K,J) + F * U(K,I)
  450    CONTINUE
!
  460    DO 470 J = I, M
  470    U(J,I) = U(J,I) / G
!
     go to 490
!
  475    DO 480 J = I, M
  480    U(J,I) = 0.0E0
!
  490    U(I,I) = U(I,I) + 1.0E0
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
        if (.NOT. MATU) go to 560
!
        DO 550 J = 1, M
           Y = U(J,L1)
           Z = U(J,I)
           U(J,L1) = Y * C + Z * S
           U(J,I) = -Y * S + Z * C
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
        if (.NOT. MATV) go to 575
!
        DO 570 J = 1, N
           X = V(J,I1)
           Z = V(J,I)
           V(J,I1) = X * C + Z * S
           V(J,I) = -X * S + Z * C
  570       CONTINUE
!
  575       Z = PYTHAG(F,H)
        W(I1) = Z
!     .......... ROTATION CAN BE ARBITRARY if Z IS ZERO ..........
        if (Z  ==  0.0E0) go to 580
        C = F / Z
        S = H / Z
  580       F = C * G + S * Y
        X = -S * G + C * Y
        if (.NOT. MATU) go to 600
!
        DO 590 J = 1, M
           Y = U(J,I1)
           Z = U(J,I)
           U(J,I1) = Y * C + Z * S
           U(J,I) = -Y * S + Z * C
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
     if (.NOT. MATV) go to 700
!
     DO 690 J = 1, N
  690    V(J,K) = -V(J,K)
!
  700 CONTINUE
!
  go to 1001
!     .......... SET ERROR -- NO CONVERGENCE TO A
!                SINGULAR VALUE AFTER 30 ITERATIONS ..........
 1000 IERR = K
 1001 RETURN
end
