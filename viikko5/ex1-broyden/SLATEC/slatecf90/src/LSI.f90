subroutine LSI (W, MDW, MA, MG, N, PRGOPT, X, RNORM, MODE, WS, IP)
!
!! LSI is subsidiary to LSEI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (LSI-S, DLSI-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!     This is a companion subprogram to LSEI.  The documentation for
!     LSEI has complete usage instructions.
!
!     Solve..
!              AX = B,  A  MA by N  (least squares equations)
!     subject to..
!
!              GX >= H, G  MG by N  (inequality constraints)
!
!     Input..
!
!      W(*,*) contains  (A B) in rows 1,...,MA+MG, cols 1,...,N+1.
!                       (G H)
!
!     MDW,MA,MG,N
!              contain (resp) var. dimension of W(*,*),
!              and matrix dimensions.
!
!     PRGOPT(*),
!              Program option vector.
!
!     OUTPUT..
!
!      X(*),RNORM
!
!              Solution vector(unless MODE=2), length of AX-B.
!
!      MODE
!              =0   Inequality constraints are compatible.
!              =2   Inequality constraints contradictory.
!
!      WS(*),
!              Working storage of dimension K+N+(MG+2)*(N+7),
!              where K=MAX(MA+MG,N).
!      IP(MG+2*N+1)
!              Integer working storage
!
!***ROUTINES CALLED  H12, HFTI, LPDP, R1MACH, SASUM, SAXPY, SCOPY, SDOT,
!                    SSCAL, SSWAP
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890618  Completely restructured and extensively revised (WRB & RWC)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   920422  Changed call to HFTI to include variable MA.  (WRB)
!***END PROLOGUE  LSI
  INTEGER IP(*), MA, MDW, MG, MODE, N
  REAL             PRGOPT(*), RNORM, W(MDW,*), WS(*), X(*)
!
  EXTERNAL H12, HFTI, LPDP, R1MACH, SASUM, SAXPY, SCOPY, SDOT, &
     SSCAL, SSWAP
  REAL             R1MACH, SASUM, SDOT
!
  REAL             ANORM, FAC, GAM, RB, SRELPR, TAU, TOL, XNORM
  INTEGER I, J, K, KEY, KRANK, KRM1, KRP1, L, LAST, LINK, M, MAP1, &
     MDLPDP, MINMAN, N1, N2, N3, NEXT, NP1
  LOGICAL COV, FIRST, SCLCOV
!
  SAVE SRELPR, FIRST
  DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  LSI
!
!     Set the nominal tolerance used in the code.
!
  if (FIRST) SRELPR = R1MACH(4)
  FIRST = .FALSE.
  TOL = SQRT(SRELPR)
!
  MODE = 0
  RNORM = 0.E0
  M = MA + MG
  NP1 = N + 1
  KRANK = 0
  if (N <= 0 .OR. M <= 0) go to 370
!
!     To process option vector.
!
  COV = .FALSE.
  SCLCOV = .TRUE.
  LAST = 1
  LINK = PRGOPT(1)
!
  100 if (LINK > 1) THEN
     KEY = PRGOPT(LAST+1)
     if (KEY == 1) COV = PRGOPT(LAST+2)  /=  0.E0
     if (KEY == 10) SCLCOV = PRGOPT(LAST+2)  ==  0.E0
     if (KEY == 5) TOL = MAX(SRELPR,PRGOPT(LAST+2))
     NEXT = PRGOPT(LINK)
     LAST = LINK
     LINK = NEXT
     go to 100
  end if
!
!     Compute matrix norm of least squares equations.
!
  ANORM = 0.E0
  DO 110 J = 1,N
     ANORM = MAX(ANORM,SASUM(MA,W(1,J),1))
  110 CONTINUE
!
!     Set tolerance for HFTI( ) rank test.
!
  TAU = TOL*ANORM
!
!     Compute Householder orthogonal decomposition of matrix.
!
  call SCOPY (N, 0.E0, 0, WS, 1)
  call SCOPY (MA, W(1, NP1), 1, WS, 1)
  K = MAX(M,N)
  MINMAN = MIN(MA,N)
  N1 = K + 1
  N2 = N1 + N
  call HFTI (W, MDW, MA, N, WS, MA, 1, TAU, KRANK, RNORM, WS(N2), &
             WS(N1), IP)
  FAC = 1.E0
  GAM = MA - KRANK
  if (KRANK < MA .AND. SCLCOV) FAC = RNORM**2/GAM
!
!     Reduce to LPDP and solve.
!
  MAP1 = MA + 1
!
!     Compute inequality rt-hand side for LPDP.
!
  if (MA < M) THEN
     if (MINMAN > 0) THEN
        DO 120 I = MAP1,M
           W(I,NP1) = W(I,NP1) - SDOT(N,W(I,1),MDW,WS,1)
  120       CONTINUE
!
!           Apply permutations to col. of inequality constraint matrix.
!
        DO 130 I = 1,MINMAN
           call SSWAP (MG, W(MAP1,I), 1, W(MAP1,IP(I)), 1)
  130       CONTINUE
!
!           Apply Householder transformations to constraint matrix.
!
        if (KRANK > 0 .AND. KRANK < N) THEN
           DO 140 I = KRANK,1,-1
              call H12 (2, I, KRANK+1, N, W(I,1), MDW, WS(N1+I-1), &
                        W(MAP1,1), MDW, 1, MG)
  140          CONTINUE
        ENDIF
!
!           Compute permuted inequality constraint matrix times r-inv.
!
        DO 160 I = MAP1,M
           DO 150 J = 1,KRANK
              W(I,J) = (W(I,J)-SDOT(J-1,W(1,J),1,W(I,1),MDW))/W(J,J)
  150          CONTINUE
  160       CONTINUE
     ENDIF
!
!        Solve the reduced problem with LPDP algorithm,
!        the least projected distance problem.
!
     call LPDP(W(MAP1,1), MDW, MG, KRANK, N-KRANK, PRGOPT, X, &
               XNORM, MDLPDP, WS(N2), IP(N+1))
!
!        Compute solution in original coordinates.
!
     if (MDLPDP == 1) THEN
        DO 170 I = KRANK,1,-1
           X(I) = (X(I)-SDOT(KRANK-I,W(I,I+1),MDW,X(I+1),1))/W(I,I)
  170       CONTINUE
!
!           Apply Householder transformation to solution vector.
!
        if (KRANK < N) THEN
           DO 180 I = 1,KRANK
              call H12 (2, I, KRANK+1, N, W(I,1), MDW, WS(N1+I-1), &
                        X, 1, 1, 1)
  180          CONTINUE
        ENDIF
!
!           Repermute variables to their input order.
!
        if (MINMAN > 0) THEN
           DO 190 I = MINMAN,1,-1
              call SSWAP (1, X(I), 1, X(IP(I)), 1)
  190          CONTINUE
!
!              Variables are now in original coordinates.
!              Add solution of unconstrained problem.
!
           DO 200 I = 1,N
              X(I) = X(I) + WS(I)
  200          CONTINUE
!
!              Compute the residual vector norm.
!
           RNORM = SQRT(RNORM**2+XNORM**2)
        ENDIF
     ELSE
        MODE = 2
     ENDIF
  ELSE
     call SCOPY (N, WS, 1, X, 1)
  end if
!
!     Compute covariance matrix based on the orthogonal decomposition
!     from HFTI( ).
!
  if (.NOT.COV .OR. KRANK <= 0) go to 370
  KRM1 = KRANK - 1
  KRP1 = KRANK + 1
!
!     Copy diagonal terms to working array.
!
  call SCOPY (KRANK, W, MDW+1, WS(N2), 1)
!
!     Reciprocate diagonal terms.
!
  DO 210 J = 1,KRANK
     W(J,J) = 1.E0/W(J,J)
  210 CONTINUE
!
!     Invert the upper triangular QR factor on itself.
!
  if (KRANK > 1) THEN
     DO 230 I = 1,KRM1
        DO 220 J = I+1,KRANK
           W(I,J) = -SDOT(J-I,W(I,I),MDW,W(I,J),1)*W(J,J)
  220       CONTINUE
  230    CONTINUE
  end if
!
!     Compute the inverted factor times its transpose.
!
  DO 250 I = 1,KRANK
     DO 240 J = I,KRANK
        W(I,J) = SDOT(KRANK+1-J,W(I,J),MDW,W(J,J),MDW)
  240    CONTINUE
  250 CONTINUE
!
!     Zero out lower trapezoidal part.
!     Copy upper triangular to lower triangular part.
!
  if (KRANK < N) THEN
     DO 260 J = 1,KRANK
        call SCOPY (J, W(1,J), 1, W(J,1), MDW)
  260    CONTINUE
!
     DO 270 I = KRP1,N
        call SCOPY (I, 0.E0, 0, W(I,1), MDW)
  270    CONTINUE
!
!        Apply right side transformations to lower triangle.
!
     N3 = N2 + KRP1
     DO 330 I = 1,KRANK
        L = N1 + I
        K = N2 + I
        RB = WS(L-1)*WS(K-1)
!
!           If RB >= 0.E0, transformation can be regarded as zero.
!
        if (RB < 0.E0) THEN
           RB = 1.E0/RB
!
!              Store unscaled rank one Householder update in work array.
!
           call SCOPY (N, 0.E0, 0, WS(N3), 1)
           L = N1 + I
           K = N3 + I
           WS(K-1) = WS(L-1)
!
           DO 280 J = KRP1,N
              WS(N3+J-1) = W(I,J)
  280          CONTINUE
!
           DO 290 J = 1,N
              WS(J) = RB*(SDOT(J-I,W(J,I),MDW,WS(N3+I-1),1)+ &
                      SDOT(N-J+1,W(J,J),1,WS(N3+J-1),1))
  290          CONTINUE
!
           L = N3 + I
           GAM = 0.5E0*RB*SDOT(N-I+1,WS(L-1),1,WS(I),1)
           call SAXPY (N-I+1, GAM, WS(L-1), 1, WS(I), 1)
           DO 320 J = I,N
              DO 300 L = 1,I-1
                 W(J,L) = W(J,L) + WS(N3+J-1)*WS(L)
  300             CONTINUE
!
              DO 310 L = I,J
                 W(J,L) = W(J,L) + WS(J)*WS(N3+L-1)+WS(L)*WS(N3+J-1)
  310             CONTINUE
  320          CONTINUE
        ENDIF
  330    CONTINUE
!
!        Copy lower triangle to upper triangle to symmetrize the
!        covariance matrix.
!
     DO 340 I = 1,N
        call SCOPY (I, W(I,1), MDW, W(1,I), 1)
  340    CONTINUE
  end if
!
!     Repermute rows and columns.
!
  DO 350 I = MINMAN,1,-1
     K = IP(I)
     if (I /= K) THEN
        call SSWAP (1, W(I,I), 1, W(K,K), 1)
        call SSWAP (I-1, W(1,I), 1, W(1,K), 1)
        call SSWAP (K-I-1, W(I,I+1), MDW, W(I+1,K), 1)
        call SSWAP (N-K, W(I, K+1), MDW, W(K, K+1), MDW)
     ENDIF
  350 CONTINUE
!
!     Put in normalized residual sum of squares scale factor
!     and symmetrize the resulting covariance matrix.
!
  DO 360 J = 1,N
     call SSCAL (J, FAC, W(1,J), 1)
     call SCOPY (J, W(1,J), 1, W(J,1), MDW)
  360 CONTINUE
!
  370 IP(1) = KRANK
  IP(2) = N + MAX(M,N) + (MG+2)*(N+7)
  return
end
