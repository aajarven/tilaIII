subroutine CINVIT (NM, N, AR, AI, WR, WI, SELECT, MM, M, ZR, ZI, &
     IERR, RM1, RM2, RV1, RV2)
!
!! CINVIT computes the eigenvectors of a complex upper Hessenberg matrix ...
!  associated with specified eigenvalues using inverse iteration.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C2B
!***TYPE      COMPLEX (INVIT-S, CINVIT-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure CXINVIT
!     by Peters and Wilkinson.
!     HANDBOOK FOR AUTO. COMP. VOL.II-LINEAR ALGEBRA, 418-439(1971).
!
!     This subroutine finds those eigenvectors of A COMPLEX UPPER
!     Hessenberg matrix corresponding to specified eigenvalues,
!     using inverse iteration.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, AR, AI, ZR and ZI, as declared in the
!          calling program dimension statement.  NM is an INTEGER
!          variable.
!
!        N is the order of the matrix A=(AR,AI).  N is an INTEGER
!          variable.  N must be less than or equal to NM.
!
!        AR and AI contain the real and imaginary parts, respectively,
!          of the complex upper Hessenberg matrix.  AR and AI are
!          two-dimensional REAL arrays, dimensioned AR(NM,N)
!          and AI(NM,N).
!
!        WR and WI contain the real and imaginary parts, respectively,
!          of the eigenvalues of the matrix.  The eigenvalues must be
!          stored in a manner identical to that of subroutine  COMLR,
!          which recognizes possible splitting of the matrix.  WR and
!          WI are one-dimensional REAL arrays, dimensioned WR(N) and
!          WI(N).
!
!        SELECT specifies the eigenvectors to be found.  The
!          eigenvector corresponding to the J-th eigenvalue is
!          specified by setting SELECT(J) to .TRUE.  SELECT is a
!          one-dimensional LOGICAL array, dimensioned SELECT(N).
!
!        MM should be set to an upper bound for the number of
!          eigenvectors to be found.  MM is an INTEGER variable.
!
!     On OUTPUT
!
!        AR, AI, WI, and SELECT are unaltered.
!
!        WR may have been altered since close eigenvalues are perturbed
!          slightly in searching for independent eigenvectors.
!
!        M is the number of eigenvectors actually found.  M is an
!          INTEGER variable.
!
!        ZR and ZI contain the real and imaginary parts, respectively,
!          of the eigenvectors corresponding to the flagged eigenvalues.
!          The eigenvectors are normalized so that the component of
!          largest magnitude is 1.  Any vector which fails the
!          acceptance test is set to zero.  ZR and ZI are
!          two-dimensional REAL arrays, dimensioned ZR(NM,MM) and
!          ZI(NM,MM).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          -(2*N+1)   if more than MM eigenvectors have been requested
!                     (the MM eigenvectors calculated to this point are
!                     in ZR and ZI),
!          -K         if the iteration corresponding to the K-th
!                     value fails (if this occurs more than once, K
!                     is the index of the last occurrence); the
!                     corresponding columns of ZR and ZI are set to
!                     zero vectors,
!          -(N+K)     if both error situations occur.
!
!        RV1 and RV2 are one-dimensional REAL arrays used for
!          temporary storage, dimensioned RV1(N) and RV2(N).
!          They hold the approximate eigenvectors during the inverse
!          iteration process.
!
!        RM1 and RM2 are two-dimensional REAL arrays used for
!          temporary storage, dimensioned RM1(N,N) and RM2(N,N).
!          These arrays hold the triangularized form of the upper
!          Hessenberg matrix used in the inverse iteration process.
!
!     The ALGOL procedure GUESSVEC appears in CINVIT in-line.
!
!     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!     Calls CDIV for complex division.
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  CDIV, PYTHAG
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CINVIT
!
  INTEGER I,J,K,M,N,S,II,MM,MP,NM,UK,IP1,ITS,KM1,IERR
  REAL AR(NM,*),AI(NM,*),WR(*),WI(*),ZR(NM,*),ZI(NM,*)
  REAL RM1(N,*),RM2(N,*),RV1(*),RV2(*)
  REAL X,Y,EPS3,NORM,NORMV,GROWTO,ILAMBD,RLAMBD,UKROOT
  REAL PYTHAG
  LOGICAL SELECT(N)
!
!***FIRST EXECUTABLE STATEMENT  CINVIT
  IERR = 0
  UK = 0
  S = 1
!
  DO 980 K = 1, N
     if (.NOT. SELECT(K)) go to 980
     if (S  >  MM) go to 1000
     if (UK  >=  K) go to 200
!     .......... CHECK FOR POSSIBLE SPLITTING ..........
     DO 120 UK = K, N
        if (UK  ==  N) go to 140
        if (AR(UK+1,UK)  ==  0.0E0 .AND. AI(UK+1,UK)  ==  0.0E0) &
           go to 140
  120    CONTINUE
!     .......... COMPUTE INFINITY NORM OF LEADING UK BY UK
!                (HESSENBERG) MATRIX ..........
  140    NORM = 0.0E0
     MP = 1
!
     DO 180 I = 1, UK
        X = 0.0E0
!
        DO 160 J = MP, UK
  160       X = X + PYTHAG(AR(I,J),AI(I,J))
!
        if (X  >  NORM) NORM = X
        MP = I
  180    CONTINUE
!     .......... EPS3 REPLACES ZERO PIVOT IN DECOMPOSITION
!                AND CLOSE ROOTS ARE MODIFIED BY EPS3 ..........
     if (NORM  ==  0.0E0) NORM = 1.0E0
     EPS3 = NORM
  190    EPS3 = 0.5E0*EPS3
     if (NORM + EPS3  >  NORM) go to 190
     EPS3 = 2.0E0*EPS3
!     .......... GROWTO IS THE CRITERION FOR GROWTH ..........
     UKROOT = SQRT(REAL(UK))
     GROWTO = 0.1E0 / UKROOT
  200    RLAMBD = WR(K)
     ILAMBD = WI(K)
     if (K  ==  1) go to 280
     KM1 = K - 1
     go to 240
!     .......... PERTURB EIGENVALUE if IT IS CLOSE
!                TO ANY PREVIOUS EIGENVALUE ..........
  220    RLAMBD = RLAMBD + EPS3
!     .......... FOR I=K-1 STEP -1 UNTIL 1 DO -- ..........
  240    DO 260 II = 1, KM1
        I = K - II
        if (SELECT(I) .AND. ABS(WR(I)-RLAMBD)  <  EPS3 .AND. &
           ABS(WI(I)-ILAMBD)  <  EPS3) go to 220
  260    CONTINUE
!
     WR(K) = RLAMBD
!     .......... FORM UPPER HESSENBERG (AR,AI)-(RLAMBD,ILAMBD)*I
!                AND INITIAL COMPLEX VECTOR ..........
  280    MP = 1
!
     DO 320 I = 1, UK
!
        DO 300 J = MP, UK
           RM1(I,J) = AR(I,J)
           RM2(I,J) = AI(I,J)
  300       CONTINUE
!
        RM1(I,I) = RM1(I,I) - RLAMBD
        RM2(I,I) = RM2(I,I) - ILAMBD
        MP = I
        RV1(I) = EPS3
  320    CONTINUE
!     .......... TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
!                REPLACING ZERO PIVOTS BY EPS3 ..........
     if (UK  ==  1) go to 420
!
     DO 400 I = 2, UK
        MP = I - 1
        if (PYTHAG(RM1(I,MP),RM2(I,MP))  <=  &
           PYTHAG(RM1(MP,MP),RM2(MP,MP))) go to 360
!
        DO 340 J = MP, UK
           Y = RM1(I,J)
           RM1(I,J) = RM1(MP,J)
           RM1(MP,J) = Y
           Y = RM2(I,J)
           RM2(I,J) = RM2(MP,J)
           RM2(MP,J) = Y
  340       CONTINUE
!
  360       if (RM1(MP,MP)  ==  0.0E0 .AND. RM2(MP,MP)  ==  0.0E0) &
           RM1(MP,MP) = EPS3
        call CDIV(RM1(I,MP),RM2(I,MP),RM1(MP,MP),RM2(MP,MP),X,Y)
        if (X  ==  0.0E0 .AND. Y  ==  0.0E0) go to 400
!
        DO 380 J = I, UK
           RM1(I,J) = RM1(I,J) - X * RM1(MP,J) + Y * RM2(MP,J)
           RM2(I,J) = RM2(I,J) - X * RM2(MP,J) - Y * RM1(MP,J)
  380       CONTINUE
!
  400    CONTINUE
!
  420    if (RM1(UK,UK)  ==  0.0E0 .AND. RM2(UK,UK)  ==  0.0E0) &
        RM1(UK,UK) = EPS3
     ITS = 0
!     .......... BACK SUBSTITUTION
!                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  660    DO 720 II = 1, UK
        I = UK + 1 - II
        X = RV1(I)
        Y = 0.0E0
        if (I  ==  UK) go to 700
        IP1 = I + 1
!
        DO 680 J = IP1, UK
           X = X - RM1(I,J) * RV1(J) + RM2(I,J) * RV2(J)
           Y = Y - RM1(I,J) * RV2(J) - RM2(I,J) * RV1(J)
  680       CONTINUE
!
  700       call CDIV(X,Y,RM1(I,I),RM2(I,I),RV1(I),RV2(I))
  720    CONTINUE
!     .......... ACCEPTANCE TEST FOR EIGENVECTOR
!                AND NORMALIZATION ..........
     ITS = ITS + 1
     NORM = 0.0E0
     NORMV = 0.0E0
!
     DO 780 I = 1, UK
        X = PYTHAG(RV1(I),RV2(I))
        if (NORMV  >=  X) go to 760
        NORMV = X
        J = I
  760       NORM = NORM + X
  780    CONTINUE
!
     if (NORM  <  GROWTO) go to 840
!     .......... ACCEPT VECTOR ..........
     X = RV1(J)
     Y = RV2(J)
!
     DO 820 I = 1, UK
        call CDIV(RV1(I),RV2(I),X,Y,ZR(I,S),ZI(I,S))
  820    CONTINUE
!
     if (UK  ==  N) go to 940
     J = UK + 1
     go to 900
!     .......... IN-LINE PROCEDURE FOR CHOOSING
!                A NEW STARTING VECTOR ..........
  840    if (ITS  >=  UK) go to 880
     X = UKROOT
     Y = EPS3 / (X + 1.0E0)
     RV1(1) = EPS3
!
     DO 860 I = 2, UK
  860    RV1(I) = Y
!
     J = UK - ITS + 1
     RV1(J) = RV1(J) - EPS3 * X
     go to 660
!     .......... SET ERROR -- UNACCEPTED EIGENVECTOR ..........
  880    J = 1
     IERR = -K
!     .......... SET REMAINING VECTOR COMPONENTS TO ZERO ..........
  900    DO 920 I = J, N
        ZR(I,S) = 0.0E0
        ZI(I,S) = 0.0E0
  920    CONTINUE
!
  940    S = S + 1
  980 CONTINUE
!
  go to 1001
!     .......... SET ERROR -- UNDERESTIMATE OF EIGENVECTOR
!                SPACE REQUIRED ..........
 1000 if (IERR  /=  0) IERR = IERR - N
  if (IERR  ==  0) IERR = -(2 * N + 1)
 1001 M = S - 1
  return
end
