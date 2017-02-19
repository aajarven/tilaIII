subroutine INVIT (NM, N, A, WR, WI, SELECT, MM, M, Z, IERR, RM1, &
     RV1, RV2)
!
!! INVIT computes the eigenvectors of a real upper Hessenberg ...
!            matrix associated with specified eigenvalues by inverse ...
!            iteration.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C2B
!***TYPE      SINGLE PRECISION (INVIT-S, CINVIT-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure INVIT
!     by Peters and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
!
!     This subroutine finds those eigenvectors of a REAL UPPER
!     Hessenberg matrix corresponding to specified eigenvalues,
!     using inverse iteration.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        A contains the upper Hessenberg matrix.  A is a two-dimensional
!          REAL array, dimensioned A(NM,N).
!
!        WR and WI contain the real and imaginary parts, respectively,
!          of the eigenvalues of the Hessenberg matrix.  The eigenvalues
!          must be stored in a manner identical to that output by
!          subroutine  HQR,  which recognizes possible splitting of the
!          matrix.  WR and WI are one-dimensional REAL arrays,
!          dimensioned WR(N) and WI(N).
!
!        SELECT specifies the eigenvectors to be found. The
!          eigenvector corresponding to the J-th eigenvalue is
!          specified by setting SELECT(J) to .TRUE.  SELECT is a
!          one-dimensional LOGICAL array, dimensioned SELECT(N).
!
!        MM should be set to an upper bound for the number of
!          columns required to store the eigenvectors to be found.
!          NOTE that two columns are required to store the
!          eigenvector corresponding to a complex eigenvalue.  One
!          column is required to store the eigenvector corresponding
!          to a real eigenvalue.  MM is an INTEGER variable.
!
!     On OUTPUT
!
!        A and WI are unaltered.
!
!        WR may have been altered since close eigenvalues are perturbed
!          slightly in searching for independent eigenvectors.
!
!        SELECT may have been altered.  If the elements corresponding
!          to a pair of conjugate complex eigenvalues were each
!          initially set to .TRUE., the program resets the second of
!          the two elements to .FALSE.
!
!        M is the number of columns actually used to store the
!          eigenvectors.  M is an INTEGER variable.
!
!        Z contains the real and imaginary parts of the eigenvectors.
!          The eigenvectors are packed into the columns of Z starting
!          at the first column.  If the next selected eigenvalue is
!          real, the next column of Z contains its eigenvector.  If the
!          eigenvalue is complex, the next two columns of Z contain the
!          real and imaginary parts of its eigenvector, with the real
!          part first.  The eigenvectors are normalized so that the
!          component of largest magnitude is 1. Any vector which fails
!          the acceptance test is set to zero.  Z is a two-dimensional
!          REAL array, dimensioned Z(NM,MM).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          -(2*N+1)   if more than MM columns of Z are necessary
!                     to store the eigenvectors corresponding to
!                     the specified eigenvalues (in this case, M is
!                     equal to the number of columns of Z containing
!                     eigenvectors already computed),
!          -K         if the iteration corresponding to the K-th
!                     value fails (if this occurs more than once, K
!                     is the index of the last occurrence); the
!                     corresponding columns of Z are set to zero
!                     vectors,
!          -(N+K)     if both error situations occur.
!
!        RM1 is a two-dimensional REAL array used for temporary storage.
!          This array holds the triangularized form of the upper
!          Hessenberg matrix used in the inverse iteration process.
!          RM1 is dimensioned RM1(N,N).
!
!        RV1 and RV2 are one-dimensional REAL arrays used for temporary
!          storage.  They hold the approximate eigenvectors during the
!          inverse iteration process.  RV1 and RV2 are dimensioned
!          RV1(N) and RV2(N).
!
!     The ALGOL procedure GUESSVEC appears in INVIT in-line.
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
!***END PROLOGUE  INVIT
!
  INTEGER I,J,K,L,M,N,S,II,IP,MM,MP,NM,NS,N1,UK,IP1,ITS,KM1,IERR
  REAL A(NM,*),WR(*),WI(*),Z(NM,*)
  REAL RM1(N,*),RV1(*),RV2(*)
  REAL T,W,X,Y,EPS3
  REAL NORM,NORMV,GROWTO,ILAMBD,RLAMBD,UKROOT
  REAL PYTHAG
  LOGICAL SELECT(N)
!
!***FIRST EXECUTABLE STATEMENT  INVIT
  IERR = 0
  UK = 0
  S = 1
!     .......... IP = 0, REAL EIGENVALUE
!                     1, FIRST OF CONJUGATE COMPLEX PAIR
!                    -1, SECOND OF CONJUGATE COMPLEX PAIR ..........
  IP = 0
  N1 = N - 1
!
  DO 980 K = 1, N
     if (WI(K)  ==  0.0E0 .OR. IP  <  0) go to 100
     IP = 1
     if (SELECT(K) .AND. SELECT(K+1)) SELECT(K+1) = .FALSE.
  100    if (.NOT. SELECT(K)) go to 960
     if (WI(K)  /=  0.0E0) S = S + 1
     if (S  >  MM) go to 1000
     if (UK  >=  K) go to 200
!     .......... CHECK FOR POSSIBLE SPLITTING ..........
     DO 120 UK = K, N
        if (UK  ==  N) go to 140
        if (A(UK+1,UK)  ==  0.0E0) go to 140
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
  160       X = X + ABS(A(I,J))
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
!     .......... GROWTO IS THE CRITERION FOR THE GROWTH ..........
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
!     .......... PERTURB CONJUGATE EIGENVALUE TO MATCH ..........
     IP1 = K + IP
     WR(IP1) = RLAMBD
!     .......... FORM UPPER HESSENBERG A-RLAMBD*I (TRANSPOSED)
!                AND INITIAL REAL VECTOR ..........
  280    MP = 1
!
     DO 320 I = 1, UK
!
        DO 300 J = MP, UK
  300       RM1(J,I) = A(I,J)
!
        RM1(I,I) = RM1(I,I) - RLAMBD
        MP = I
        RV1(I) = EPS3
  320    CONTINUE
!
     ITS = 0
     if (ILAMBD  /=  0.0E0) go to 520
!     .......... REAL EIGENVALUE.
!                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
!                REPLACING ZERO PIVOTS BY EPS3 ..........
     if (UK  ==  1) go to 420
!
     DO 400 I = 2, UK
        MP = I - 1
        if (ABS(RM1(MP,I))  <=  ABS(RM1(MP,MP))) go to 360
!
        DO 340 J = MP, UK
           Y = RM1(J,I)
           RM1(J,I) = RM1(J,MP)
           RM1(J,MP) = Y
  340       CONTINUE
!
  360       if (RM1(MP,MP)  ==  0.0E0) RM1(MP,MP) = EPS3
        X = RM1(MP,I) / RM1(MP,MP)
        if (X  ==  0.0E0) go to 400
!
        DO 380 J = I, UK
  380       RM1(J,I) = RM1(J,I) - X * RM1(J,MP)
!
  400    CONTINUE
!
  420    if (RM1(UK,UK)  ==  0.0E0) RM1(UK,UK) = EPS3
!     .......... BACK SUBSTITUTION FOR REAL VECTOR
!                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  440    DO 500 II = 1, UK
        I = UK + 1 - II
        Y = RV1(I)
        if (I  ==  UK) go to 480
        IP1 = I + 1
!
        DO 460 J = IP1, UK
  460       Y = Y - RM1(J,I) * RV1(J)
!
  480       RV1(I) = Y / RM1(I,I)
  500    CONTINUE
!
     go to 740
!     .......... COMPLEX EIGENVALUE.
!                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
!                REPLACING ZERO PIVOTS BY EPS3.  STORE IMAGINARY
!                PARTS IN UPPER TRIANGLE STARTING AT (1,3) ..........
  520    NS = N - S
     Z(1,S-1) = -ILAMBD
     Z(1,S) = 0.0E0
     if (N  ==  2) go to 550
     RM1(1,3) = -ILAMBD
     Z(1,S-1) = 0.0E0
     if (N  ==  3) go to 550
!
     DO 540 I = 4, N
  540    RM1(1,I) = 0.0E0
!
  550    DO 640 I = 2, UK
        MP = I - 1
        W = RM1(MP,I)
        if (I  <  N) T = RM1(MP,I+1)
        if (I  ==  N) T = Z(MP,S-1)
        X = RM1(MP,MP) * RM1(MP,MP) + T * T
        if (W * W  <=  X) go to 580
        X = RM1(MP,MP) / W
        Y = T / W
        RM1(MP,MP) = W
        if (I  <  N) RM1(MP,I+1) = 0.0E0
        if (I  ==  N) Z(MP,S-1) = 0.0E0
!
        DO 560 J = I, UK
           W = RM1(J,I)
           RM1(J,I) = RM1(J,MP) - X * W
           RM1(J,MP) = W
           if (J  <  N1) go to 555
           L = J - NS
           Z(I,L) = Z(MP,L) - Y * W
           Z(MP,L) = 0.0E0
           go to 560
  555          RM1(I,J+2) = RM1(MP,J+2) - Y * W
           RM1(MP,J+2) = 0.0E0
  560       CONTINUE
!
        RM1(I,I) = RM1(I,I) - Y * ILAMBD
        if (I  <  N1) go to 570
        L = I - NS
        Z(MP,L) = -ILAMBD
        Z(I,L) = Z(I,L) + X * ILAMBD
        go to 640
  570       RM1(MP,I+2) = -ILAMBD
        RM1(I,I+2) = RM1(I,I+2) + X * ILAMBD
        go to 640
  580       if (X  /=  0.0E0) go to 600
        RM1(MP,MP) = EPS3
        if (I  <  N) RM1(MP,I+1) = 0.0E0
        if (I  ==  N) Z(MP,S-1) = 0.0E0
        T = 0.0E0
        X = EPS3 * EPS3
  600       W = W / X
        X = RM1(MP,MP) * W
        Y = -T * W
!
        DO 620 J = I, UK
           if (J  <  N1) go to 610
           L = J - NS
           T = Z(MP,L)
           Z(I,L) = -X * T - Y * RM1(J,MP)
           go to 615
  610          T = RM1(MP,J+2)
           RM1(I,J+2) = -X * T - Y * RM1(J,MP)
  615          RM1(J,I) = RM1(J,I) - X * RM1(J,MP) + Y * T
  620       CONTINUE
!
        if (I  <  N1) go to 630
        L = I - NS
        Z(I,L) = Z(I,L) - ILAMBD
        go to 640
  630       RM1(I,I+2) = RM1(I,I+2) - ILAMBD
  640    CONTINUE
!
     if (UK  <  N1) go to 650
     L = UK - NS
     T = Z(UK,L)
     go to 655
  650    T = RM1(UK,UK+2)
  655    if (RM1(UK,UK)  ==  0.0E0 .AND. T  ==  0.0E0) RM1(UK,UK) = EPS3
!     .......... BACK SUBSTITUTION FOR COMPLEX VECTOR
!                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  660    DO 720 II = 1, UK
        I = UK + 1 - II
        X = RV1(I)
        Y = 0.0E0
        if (I  ==  UK) go to 700
        IP1 = I + 1
!
        DO 680 J = IP1, UK
           if (J  <  N1) go to 670
           L = J - NS
           T = Z(I,L)
           go to 675
  670          T = RM1(I,J+2)
  675          X = X - RM1(J,I) * RV1(J) + T * RV2(J)
           Y = Y - RM1(J,I) * RV2(J) - T * RV1(J)
  680       CONTINUE
!
  700       if (I  <  N1) go to 710
        L = I - NS
        T = Z(I,L)
        go to 715
  710       T = RM1(I,I+2)
  715       call CDIV(X,Y,RM1(I,I),T,RV1(I),RV2(I))
  720    CONTINUE
!     .......... ACCEPTANCE TEST FOR REAL OR COMPLEX
!                EIGENVECTOR AND NORMALIZATION ..........
  740    ITS = ITS + 1
     NORM = 0.0E0
     NORMV = 0.0E0
!
     DO 780 I = 1, UK
        if (ILAMBD  ==  0.0E0) X = ABS(RV1(I))
        if (ILAMBD  /=  0.0E0) X = PYTHAG(RV1(I),RV2(I))
        if (NORMV  >=  X) go to 760
        NORMV = X
        J = I
  760       NORM = NORM + X
  780    CONTINUE
!
     if (NORM  <  GROWTO) go to 840
!     .......... ACCEPT VECTOR ..........
     X = RV1(J)
     if (ILAMBD  ==  0.0E0) X = 1.0E0 / X
     if (ILAMBD  /=  0.0E0) Y = RV2(J)
!
     DO 820 I = 1, UK
        if (ILAMBD  /=  0.0E0) go to 800
        Z(I,S) = RV1(I) * X
        go to 820
  800       call CDIV(RV1(I),RV2(I),X,Y,Z(I,S-1),Z(I,S))
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
     if (ILAMBD  ==  0.0E0) go to 440
     go to 660
!     .......... SET ERROR -- UNACCEPTED EIGENVECTOR ..........
  880    J = 1
     IERR = -K
!     .......... SET REMAINING VECTOR COMPONENTS TO ZERO ..........
  900    DO 920 I = J, N
        Z(I,S) = 0.0E0
        if (ILAMBD  /=  0.0E0) Z(I,S-1) = 0.0E0
  920    CONTINUE
!
  940    S = S + 1
  960    if (IP  ==  (-1)) IP = 0
     if (IP  ==  1) IP = -1
  980 CONTINUE
!
  go to 1001
!     .......... SET ERROR -- UNDERESTIMATE OF EIGENVECTOR
!                SPACE REQUIRED ..........
 1000 if (IERR  /=  0) IERR = IERR - N
  if (IERR  ==  0) IERR = -(2 * N + 1)
 1001 M = S - 1 - ABS(IP)
  return
end
