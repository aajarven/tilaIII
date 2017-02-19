subroutine HQR2 (NM, N, LOW, IGH, H, WR, WI, Z, IERR)
!
!! HQR2 computes the eigenvalues and eigenvectors of a real upper ...
!            Hessenberg matrix using QR method.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C2B
!***TYPE      SINGLE PRECISION (HQR2-S, COMQR2-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure HQR2,
!     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!
!     This subroutine finds the eigenvalues and eigenvectors
!     of a REAL UPPER Hessenberg matrix by the QR method.  The
!     eigenvectors of a REAL GENERAL matrix can also be found
!     if  ELMHES  and  ELTRAN  or  ORTHES  and  ORTRAN  have
!     been used to reduce this general matrix to Hessenberg form
!     and to accumulate the similarity transformations.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, H and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix H.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  BALANC.  If  BALANC  has not been
!          used, set LOW=1 and IGH equal to the order of the matrix, N.
!
!        H contains the upper Hessenberg matrix.  H is a two-dimensional
!          REAL array, dimensioned H(NM,N).
!
!        Z contains the transformation matrix produced by  ELTRAN
!          after the reduction by  ELMHES, or by  ORTRAN  after the
!          reduction by  ORTHES, if performed.  If the eigenvectors
!          of the Hessenberg matrix are desired, Z must contain the
!          identity matrix.  Z is a two-dimensional REAL array,
!          dimensioned Z(NM,M).
!
!     On OUTPUT
!
!        H has been destroyed.
!
!        WR and WI contain the real and imaginary parts, respectively,
!          of the eigenvalues.  The eigenvalues are unordered except
!          that complex conjugate pairs of values appear consecutively
!          with the eigenvalue having the positive imaginary part first.
!          If an error exit is made, the eigenvalues should be correct
!          for indices IERR+1, IERR+2, ..., N.  WR and WI are one-
!          dimensional REAL arrays, dimensioned WR(N) and WI(N).
!
!        Z contains the real and imaginary parts of the eigenvectors.
!          If the J-th eigenvalue is real, the J-th column of Z
!          contains its eigenvector.  If the J-th eigenvalue is complex
!          with positive imaginary part, the J-th and (J+1)-th
!          columns of Z contain the real and imaginary parts of its
!          eigenvector.  The eigenvectors are unnormalized.  If an
!          error exit is made, none of the eigenvectors has been found.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after a total of 30*N iterations.
!                     The eigenvalues should be correct for indices
!                     IERR+1, IERR+2, ..., N, but no eigenvectors are
!                     computed.
!
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
!***ROUTINES CALLED  CDIV
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  HQR2
!
  INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN
  INTEGER IGH,ITN,ITS,LOW,MP2,ENM2,IERR
  REAL H(NM,*),WR(*),WI(*),Z(NM,*)
  REAL P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,S1,S2
  LOGICAL NOTLAS
!
!***FIRST EXECUTABLE STATEMENT  HQR2
  IERR = 0
  NORM = 0.0E0
  K = 1
!     .......... STORE ROOTS ISOLATED BY BALANC
!                AND COMPUTE MATRIX NORM ..........
  DO 50 I = 1, N
!
     DO 40 J = K, N
   40    NORM = NORM + ABS(H(I,J))
!
     K = I
     if (I  >=  LOW .AND. I  <=  IGH) go to 50
     WR(I) = H(I,I)
     WI(I) = 0.0E0
   50 CONTINUE
!
  EN = IGH
  T = 0.0E0
  ITN = 30*N
!     .......... SEARCH FOR NEXT EIGENVALUES ..........
   60 if (EN  <  LOW) go to 340
  ITS = 0
  NA = EN - 1
  ENM2 = NA - 1
!     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
!                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
   70 DO 80 LL = LOW, EN
     L = EN + LOW - LL
     if (L  ==  LOW) go to 100
     S = ABS(H(L-1,L-1)) + ABS(H(L,L))
     if (S  ==  0.0E0) S = NORM
     S2 = S + ABS(H(L,L-1))
     if (S2  ==  S) go to 100
   80 CONTINUE
!     .......... FORM SHIFT ..........
  100 X = H(EN,EN)
  if (L  ==  EN) go to 270
  Y = H(NA,NA)
  W = H(EN,NA) * H(NA,EN)
  if (L  ==  NA) go to 280
  if (ITN  ==  0) go to 1000
  if (ITS  /=  10 .AND. ITS  /=  20) go to 130
!     .......... FORM EXCEPTIONAL SHIFT ..........
  T = T + X
!
  DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
!
  S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))
  X = 0.75E0 * S
  Y = X
  W = -0.4375E0 * S * S
  130 ITS = ITS + 1
  ITN = ITN - 1
!     .......... LOOK FOR TWO CONSECUTIVE SMALL
!                SUB-DIAGONAL ELEMENTS.
!                FOR M=EN-2 STEP -1 UNTIL L DO -- ..........
  DO 140 MM = L, ENM2
     M = ENM2 + L - MM
     ZZ = H(M,M)
     R = X - ZZ
     S = Y - ZZ
     P = (R * S - W) / H(M+1,M) + H(M,M+1)
     Q = H(M+1,M+1) - ZZ - R - S
     R = H(M+2,M+1)
     S = ABS(P) + ABS(Q) + ABS(R)
     P = P / S
     Q = Q / S
     R = R / S
     if (M  ==  L) go to 150
     S1 = ABS(P) * (ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))
     S2 = S1 + ABS(H(M,M-1)) * (ABS(Q) + ABS(R))
     if (S2  ==  S1) go to 150
  140 CONTINUE
!
  150 MP2 = M + 2
!
  DO 160 I = MP2, EN
     H(I,I-2) = 0.0E0
     if (I  ==  MP2) go to 160
     H(I,I-3) = 0.0E0
  160 CONTINUE
!     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND
!                COLUMNS M TO EN ..........
  DO 260 K = M, NA
     NOTLAS = K  /=  NA
     if (K  ==  M) go to 170
     P = H(K,K-1)
     Q = H(K+1,K-1)
     R = 0.0E0
     if (NOTLAS) R = H(K+2,K-1)
     X = ABS(P) + ABS(Q) + ABS(R)
     if (X  ==  0.0E0) go to 260
     P = P / X
     Q = Q / X
     R = R / X
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)
     if (K  ==  M) go to 180
     H(K,K-1) = -S * X
     go to 190
  180    if (L  /=  M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
     X = P / S
     Y = Q / S
     ZZ = R / S
     Q = Q / P
     R = R / P
!     .......... ROW MODIFICATION ..........
     DO 210 J = K, N
        P = H(K,J) + Q * H(K+1,J)
        if (.NOT. NOTLAS) go to 200
        P = P + R * H(K+2,J)
        H(K+2,J) = H(K+2,J) - P * ZZ
  200       H(K+1,J) = H(K+1,J) - P * Y
        H(K,J) = H(K,J) - P * X
  210    CONTINUE
!
     J = MIN(EN,K+3)
!     .......... COLUMN MODIFICATION ..........
     DO 230 I = 1, J
        P = X * H(I,K) + Y * H(I,K+1)
        if (.NOT. NOTLAS) go to 220
        P = P + ZZ * H(I,K+2)
        H(I,K+2) = H(I,K+2) - P * R
  220       H(I,K+1) = H(I,K+1) - P * Q
        H(I,K) = H(I,K) - P
  230    CONTINUE
!     .......... ACCUMULATE TRANSFORMATIONS ..........
     DO 250 I = LOW, IGH
        P = X * Z(I,K) + Y * Z(I,K+1)
        if (.NOT. NOTLAS) go to 240
        P = P + ZZ * Z(I,K+2)
        Z(I,K+2) = Z(I,K+2) - P * R
  240       Z(I,K+1) = Z(I,K+1) - P * Q
        Z(I,K) = Z(I,K) - P
  250    CONTINUE
!
  260 CONTINUE
!
  go to 70
!     .......... ONE ROOT FOUND ..........
  270 H(EN,EN) = X + T
  WR(EN) = H(EN,EN)
  WI(EN) = 0.0E0
  EN = NA
  go to 60
!     .......... TWO ROOTS FOUND ..........
  280 P = (Y - X) / 2.0E0
  Q = P * P + W
  ZZ = SQRT(ABS(Q))
  H(EN,EN) = X + T
  X = H(EN,EN)
  H(NA,NA) = Y + T
  if (Q  <  0.0E0) go to 320
!     .......... REAL PAIR ..........
  ZZ = P + SIGN(ZZ,P)
  WR(NA) = X + ZZ
  WR(EN) = WR(NA)
  if (ZZ  /=  0.0E0) WR(EN) = X - W / ZZ
  WI(NA) = 0.0E0
  WI(EN) = 0.0E0
  X = H(EN,NA)
  S = ABS(X) + ABS(ZZ)
  P = X / S
  Q = ZZ / S
  R = SQRT(P*P+Q*Q)
  P = P / R
  Q = Q / R
!     .......... ROW MODIFICATION ..........
  DO 290 J = NA, N
     ZZ = H(NA,J)
     H(NA,J) = Q * ZZ + P * H(EN,J)
     H(EN,J) = Q * H(EN,J) - P * ZZ
  290 CONTINUE
!     .......... COLUMN MODIFICATION ..........
  DO 300 I = 1, EN
     ZZ = H(I,NA)
     H(I,NA) = Q * ZZ + P * H(I,EN)
     H(I,EN) = Q * H(I,EN) - P * ZZ
  300 CONTINUE
!     .......... ACCUMULATE TRANSFORMATIONS ..........
  DO 310 I = LOW, IGH
     ZZ = Z(I,NA)
     Z(I,NA) = Q * ZZ + P * Z(I,EN)
     Z(I,EN) = Q * Z(I,EN) - P * ZZ
  310 CONTINUE
!
  go to 330
!     .......... COMPLEX PAIR ..........
  320 WR(NA) = X + P
  WR(EN) = X + P
  WI(NA) = ZZ
  WI(EN) = -ZZ
  330 EN = ENM2
  go to 60
!     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
!                VECTORS OF UPPER TRIANGULAR FORM ..........
  340 if (NORM  ==  0.0E0) go to 1001
!     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
  DO 800 NN = 1, N
     EN = N + 1 - NN
     P = WR(EN)
     Q = WI(EN)
     NA = EN - 1
     if (Q) 710, 600, 800
!     .......... REAL VECTOR ..........
  600    M = EN
     H(EN,EN) = 1.0E0
     if (NA  ==  0) go to 800
!     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
     DO 700 II = 1, NA
        I = EN - II
        W = H(I,I) - P
        R = H(I,EN)
        if (M  >  NA) go to 620
!
        DO 610 J = M, NA
  610       R = R + H(I,J) * H(J,EN)
!
  620       if (WI(I)  >=  0.0E0) go to 630
        ZZ = W
        S = R
        go to 700
  630       M = I
        if (WI(I)  /=  0.0E0) go to 640
        T = W
        if (T  /=  0.0E0) go to 635
        T = NORM
  632       T = 0.5E0*T
        if (NORM + T  >  NORM) go to 632
        T = 2.0E0*T
  635       H(I,EN) = -R / T
        go to 700
!     .......... SOLVE REAL EQUATIONS ..........
  640       X = H(I,I+1)
        Y = H(I+1,I)
        Q = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I)
        T = (X * S - ZZ * R) / Q
        H(I,EN) = T
        if (ABS(X)  <=  ABS(ZZ)) go to 650
        H(I+1,EN) = (-R - W * T) / X
        go to 700
  650       H(I+1,EN) = (-S - Y * T) / ZZ
  700    CONTINUE
!     .......... END REAL VECTOR ..........
     go to 800
!     .......... COMPLEX VECTOR ..........
  710    M = NA
!     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
!                EIGENVECTOR MATRIX IS TRIANGULAR ..........
     if (ABS(H(EN,NA))  <=  ABS(H(NA,EN))) go to 720
     H(NA,NA) = Q / H(EN,NA)
     H(NA,EN) = -(H(EN,EN) - P) / H(EN,NA)
     go to 730
  720    call CDIV(0.0E0,-H(NA,EN),H(NA,NA)-P,Q,H(NA,NA),H(NA,EN))
  730    H(EN,NA) = 0.0E0
     H(EN,EN) = 1.0E0
     ENM2 = NA - 1
     if (ENM2  ==  0) go to 800
!     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
     DO 790 II = 1, ENM2
        I = NA - II
        W = H(I,I) - P
        RA = 0.0E0
        SA = H(I,EN)
!
        DO 760 J = M, NA
           RA = RA + H(I,J) * H(J,NA)
           SA = SA + H(I,J) * H(J,EN)
  760       CONTINUE
!
        if (WI(I)  >=  0.0E0) go to 770
        ZZ = W
        R = RA
        S = SA
        go to 790
  770       M = I
        if (WI(I)  /=  0.0E0) go to 780
        call CDIV(-RA,-SA,W,Q,H(I,NA),H(I,EN))
        go to 790
!     .......... SOLVE COMPLEX EQUATIONS ..........
  780       X = H(I,I+1)
        Y = H(I+1,I)
        VR = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I) - Q * Q
        VI = (WR(I) - P) * 2.0E0 * Q
        if (VR  /=  0.0E0 .OR. VI  /=  0.0E0) go to 783
        S1 = NORM * (ABS(W)+ABS(Q)+ABS(X)+ABS(Y)+ABS(ZZ))
        VR = S1
  782       VR = 0.5E0*VR
        if (S1 + VR  >  S1) go to 782
        VR = 2.0E0*VR
  783       call CDIV(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA,VR,VI, &
                  H(I,NA),H(I,EN))
        if (ABS(X)  <=  ABS(ZZ) + ABS(Q)) go to 785
        H(I+1,NA) = (-RA - W * H(I,NA) + Q * H(I,EN)) / X
        H(I+1,EN) = (-SA - W * H(I,EN) - Q * H(I,NA)) / X
        go to 790
  785       call CDIV(-R-Y*H(I,NA),-S-Y*H(I,EN),ZZ,Q, &
                  H(I+1,NA),H(I+1,EN))
  790    CONTINUE
!     .......... END COMPLEX VECTOR ..........
  800 CONTINUE
!     .......... END BACK SUBSTITUTION.
!                VECTORS OF ISOLATED ROOTS ..........
  DO 840 I = 1, N
     if (I  >=  LOW .AND. I  <=  IGH) go to 840
!
     DO 820 J = I, N
  820    Z(I,J) = H(I,J)
!
  840 CONTINUE
!     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
!                VECTORS OF ORIGINAL FULL MATRIX.
!                FOR J=N STEP -1 UNTIL LOW DO -- ..........
  DO 880 JJ = LOW, N
     J = N + LOW - JJ
     M = MIN(J,IGH)
!
     DO 880 I = LOW, IGH
        ZZ = 0.0E0
!
        DO 860 K = LOW, M
  860       ZZ = ZZ + Z(I,K) * H(K,J)
!
        Z(I,J) = ZZ
  880 CONTINUE
!
  go to 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
end
