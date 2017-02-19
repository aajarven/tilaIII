subroutine HQR (NM, N, LOW, IGH, H, WR, WI, IERR)
!
!! HQR computes the eigenvalues of a real upper Hessenberg matrix ...
!            using the QR method.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C2B
!***TYPE      SINGLE PRECISION (HQR-S, COMQR-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure HQR,
!     NUM. MATH. 14, 219-231(1970) by Martin, Peters, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 359-371(1971).
!
!     This subroutine finds the eigenvalues of a REAL
!     UPPER Hessenberg matrix by the QR method.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, H, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix H.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  BALANC.  If  BALANC  has not been
!          used, set LOW=1 and IGH equal to the order of the matrix, N.
!
!        H contains the upper Hessenberg matrix.  Information about
!          the transformations used in the reduction to Hessenberg
!          form by  ELMHES  or  ORTHES, if performed, is stored
!          in the remaining triangle under the Hessenberg matrix.
!          H is a two-dimensional REAL array, dimensioned H(NM,N).
!
!     On OUTPUT
!
!        H has been destroyed.  Therefore, it must be saved before
!          calling  HQR  if subsequent calculation and back
!          transformation of eigenvectors is to be performed.
!
!        WR and WI contain the real and imaginary parts, respectively,
!          of the eigenvalues.  The eigenvalues are unordered except
!          that complex conjugate pairs of values appear consecutively
!          with the eigenvalue having the positive imaginary part first.
!          If an error exit is made, the eigenvalues should be correct
!          for indices IERR+1, IERR+2, ..., N.  WR and WI are one-
!          dimensional REAL arrays, dimensioned WR(N) and WI(N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after a total of 30*N iterations.
!                     The eigenvalues should be correct for indices
!                     IERR+1, IERR+2, ..., N.
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  HQR
!
  INTEGER I,J,K,L,M,N,EN,LL,MM,NA,NM,IGH,ITN,ITS,LOW,MP2,ENM2,IERR
  REAL H(NM,*),WR(*),WI(*)
  REAL P,Q,R,S,T,W,X,Y,ZZ,NORM,S1,S2
  LOGICAL NOTLAS
!
!***FIRST EXECUTABLE STATEMENT  HQR
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
   60 if (EN  <  LOW) go to 1001
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
     DO 210 J = K, EN
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
     DO 230 I = L, J
        P = X * H(I,K) + Y * H(I,K+1)
        if (.NOT. NOTLAS) go to 220
        P = P + ZZ * H(I,K+2)
        H(I,K+2) = H(I,K+2) - P * R
  220       H(I,K+1) = H(I,K+1) - P * Q
        H(I,K) = H(I,K) - P
  230    CONTINUE
!
  260 CONTINUE
!
  go to 70
!     .......... ONE ROOT FOUND ..........
  270 WR(EN) = X + T
  WI(EN) = 0.0E0
  EN = NA
  go to 60
!     .......... TWO ROOTS FOUND ..........
  280 P = (Y - X) / 2.0E0
  Q = P * P + W
  ZZ = SQRT(ABS(Q))
  X = X + T
  if (Q  <  0.0E0) go to 320
!     .......... REAL PAIR ..........
  ZZ = P + SIGN(ZZ,P)
  WR(NA) = X + ZZ
  WR(EN) = WR(NA)
  if (ZZ  /=  0.0E0) WR(EN) = X - W / ZZ
  WI(NA) = 0.0E0
  WI(EN) = 0.0E0
  go to 330
!     .......... COMPLEX PAIR ..........
  320 WR(NA) = X + P
  WR(EN) = X + P
  WI(NA) = ZZ
  WI(EN) = -ZZ
  330 EN = ENM2
  go to 60
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
end
