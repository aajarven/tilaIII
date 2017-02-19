subroutine DMGSBV (M, N, A, IA, NIV, IFLAG, S, P, IP, INHOMO, V, &
     W, WCND)
!
!! DMGSBV orthogonalizes a set of vectors and determines their rank.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DBVSUP
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (MGSBV-S, DMGSBV-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
! Orthogonalize a set of N double precision vectors and determine their
! rank.
!
! **********************************************************************
! INPUT
! **********************************************************************
!   M = dimension of vectors.
!   N = no. of vectors.
!   A = array whose first N cols contain the vectors.
!   IA = first dimension of array A (col length).
!   NIV = number of independent vectors needed.
!   INHOMO = 1 corresponds to having a non-zero particular solution.
!   V = particular solution vector (not included in the pivoting).
!   INDPVT = 1 means pivoting will not be used.
!
! **********************************************************************
! OUTPUT
! **********************************************************************
!   NIV = no. of linear independent vectors in input set.
!     A = matrix whose first NIV cols. contain NIV orthogonal vectors
!         which span the vector space determined by the input vectors.
!   IFLAG
!          = 0 success
!          = 1 incorrect input
!          = 2 rank of new vectors less than N
!   P = decomposition matrix.  P is upper triangular and
!             (old vectors) = (new vectors) * P.
!         The old vectors will be reordered due to pivoting.
!         The dimension of P must be  >=  N*(N+1)/2.
!             (  N*(2*N+1) when N  /=  NFCC )
!   IP = pivoting vector. The dimension of IP must be  >=  N.
!             (  2*N when N  /=  NFCC )
!   S = square of norms of incoming vectors.
!   V = vector which is orthogonal to the vectors of A.
!   W = orthogonalization information for the vector V.
!   WCND = worst case (smallest) norm decrement value of the
!          vectors being orthogonalized  (represents a test
!          for linear dependence of the vectors).
! **********************************************************************
!
!***SEE ALSO  DBVSUP
!***ROUTINES CALLED  DDOT, DPRVEC
!***COMMON BLOCKS    DML18J, DML5MC
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   890921  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DMGSBV
!
  DOUBLE PRECISION DDOT, DPRVEC
  INTEGER I, IA, ICOCO, IFLAG, INDPVT, INHOMO, INTEG, IP(*), IP1, &
       IX, IZ, J, JK, JP, JQ, JY, JZ, K, KD, KJ, KP, L, LIX, LPAR, &
       LR, M, M2, MXNON, N, NDISK, NEQ, NEQIVP, NFCC, NIC, NIV, &
       NIVN, NMNR, NN, NOPG, NP1, NPS, NR, NRM1, NTAPE, NTP, &
       NUMORT, NXPTS
  DOUBLE PRECISION A(IA,*), AE, DOT, EPS, FOURU, P(*), PJP, PSAVE, &
       RE, RY, S(*), SQOVFL, SRU, SV, T, TOL, TWOU, URO, V(*), VL, &
       VNORM, W(*), WCND, Y
!
!
  COMMON /DML18J/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ, &
                  INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,NFCC, &
                  ICOCO
!
  COMMON /DML5MC/ URO,SRU,EPS,SQOVFL,TWOU,FOURU,LPAR
!
!***FIRST EXECUTABLE STATEMENT  DMGSBV
  if (M  >  0 .AND. N  >  0 .AND. IA  >=  M) go to 10
     IFLAG = 1
  go to 280
   10 CONTINUE
!        BEGIN BLOCK PERMITTING ...EXITS TO 270
!           BEGIN BLOCK PERMITTING ...EXITS TO 260
!
           JP = 0
           IFLAG = 0
           NP1 = N + 1
           Y = 0.0D0
           M2 = M/2
!
!              CALCULATE SQUARE OF NORMS OF INCOMING VECTORS AND SEARCH
!              FOR VECTOR WITH LARGEST MAGNITUDE
!
           J = 0
           DO 40 I = 1, N
              VL = DDOT(M,A(1,I),1,A(1,I),1)
              S(I) = VL
              if (N  ==  NFCC) go to 20
                 J = 2*I - 1
                 P(J) = VL
                 IP(J) = J
   20             CONTINUE
              J = J + 1
              P(J) = VL
              IP(J) = J
              if (VL  <=  Y) go to 30
                 Y = VL
                 IX = I
   30             CONTINUE
   40          CONTINUE
           if (INDPVT  /=  1) go to 50
              IX = 1
              Y = P(1)
   50          CONTINUE
           LIX = IX
           if (N  /=  NFCC) LIX = 2*IX - 1
           P(LIX) = P(1)
           S(NP1) = 0.0D0
           if (INHOMO  ==  1) S(NP1) = DDOT(M,V,1,V,1)
           WCND = 1.0D0
           NIVN = NIV
           NIV = 0
!
!           ...EXIT
           if (Y  ==  0.0D0) go to 260
!              *********************************************************
           DO 240 NR = 1, N
!                 BEGIN BLOCK PERMITTING ...EXITS TO 230
!              ......EXIT
                 if (NIVN  ==  NIV) go to 250
                 NIV = NR
                 if (IX  ==  NR) go to 130
!
!                       PIVOTING OF COLUMNS OF P MATRIX
!
                    NN = N
                    LIX = IX
                    LR = NR
                    if (N  ==  NFCC) go to 60
                       NN = NFCC
                       LIX = 2*IX - 1
                       LR = 2*NR - 1
   60                   CONTINUE
                    if (NR  ==  1) go to 80
                       KD = LIX - LR
                       KJ = LR
                       NRM1 = LR - 1
                       DO 70 J = 1, NRM1
                          PSAVE = P(KJ)
                          JK = KJ + KD
                          P(KJ) = P(JK)
                          P(JK) = PSAVE
                          KJ = KJ + NN - J
   70                      CONTINUE
                       JY = JK + NMNR
                       JZ = JY - KD
                       P(JY) = P(JZ)
   80                   CONTINUE
                    IZ = IP(LIX)
                    IP(LIX) = IP(LR)
                    IP(LR) = IZ
                    SV = S(IX)
                    S(IX) = S(NR)
                    S(NR) = SV
                    if (N  ==  NFCC) go to 110
                       if (NR  ==  1) go to 100
                          KJ = LR + 1
                          DO 90 K = 1, NRM1
                             PSAVE = P(KJ)
                             JK = KJ + KD
                             P(KJ) = P(JK)
                             P(JK) = PSAVE
                             KJ = KJ + NFCC - K
   90                         CONTINUE
  100                      CONTINUE
                       IZ = IP(LIX+1)
                       IP(LIX+1) = IP(LR+1)
                       IP(LR+1) = IZ
  110                   CONTINUE
!
!                       PIVOTING OF COLUMNS OF VECTORS
!
                    DO 120 L = 1, M
                       T = A(L,IX)
                       A(L,IX) = A(L,NR)
                       A(L,NR) = T
  120                   CONTINUE
  130                CONTINUE
!
!                    CALCULATE P(NR,NR) AS NORM SQUARED OF PIVOTAL
!                    VECTOR
!
                 JP = JP + 1
                 P(JP) = Y
                 RY = 1.0D0/Y
                 NMNR = N - NR
                 if (N  ==  NFCC) go to 140
                    NMNR = NFCC - (2*NR - 1)
                    JP = JP + 1
                    P(JP) = 0.0D0
                    KP = JP + NMNR
                    P(KP) = Y
  140                CONTINUE
                 if (NR  ==  N .OR. NIVN  ==  NIV) go to 200
!
!                       CALCULATE ORTHOGONAL PROJECTION VECTORS AND
!                       SEARCH FOR LARGEST NORM
!
                    Y = 0.0D0
                    IP1 = NR + 1
                    IX = IP1
!                       ************************************************
                    DO 190 J = IP1, N
                       DOT = DDOT(M,A(1,NR),1,A(1,J),1)
                       JP = JP + 1
                       JQ = JP + NMNR
                       if (N  /=  NFCC) JQ = JQ + NMNR - 1
                       P(JQ) = P(JP) - DOT*(DOT*RY)
                       P(JP) = DOT*RY
                       DO 150 I = 1, M
                          A(I,J) = A(I,J) - P(JP)*A(I,NR)
  150                      CONTINUE
                       if (N  ==  NFCC) go to 170
                          KP = JP + NMNR
                          JP = JP + 1
                          PJP = RY*DPRVEC(M,A(1,NR),A(1,J))
                          P(JP) = PJP
                          P(KP) = -PJP
                          KP = KP + 1
                          P(KP) = RY*DOT
                          DO 160 K = 1, M2
                             L = M2 + K
                             A(K,J) = A(K,J) - PJP*A(L,NR)
                             A(L,J) = A(L,J) + PJP*A(K,NR)
  160                         CONTINUE
                          P(JQ) = P(JQ) - PJP*(PJP/RY)
  170                      CONTINUE
!
!                          TEST FOR CANCELLATION IN RECURRENCE RELATION
!
                       if (P(JQ)  <=  S(J)*SRU) &
                          P(JQ) = DDOT(M,A(1,J),1,A(1,J),1)
                       if (P(JQ)  <=  Y) go to 180
                          Y = P(JQ)
                          IX = J
  180                      CONTINUE
  190                   CONTINUE
                    if (N  /=  NFCC) JP = KP
!                       ************************************************
                    if (INDPVT  ==  1) IX = IP1
!
!                       RECOMPUTE NORM SQUARED OF PIVOTAL VECTOR WITH
!                       SCALAR PRODUCT
!
                    Y = DDOT(M,A(1,IX),1,A(1,IX),1)
!           ............EXIT
                    if (Y  <=  EPS*S(IX)) go to 260
                    WCND = MIN(WCND,Y/S(IX))
  200                CONTINUE
!
!                    COMPUTE ORTHOGONAL PROJECTION OF PARTICULAR
!                    SOLUTION
!
!                 ...EXIT
                 if (INHOMO  /=  1) go to 230
                 LR = NR
                 if (N  /=  NFCC) LR = 2*NR - 1
                 W(LR) = DDOT(M,A(1,NR),1,V,1)*RY
                 DO 210 I = 1, M
                    V(I) = V(I) - W(LR)*A(I,NR)
  210                CONTINUE
!                 ...EXIT
                 if (N  ==  NFCC) go to 230
                 LR = 2*NR
                 W(LR) = RY*DPRVEC(M,V,A(1,NR))
                 DO 220 K = 1, M2
                    L = M2 + K
                    V(K) = V(K) + W(LR)*A(L,NR)
                    V(L) = V(L) - W(LR)*A(K,NR)
  220                CONTINUE
  230             CONTINUE
  240          CONTINUE
  250          CONTINUE
!              *********************************************************
!
!                  TEST FOR LINEAR DEPENDENCE OF PARTICULAR SOLUTION
!
!        ......EXIT
           if (INHOMO  /=  1) go to 270
           if ((N  >  1) .AND. (S(NP1)  <  1.0)) go to 270
           VNORM = DDOT(M,V,1,V,1)
           if (S(NP1)  /=  0.0D0) WCND = MIN(WCND,VNORM/S(NP1))
!        ......EXIT
           if (VNORM  >=  EPS*S(NP1)) go to 270
  260       CONTINUE
        IFLAG = 2
        WCND = EPS
  270    CONTINUE
  280 CONTINUE
  return
end
