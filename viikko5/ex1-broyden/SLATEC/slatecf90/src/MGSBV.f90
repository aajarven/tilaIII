subroutine MGSBV (M, N, A, IA, NIV, IFLAG, S, P, IP, INHOMO, V, W, &
     WCND)
!
!! MGSBV is subsidiary to BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (MGSBV-S, DMGSBV-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
! Orthogonalize a set of N real vectors and determine their rank
!
! **********************************************************************
! INPUT
! **********************************************************************
!   M = Dimension of vectors
!   N = No. of vectors
!   A = Array whose first N cols contain the vectors
!   IA = First dimension of array A (col length)
!   NIV = Number of independent vectors needed
!   INHOMO = 1 Corresponds to having a non-zero particular solution
!   V = Particular solution vector (not included in the pivoting)
!   INDPVT = 1 Means pivoting will not be used
!
! **********************************************************************
! OUTPUT
! **********************************************************************
!   NIV = No. of linear independent vectors in input set
!     A = Matrix whose first NIV cols. contain NIV orthogonal vectors
!         which span the vector space determined by the input vectors
!   IFLAG
!          = 0 success
!          = 1 incorrect input
!          = 2 rank of new vectors less than N
!   P = Decomposition matrix.  P is upper triangular and
!             (old vectors) = (new vectors) * P.
!         The old vectors will be reordered due to pivoting
!         The dimension of p must be  >=  N*(N+1)/2.
!             (  N*(2*N+1) when N  /=  NFCC )
!   IP = Pivoting vector. The dimension of IP must be  >=  N.
!             (  2*N when N  /=  NFCC )
!   S = Square of norms of incoming vectors
!   V = Vector which is orthogonal to the vectors of A
!   W = Orthogonalization information for the vector V
!   WCND = Worst case (smallest) norm decrement value of the
!          vectors being orthogonalized  (represents a test
!          for linear dependence of the vectors)
! **********************************************************************
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  PRVEC, SDOT
!***COMMON BLOCKS    ML18JR, ML5MCO
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  MGSBV
!
  DIMENSION A(IA,*),V(*),W(*),P(*),IP(*),S(*)
!
!
  COMMON /ML18JR/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ, &
                  INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,NFCC, &
                  ICOCO
!
  COMMON /ML5MCO/ URO,SRU,EPS,SQOVFL,TWOU,FOURU,LPAR
!
!***FIRST EXECUTABLE STATEMENT  MGSBV
  if ( M  >  0  .AND.  N  >  0  .AND.  IA  >=  M) go to 10
  IFLAG=1
  return
!
   10 JP=0
  IFLAG=0
  NP1=N+1
  Y=0.0
  M2=M/2
!
!     CALCULATE SQUARE OF NORMS OF INCOMING VECTORS AND SEARCH FOR
!     VECTOR WITH LARGEST MAGNITUDE
!
  J=0
  DO 30 I=1,N
  VL=SDOT(M,A(1,I),1,A(1,I),1)
  S(I)=VL
  if (N  ==  NFCC) go to 25
  J=2*I-1
  P(J)=VL
  IP(J)=J
   25 J=J+1
  P(J)=VL
  IP(J)=J
  if ( VL  <=  Y) go to 30
  Y=VL
  IX=I
   30 CONTINUE
  if (INDPVT  /=  1) go to 33
  IX=1
  Y=P(1)
   33 LIX=IX
  if (N  /=  NFCC) LIX=2*IX-1
  P(LIX)=P(1)
  S(NP1)=0.
  if (INHOMO  ==  1) S(NP1)=SDOT(M,V,1,V,1)
  WCND=1.
  NIVN=NIV
  NIV=0
!
  if ( Y  ==  0.0) go to 170
! **********************************************************************
  DO 140 NR=1,N
  if (NIVN  ==  NIV) go to 150
  NIV=NR
  if ( IX  ==  NR) go to 80
!
!     PIVOTING OF COLUMNS OF P MATRIX
!
  NN=N
  LIX=IX
  LR=NR
  if (N  ==  NFCC) go to 40
  NN=NFCC
  LIX=2*IX-1
  LR=2*NR-1
   40 if ( NR  ==  1) go to 60
  KD=LIX-LR
  KJ=LR
  NRM1=LR-1
  DO 50 J=1,NRM1
  PSAVE=P(KJ)
  JK=KJ+KD
  P(KJ)=P(JK)
  P(JK)=PSAVE
   50 KJ=KJ+NN-J
  JY=JK+NMNR
  JZ=JY-KD
  P(JY)=P(JZ)
   60 IZ=IP(LIX)
  IP(LIX)=IP(LR)
  IP(LR)=IZ
  SV=S(IX)
  S(IX)=S(NR)
  S(NR)=SV
  if (N  ==  NFCC) go to 69
  if (NR  ==  1) go to 67
  KJ=LR+1
  DO 65 K=1,NRM1
  PSAVE=P(KJ)
  JK=KJ+KD
  P(KJ)=P(JK)
  P(JK)=PSAVE
   65 KJ=KJ+NFCC-K
   67 IZ=IP(LIX+1)
  IP(LIX+1)=IP(LR+1)
  IP(LR+1)=IZ
!
!     PIVOTING OF COLUMNS OF VECTORS
!
   69 DO 70 L=1,M
  T=A(L,IX)
  A(L,IX)=A(L,NR)
   70 A(L,NR)=T
!
!     CALCULATE P(NR,NR) AS NORM SQUARED OF PIVOTAL VECTOR
!
   80 JP=JP+1
  P(JP)=Y
  RY=1.0/Y
  NMNR=N-NR
  if (N  ==  NFCC) go to 85
  NMNR=NFCC-(2*NR-1)
  JP=JP+1
  P(JP)=0.
  KP=JP+NMNR
  P(KP)=Y
   85 if ( NR  ==  N  .OR.  NIVN  ==  NIV) go to 125
!
!    CALCULATE ORTHOGONAL PROJECTION VECTORS AND SEARCH FOR LARGEST NORM
!
  Y=0.0
  IP1=NR+1
  IX=IP1
!     ****************************************
  DO 120 J=IP1,N
  DOT=SDOT(M,A(1,NR),1,A(1,J),1)
  JP=JP+1
  JQ=JP+NMNR
  if (N  /=  NFCC) JQ=JQ+NMNR-1
  P(JQ)=P(JP)-DOT*(DOT*RY)
  P(JP)=DOT*RY
  DO 90 I = 1,M
   90 A(I,J)=A(I,J)-P(JP)*A(I,NR)
  if (N  ==  NFCC) go to 99
  KP=JP+NMNR
  JP=JP+1
  PJP=RY*PRVEC(M,A(1,NR),A(1,J))
  P(JP)=PJP
  P(KP)=-PJP
  KP=KP+1
  P(KP)=RY*DOT
  DO 95 K=1,M2
  L=M2+K
  A(K,J)=A(K,J)-PJP*A(L,NR)
   95 A(L,J)=A(L,J)+PJP*A(K,NR)
  P(JQ)=P(JQ)-PJP*(PJP/RY)
!
!     TEST FOR CANCELLATION IN RECURRENCE RELATION
!
   99 if ( P(JQ)  >  S(J)*SRU) go to 100
  P(JQ)=SDOT(M,A(1,J),1,A(1,J),1)
  100 if ( P(JQ)  <=  Y) go to 120
  Y=P(JQ)
  IX=J
  120 CONTINUE
  if (N  /=  NFCC) JP=KP
!     ****************************************
  if ( INDPVT  ==  1) IX=IP1
!
!     RECOMPUTE NORM SQUARED OF PIVOTAL VECTOR WITH SCALAR PRODUCT
!
  Y=SDOT(M,A(1,IX),1,A(1,IX),1)
  if ( Y   <=   EPS*S(IX))  go to 170
  WCND=MIN(WCND,Y/S(IX))
!
!     COMPUTE ORTHOGONAL PROJECTION OF PARTICULAR SOLUTION
!
  125 if ( INHOMO  /=  1) go to 140
  LR=NR
  if (N  /=  NFCC) LR=2*NR-1
  W(LR)=SDOT(M,A(1,NR),1,V,1)*RY
  DO 130 I=1,M
  130 V(I)=V(I)-W(LR)*A(I,NR)
  if (N  ==  NFCC) go to 140
  LR=2*NR
  W(LR)=RY*PRVEC(M,V,A(1,NR))
  DO 135 K=1,M2
  L=M2+K
  V(K)=V(K)+W(LR)*A(L,NR)
  135 V(L)=V(L)-W(LR)*A(K,NR)
  140 CONTINUE
! **********************************************************************
!
!     TEST FOR LINEAR DEPENDENCE OF PARTICULAR SOLUTION
!
  150 if ( INHOMO  /=  1) RETURN
  if ((N  >  1) .AND. (S(NP1)  <  1.0)) RETURN
  VNORM=SDOT(M,V,1,V,1)
  if (S(NP1)  /=  0.) WCND=MIN(WCND,VNORM/S(NP1))
  if ( VNORM  >=  EPS*S(NP1)) RETURN
  170 IFLAG=2
  WCND=EPS
  return
end
