subroutine DU11LS (A, MDA, M, N, UB, DB, MODE, NP, KRANK, KSURE, &
     H, W, EB, IC, IR)
!
!! DU11LS performs QR factorization for DLLSIA.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DLLSIA
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (U11LS-S, DU11LS-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   **** Double Precision version of U11LS ****
!
!       This routine performs a QR factorization of A
!       using Householder transformations. Row and
!       column pivots are chosen to reduce the growth
!       of round-off and to help detect possible rank
!       deficiency.
!
!***SEE ALSO  DLLSIA
!***ROUTINES CALLED  DAXPY, DDOT, DNRM2, DSCAL, DSWAP, IDAMAX, ISWAP,
!                    XERMSG
!***REVISION HISTORY  (YYMMDD)
!   810801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced variable.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DU11LS
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DOUBLE PRECISION DDOT,DNRM2
  DIMENSION A(MDA,*),UB(*),DB(*),H(*),W(*),EB(*)
  INTEGER IC(*),IR(*)
!
!        INITIALIZATION
!
!***FIRST EXECUTABLE STATEMENT  DU11LS
  J=0
  KRANK=N
  DO 10 I=1,N
  IC(I)=I
   10 CONTINUE
  DO 12 I=1,M
  IR(I)=I
   12 CONTINUE
!
!        DETERMINE REL AND ABS ERROR VECTORS
!
!
!
!        CALCULATE COL LENGTH
!
  DO 30 I=1,N
  H(I)=DNRM2(M,A(1,I),1)
  W(I)=H(I)
   30 CONTINUE
!
!         INITIALIZE ERROR BOUNDS
!
  DO  40 I=1,N
  EB(I)=MAX(DB(I),UB(I)*H(I))
  UB(I)=EB(I)
  DB(I)=0.0D0
   40 CONTINUE
!
!          DISCARD SELF DEPENDENT COLUMNS
!
  I=1
   50 if ( EB(I) >= H(I)) go to 60
  if ( I == KRANK) go to 70
  I=I+1
  go to 50
!
!          MATRIX REDUCTION
!
   60 CONTINUE
  KK=KRANK
  KRANK=KRANK-1
  if ( MODE == 0) RETURN
  if ( I > NP) go to  64
  call XERMSG ('SLATEC', 'DU11LS', &
     'FIRST NP COLUMNS ARE LINEARLY DEPENDENT', 8, 0)
  KRANK=I-1
  return
   64 CONTINUE
  if ( I > KRANK) go to 70
  call DSWAP(1,EB(I),1,EB(KK),1)
  call DSWAP(1,UB(I),1,UB(KK),1)
  call DSWAP(1,W(I),1,W(KK),1)
  call DSWAP(1,H(I),1,H(KK),1)
  call ISWAP(1,IC(I),1,IC(KK),1)
  call DSWAP(M,A(1,I),1,A(1,KK),1)
  go to 50
!
!           TEST FOR ZERO RANK
!
   70 if ( KRANK > 0) go to 80
  KRANK=0
  KSURE=0
  return
   80 CONTINUE
!
!        M A I N    L O O P
!
  110 CONTINUE
  J=J+1
  JP1=J+1
  JM1=J-1
  KZ=KRANK
  if ( J <= NP) KZ=J
!
!        EACH COL HAS MM=M-J+1 COMPONENTS
!
  MM=M-J+1
!
!         UB DETERMINES COLUMN PIVOT
!
  115 IMIN=J
  if ( H(J) == 0.D0) go to 170
  RMIN=UB(J)/H(J)
  DO 120 I=J,KZ
  if ( UB(I) >= H(I)*RMIN) go to 120
  RMIN=UB(I)/H(I)
  IMIN=I
  120 CONTINUE
!
!     TEST FOR RANK DEFICIENCY
!
  if ( RMIN < 1.0D0) go to 200
  TT=(EB(IMIN)+ABS(DB(IMIN)))/H(IMIN)
  if ( TT >= 1.0D0) go to 170
!     COMPUTE EXACT UB
  DO 125 I=1,JM1
  W(I)=A(I,IMIN)
  125 CONTINUE
  L=JM1
  130 W(L)=W(L)/A(L,L)
  if ( L == 1) go to 150
  LM1=L-1
  DO 140 I=L,JM1
  W(LM1)=W(LM1)-A(LM1,I)*W(I)
  140 CONTINUE
  L=LM1
  go to 130
  150 TT=EB(IMIN)
  DO 160 I=1,JM1
  TT=TT+ABS(W(I))*EB(I)
  160 CONTINUE
  UB(IMIN)=TT
  if ( UB(IMIN)/H(IMIN) >= 1.0D0) go to 170
  go to 200
!
!        MATRIX REDUCTION
!
  170 CONTINUE
  KK=KRANK
  KRANK=KRANK-1
  KZ=KRANK
  if ( MODE == 0) RETURN
  if ( J > NP) go to 172
  call XERMSG ('SLATEC', 'DU11LS', &
     'FIRST NP COLUMNS ARE LINEARLY DEPENDENT', 8, 0)
  KRANK=J-1
  return
  172 CONTINUE
  if ( IMIN > KRANK) go to 180
  call ISWAP(1,IC(IMIN),1,IC(KK),1)
  call DSWAP(M,A(1,IMIN),1,A(1,KK),1)
  call DSWAP(1,EB(IMIN),1,EB(KK),1)
  call DSWAP(1,UB(IMIN),1,UB(KK),1)
  call DSWAP(1,DB(IMIN),1,DB(KK),1)
  call DSWAP(1,W(IMIN),1,W(KK),1)
  call DSWAP(1,H(IMIN),1,H(KK),1)
  180 if ( J > KRANK) go to 300
  go to 115
!
!        COLUMN PIVOT
!
  200 if ( IMIN == J) go to 230
  call DSWAP(1,H(J),1,H(IMIN),1)
  call DSWAP(M,A(1,J),1,A(1,IMIN),1)
  call DSWAP(1,EB(J),1,EB(IMIN),1)
  call DSWAP(1,UB(J),1,UB(IMIN),1)
  call DSWAP(1,DB(J),1,DB(IMIN),1)
  call DSWAP(1,W(J),1,W(IMIN),1)
  call ISWAP(1,IC(J),1,IC(IMIN),1)
!
!        ROW PIVOT
!
  230 CONTINUE
  JMAX=IDAMAX(MM,A(J,J),1)
  JMAX=JMAX+J-1
  if ( JMAX == J) go to 240
  call DSWAP(N,A(J,1),MDA,A(JMAX,1),MDA)
  call ISWAP(1,IR(J),1,IR(JMAX),1)
  240 CONTINUE
!
!     APPLY HOUSEHOLDER TRANSFORMATION
!
  TN=DNRM2(MM,A(J,J),1)
  if ( TN == 0.0D0) go to 170
  if ( A(J,J) /= 0.0D0) TN=SIGN(TN,A(J,J))
  call DSCAL(MM,1.0D0/TN,A(J,J),1)
  A(J,J)=A(J,J)+1.0D0
  if ( J == N) go to 250
  DO 248 I=JP1,N
  BB=-DDOT(MM,A(J,J),1,A(J,I),1)/A(J,J)
  call DAXPY(MM,BB,A(J,J),1,A(J,I),1)
  if ( I <= NP) go to 248
  if ( H(I) == 0.0D0) go to 248
  TT=1.0D0-(ABS(A(J,I))/H(I))**2
  TT=MAX(TT,0.0D0)
  T=TT
  TT=1.0D0+.05D0*TT*(H(I)/W(I))**2
  if ( TT == 1.0D0) go to 244
  H(I)=H(I)*SQRT(T)
  go to 246
  244 CONTINUE
  H(I)=DNRM2(M-J,A(J+1,I),1)
  W(I)=H(I)
  246 CONTINUE
  248 CONTINUE
  250 CONTINUE
  H(J)=A(J,J)
  A(J,J)=-TN
!
!
!          UPDATE UB, DB
!
  UB(J)=UB(J)/ABS(A(J,J))
  DB(J)=(SIGN(EB(J),DB(J))+DB(J))/A(J,J)
  if ( J == KRANK) go to 300
  DO 260 I=JP1,KRANK
  UB(I)=UB(I)+ABS(A(J,I))*UB(J)
  DB(I)=DB(I)-A(J,I)*DB(J)
  260 CONTINUE
  go to 110
!
!        E N D    M A I N    L O O P
!
  300 CONTINUE
!
!        COMPUTE KSURE
!
  KM1=KRANK-1
  DO 318 I=1,KM1
  IS=0
  KMI=KRANK-I
  DO 315 II=1,KMI
  if ( UB(II) <= UB(II+1)) go to 315
  IS=1
  TEMP=UB(II)
  UB(II)=UB(II+1)
  UB(II+1)=TEMP
  315 CONTINUE
  if ( IS == 0) go to 320
  318 CONTINUE
  320 CONTINUE
  KSURE=0
  SUM=0.0D0
  DO 328 I=1,KRANK
  R2=UB(I)*UB(I)
  if ( R2+SUM >= 1.0D0) go to 330
  SUM=SUM+R2
  KSURE=KSURE+1
  328 CONTINUE
  330 CONTINUE
!
!     if SYSTEM IS OF REDUCED RANK AND MODE = 2
!     COMPLETE THE DECOMPOSITION FOR SHORTEST LEAST SQUARES SOLUTION
!
  if ( KRANK == N .OR. MODE < 2) go to 360
  NMK=N-KRANK
  KP1=KRANK+1
  I=KRANK
  340 TN=DNRM2(NMK,A(I,KP1),MDA)/A(I,I)
  TN=A(I,I)*SQRT(1.0D0+TN*TN)
  call DSCAL(NMK,1.0D0/TN,A(I,KP1),MDA)
  W(I)=A(I,I)/TN+1.0D0
  A(I,I)=-TN
  if ( I == 1) go to 350
  IM1=I-1
  DO 345 II=1,IM1
  TT=-DDOT(NMK,A(II,KP1),MDA,A(I,KP1),MDA)/W(I)
  TT=TT-A(II,I)
  call DAXPY(NMK,TT,A(I,KP1),MDA,A(II,KP1),MDA)
  A(II,I)=A(II,I)+TT*W(I)
  345 CONTINUE
  I=I-1
  go to 340
  350 CONTINUE
  360 CONTINUE
  return
end
