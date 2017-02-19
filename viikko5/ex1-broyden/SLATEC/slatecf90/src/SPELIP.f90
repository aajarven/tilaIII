subroutine SPELIP (INTL, IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, &
     BETA, C, D, N, NBDCND, BDC, GAMA, BDD, XNU, COFX, COFY, AN, BN, &
     CN, DN, UN, ZN, AM, BM, CM, DM, UM, ZM, GRHS, USOL, IDMN, W, &
     PERTRB, IERROR)
!
!! SPELIP is subsidiary to SEPELI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (SPELIP-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     SPELIP sets up vectors and arrays for input to BLKTRI
!     and computes a second order solution in USOL.  A return jump to
!     SEPELI occurs if IORDER=2.  If IORDER=4 a fourth order
!     solution is generated in USOL.
!
!***SEE ALSO  SEPELI
!***ROUTINES CALLED  BLKTRI, CHKSNG, DEFER, MINSOL, ORTHOG, TRISP
!***COMMON BLOCKS    SPLPCM
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  SPELIP
!
  DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     , &
                  W(*)
  DIMENSION       GRHS(IDMN,*)           ,USOL(IDMN,*)
  DIMENSION       AN(*)      ,BN(*)      ,CN(*)      ,DN(*)      , &
                  UN(*)      ,ZN(*)
  DIMENSION       AM(*)      ,BM(*)      ,CM(*)      ,DM(*)      , &
                  UM(*)      ,ZM(*)
  COMMON /SPLPCM/ KSWX       ,KSWY       ,K          ,L          , &
                  AIT        ,BIT        ,CIT        ,DIT        , &
                  MIT        ,NIT        ,IS         ,MS         , &
                  JS         ,NS         ,DLX        ,DLY        , &
                  TDLX3      ,TDLY3      ,DLX4       ,DLY4
  LOGICAL         SINGLR
  EXTERNAL        COFX       ,COFY
!***FIRST EXECUTABLE STATEMENT  SPELIP
  KSWX = MBDCND+1
  KSWY = NBDCND+1
  K = M+1
  L = N+1
  AIT = A
  BIT = B
  CIT = C
  DIT = D
!
!     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR
!     AND NON-SPECIFIED BOUNDARIES.
!
  DO  20 I=2,M
     DO  10 J=2,N
        USOL(I,J) = GRHS(I,J)
   10    CONTINUE
   20 CONTINUE
  if (KSWX == 2 .OR. KSWX == 3) go to  40
  DO  30 J=2,N
     USOL(1,J) = GRHS(1,J)
   30 CONTINUE
   40 CONTINUE
  if (KSWX == 2 .OR. KSWX == 5) go to  60
  DO  50 J=2,N
     USOL(K,J) = GRHS(K,J)
   50 CONTINUE
   60 CONTINUE
  if (KSWY == 2 .OR. KSWY == 3) go to  80
  DO  70 I=2,M
     USOL(I,1) = GRHS(I,1)
   70 CONTINUE
   80 CONTINUE
  if (KSWY == 2 .OR. KSWY == 5) go to 100
  DO  90 I=2,M
     USOL(I,L) = GRHS(I,L)
   90 CONTINUE
  100 CONTINUE
  if (KSWX /= 2 .AND. KSWX /= 3 .AND. KSWY /= 2 .AND. KSWY /= 3) &
      USOL(1,1) = GRHS(1,1)
  if (KSWX /= 2 .AND. KSWX /= 5 .AND. KSWY /= 2 .AND. KSWY /= 3) &
      USOL(K,1) = GRHS(K,1)
  if (KSWX /= 2 .AND. KSWX /= 3 .AND. KSWY /= 2 .AND. KSWY /= 5) &
      USOL(1,L) = GRHS(1,L)
  if (KSWX /= 2 .AND. KSWX /= 5 .AND. KSWY /= 2 .AND. KSWY /= 5) &
      USOL(K,L) = GRHS(K,L)
  I1 = 1
!
!     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
!
  MP = 1
  NP = 1
  if (KSWX  ==  1) MP = 0
  if (KSWY  ==  1) NP = 0
!
!     SET DLX,DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED
!     IN NINT,MINT
!
  DLX = (BIT-AIT)/M
  MIT = K-1
  if (KSWX  ==  2) MIT = K-2
  if (KSWX  ==  4) MIT = K
  DLY = (DIT-CIT)/N
  NIT = L-1
  if (KSWY  ==  2) NIT = L-2
  if (KSWY  ==  4) NIT = L
  TDLX3 = 2.0*DLX**3
  DLX4 = DLX**4
  TDLY3 = 2.0*DLY**3
  DLY4 = DLY**4
!
!     SET SUBSCRIPT LIMITS FOR PORTION OF ARRAY TO INPUT TO BLKTRI
!
  IS = 1
  JS = 1
  if (KSWX == 2 .OR. KSWX == 3) IS = 2
  if (KSWY == 2 .OR. KSWY == 3) JS = 2
  NS = NIT+JS-1
  MS = MIT+IS-1
!
!     SET X - DIRECTION
!
  DO 110 I=1,MIT
     XI = AIT+(IS+I-2)*DLX
     call COFX (XI,AI,BI,CI)
     AXI = (AI/DLX-0.5*BI)/DLX
     BXI = -2.*AI/DLX**2+CI
     CXI = (AI/DLX+0.5*BI)/DLX
     AM(I) = AXI
     BM(I) = BXI
     CM(I) = CXI
  110 CONTINUE
!
!     SET Y DIRECTION
!
  DO 120 J=1,NIT
     YJ = CIT+(JS+J-2)*DLY
     call COFY (YJ,DJ,EJ,FJ)
     DYJ = (DJ/DLY-0.5*EJ)/DLY
     EYJ = (-2.*DJ/DLY**2+FJ)
     FYJ = (DJ/DLY+0.5*EJ)/DLY
     AN(J) = DYJ
     BN(J) = EYJ
     CN(J) = FYJ
  120 CONTINUE
!
!     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC
!
  AX1 = AM(1)
  CXM = CM(MIT)
  go to (170,130,150,160,140),KSWX
!
!     DIRICHLET-DIRICHLET IN X DIRECTION
!
  130 AM(1) = 0.0
  CM(MIT) = 0.0
  go to 170
!
!     MIXED-DIRICHLET IN X DIRECTION
!
  140 AM(1) = 0.0
  BM(1) = BM(1)+2.*ALPHA*DLX*AX1
  CM(1) = CM(1)+AX1
  CM(MIT) = 0.0
  go to 170
!
!     DIRICHLET-MIXED IN X DIRECTION
!
  150 AM(1) = 0.0
  AM(MIT) = AM(MIT)+CXM
  BM(MIT) = BM(MIT)-2.*BETA*DLX*CXM
  CM(MIT) = 0.0
  go to 170
!
!     MIXED - MIXED IN X DIRECTION
!
  160 CONTINUE
  AM(1) = 0.0
  BM(1) = BM(1)+2.*DLX*ALPHA*AX1
  CM(1) = CM(1)+AX1
  AM(MIT) = AM(MIT)+CXM
  BM(MIT) = BM(MIT)-2.*DLX*BETA*CXM
  CM(MIT) = 0.0
  170 CONTINUE
!
!     ADJUST IN Y DIRECTION UNLESS PERIODIC
!
  DY1 = AN(1)
  FYN = CN(NIT)
  go to (220,180,200,210,190),KSWY
!
!     DIRICHLET-DIRICHLET IN Y DIRECTION
!
  180 CONTINUE
  AN(1) = 0.0
  CN(NIT) = 0.0
  go to 220
!
!     MIXED-DIRICHLET IN Y DIRECTION
!
  190 CONTINUE
  AN(1) = 0.0
  BN(1) = BN(1)+2.*DLY*GAMA*DY1
  CN(1) = CN(1)+DY1
  CN(NIT) = 0.0
  go to 220
!
!     DIRICHLET-MIXED IN Y DIRECTION
!
  200 AN(1) = 0.0
  AN(NIT) = AN(NIT)+FYN
  BN(NIT) = BN(NIT)-2.*DLY*XNU*FYN
  CN(NIT) = 0.0
  go to 220
!
!     MIXED - MIXED DIRECTION IN Y DIRECTION
!
  210 CONTINUE
  AN(1) = 0.0
  BN(1) = BN(1)+2.*DLY*GAMA*DY1
  CN(1) = CN(1)+DY1
  AN(NIT) = AN(NIT)+FYN
  BN(NIT) = BN(NIT)-2.0*DLY*XNU*FYN
  CN(NIT) = 0.0
  220 if (KSWX  ==  1) go to 270
!
!     ADJUST USOL ALONG X EDGE
!
  DO 260 J=JS,NS
     if (KSWX /= 2 .AND. KSWX /= 3) go to 230
     USOL(IS,J) = USOL(IS,J)-AX1*USOL(1,J)
     go to 240
  230    USOL(IS,J) = USOL(IS,J)+2.0*DLX*AX1*BDA(J)
  240    if (KSWX /= 2 .AND. KSWX /= 5) go to 250
     USOL(MS,J) = USOL(MS,J)-CXM*USOL(K,J)
     go to 260
  250    USOL(MS,J) = USOL(MS,J)-2.0*DLX*CXM*BDB(J)
  260 CONTINUE
  270 if (KSWY  ==  1) go to 320
!
!     ADJUST USOL ALONG Y EDGE
!
  DO 310 I=IS,MS
     if (KSWY /= 2 .AND. KSWY /= 3) go to 280
     USOL(I,JS) = USOL(I,JS)-DY1*USOL(I,1)
     go to 290
  280    USOL(I,JS) = USOL(I,JS)+2.0*DLY*DY1*BDC(I)
  290    if (KSWY /= 2 .AND. KSWY /= 5) go to 300
     USOL(I,NS) = USOL(I,NS)-FYN*USOL(I,L)
     go to 310
  300    USOL(I,NS) = USOL(I,NS)-2.0*DLY*FYN*BDD(I)
  310 CONTINUE
  320 CONTINUE
!
!     SAVE ADJUSTED EDGES IN GRHS if IORDER=4
!
  if (IORDER  /=  4) go to 350
  DO 330 J=JS,NS
     GRHS(IS,J) = USOL(IS,J)
     GRHS(MS,J) = USOL(MS,J)
  330 CONTINUE
  DO 340 I=IS,MS
     GRHS(I,JS) = USOL(I,JS)
     GRHS(I,NS) = USOL(I,NS)
  340 CONTINUE
  350 CONTINUE
  IORD = IORDER
  PERTRB = 0.0
!
!     CHECK if OPERATOR IS SINGULAR
!
  call CHKSNG (MBDCND,NBDCND,ALPHA,BETA,GAMA,XNU,COFX,COFY,SINGLR)
!
!     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE
!     if SINGULAR
!
  if (SINGLR) call TRISP (MIT,AM,BM,CM,DM,UM,ZM)
  if (SINGLR) call TRISP (NIT,AN,BN,CN,DN,UN,ZN)
!
!     MAKE INITIALIZATION call TO BLKTRI
!
  if (INTL  ==  0) &
      call BLKTRI (INTL,NP,NIT,AN,BN,CN,MP,MIT,AM,BM,CM,IDMN, &
                   USOL(IS,JS),IERROR,W)
  if (IERROR  /=  0) RETURN
!
!     ADJUST RIGHT HAND SIDE if NECESSARY
!
  360 CONTINUE
  if (SINGLR) call ORTHOG (USOL,IDMN,ZN,ZM,PERTRB)
!
!     COMPUTE SOLUTION
!
  call BLKTRI (I1,NP,NIT,AN,BN,CN,MP,MIT,AM,BM,CM,IDMN,USOL(IS,JS), &
               IERROR,W)
  if (IERROR  /=  0) RETURN
!
!     SET PERIODIC BOUNDARIES if NECESSARY
!
  if (KSWX  /=  1) go to 380
  DO 370 J=1,L
     USOL(K,J) = USOL(1,J)
  370 CONTINUE
  380 if (KSWY  /=  1) go to 400
  DO 390 I=1,K
     USOL(I,L) = USOL(I,1)
  390 CONTINUE
  400 CONTINUE
!
!     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
!     NORM if OPERATOR IS SINGULAR
!
  if (SINGLR) call MINSOL (USOL,IDMN,ZN,ZM,PRTRB)
!
!     return if DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE
!     NOT FLAGGED
!
  if (IORD  ==  2) RETURN
  IORD = 2
!
!     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION
!
  call DEFER (COFX,COFY,IDMN,USOL,GRHS)
  go to 360
end
