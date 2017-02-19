subroutine LA05CD (A, IND, IA, N, IP, IW, W, G, U, MM)
!
!! LA05CD is subsidiary to DSPLP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (LA05CS-D, LA05CD-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM
!     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE
!     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING
!     THE FINAL LETTER =D= IN THE NAMES USED HERE.
!     REVISED SEP. 13, 1979.
!
!     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES
!     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL
!     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN
!     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES
!     DSPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED.
!
!***SEE ALSO  DSPLP
!***ROUTINES CALLED  LA05ED, XERMSG, XSETUN
!***COMMON BLOCKS    LA05DD
!***REVISION HISTORY  (YYMMDD)
!   811215  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900402  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!   920410  Corrected second dimension on IW declaration.  (WRB)
!   920422  Changed upper limit on DO from LAST to LAST-1.  (WRB)
!***END PROLOGUE  LA05CD
  DOUBLE PRECISION A(*), G, U, AM, W(*), SMALL, AU
  INTEGER IND(IA,2), IW(N,8)
  INTEGER IP(N,2)
  CHARACTER*8 XERN1
!
  COMMON /LA05DD/ SMALL, LP, LENL, LENU, NCP, LROW, LCOL
!***FIRST EXECUTABLE STATEMENT  LA05CD
  call XSETUN(LP)
  if (G < 0.0D0) go to 620
  JM = MM
! MCP LIMITS THE VALUE OF NCP PERMITTED BEFORE AN ERROR RETURN RESULTS.
  MCP = NCP + 20
! REMOVE OLD COLUMN
  LENU = LENU - IW(JM,2)
  KP = IP(JM,2)
  IM = IND(KP,1)
  KL = KP + IW(JM,2) - 1
  IW(JM,2) = 0
  DO 30 K=KP,KL
     I = IND(K,1)
     IND(K,1) = 0
     KR = IP(I,1)
     NZ = IW(I,1) - 1
     IW(I,1) = NZ
     KRL = KR + NZ
     DO 10 KM=KR,KRL
        if (IND(KM,2) == JM) go to 20
   10    CONTINUE
   20    A(KM) = A(KRL)
     IND(KM,2) = IND(KRL,2)
     IND(KRL,2) = 0
   30 CONTINUE
!
! INSERT NEW COLUMN
  DO 110 II=1,N
     I = IW(II,3)
     if (I == IM) M = II
     if (ABS(W(I)) <= SMALL) go to 100
     LENU = LENU + 1
     LAST = II
     if (LCOL+LENL < IA) go to 40
! COMPRESS COLUMN FILE if NECESSARY.
     if (NCP >= MCP .OR. LENL+LENU >= IA) go to 610
     call LA05ED(A, IND, IP(1,2), N, IW(1,2), IA, .FALSE.)
   40    LCOL = LCOL + 1
     NZ = IW(JM,2)
     if (NZ == 0) IP(JM,2) = LCOL
     IW(JM,2) = NZ + 1
     IND(LCOL,1) = I
     NZ = IW(I,1)
     KPL = IP(I,1) + NZ
     if (KPL > LROW) go to 50
     if (IND(KPL,2) == 0) go to 90
! NEW ENTRY HAS TO BE CREATED.
   50    if (LENL+LROW+NZ < IA) go to 60
     if (NCP >= MCP .OR. LENL+LENU+NZ >= IA) go to 610
! COMPRESS ROW FILE if NECESSARY.
     call LA05ED(A, IND(1,2), IP, N, IW, IA, .TRUE.)
   60    KP = IP(I,1)
     IP(I,1) = LROW + 1
     if (NZ == 0) go to 80
     KPL = KP + NZ - 1
     DO 70 K=KP,KPL
        LROW = LROW + 1
        A(LROW) = A(K)
        IND(LROW,2) = IND(K,2)
        IND(K,2) = 0
   70    CONTINUE
   80    LROW = LROW + 1
     KPL = LROW
! PLACE NEW ELEMENT AT END OF ROW.
   90    IW(I,1) = NZ + 1
     A(KPL) = W(I)
     IND(KPL,2) = JM
  100    W(I) = 0.0D0
  110 CONTINUE
  if (IW(IM,1) == 0 .OR. IW(JM,2) == 0 .OR. M > LAST) go to 590
!
! FIND COLUMN SINGLETONS, OTHER THAN THE SPIKE. NON-SINGLETONS ARE
!     MARKED WITH W(J)=1. ONLY IW(.,3) IS REVISED AND IW(.,4) IS USED
!     FOR WORKSPACE.
  INS = M
  M1 = M
  W(JM) = 1.0D0
  DO 140 II=M,LAST
     I = IW(II,3)
     J = IW(II,4)
     if (W(J) == 0.) go to 130
     KP = IP(I,1)
     KL = KP + IW(I,1) - 1
     DO 120 K=KP,KL
        J = IND(K,2)
        W(J) = 1.0D0
  120    CONTINUE
     IW(INS,4) = I
     INS = INS + 1
     go to 140
! PLACE SINGLETONS IN NEW POSITION.
  130    IW(M1,3) = I
     M1 = M1 + 1
  140 CONTINUE
! PLACE NON-SINGLETONS IN NEW POSITION.
  IJ = M + 1
  DO 150 II=M1,LAST-1
     IW(II,3) = IW(IJ,4)
     IJ = IJ + 1
  150 CONTINUE
! PLACE SPIKE AT END.
  IW(LAST,3) = IM
!
! FIND ROW SINGLETONS, APART FROM SPIKE ROW. NON-SINGLETONS ARE MARKED
!     WITH W(I)=2. AGAIN ONLY IW(.,3) IS REVISED AND IW(.,4) IS USED
!     FOR WORKSPACE.
  LAST1 = LAST
  JNS = LAST
  W(IM) = 2.0D0
  J = JM
  DO 180 IJ=M1,LAST
     II = LAST + M1 - IJ
     I = IW(II,3)
     if (W(I) /= 2.0D0) go to 170
     K = IP(I,1)
     if (II /= LAST) J = IND(K,2)
     KP = IP(J,2)
     KL = KP + IW(J,2) - 1
     IW(JNS,4) = I
     JNS = JNS - 1
     DO 160 K=KP,KL
        I = IND(K,1)
        W(I) = 2.0D0
  160    CONTINUE
     go to 180
  170    IW(LAST1,3) = I
     LAST1 = LAST1 - 1
  180 CONTINUE
  DO 190 II=M1,LAST1
     JNS = JNS + 1
     I = IW(JNS,4)
     W(I) = 3.0D0
     IW(II,3) = I
  190 CONTINUE
!
! DEAL WITH SINGLETON SPIKE COLUMN. NOTE THAT BUMP ROWS ARE MARKED BY
!    W(I)=3.
  DO 230 II=M1,LAST1
     KP = IP(JM,2)
     KL = KP + IW(JM,2) - 1
     IS = 0
     DO 200 K=KP,KL
        L = IND(K,1)
        if (W(L) /= 3.0D0) go to 200
        if (IS /= 0) go to 240
        I = L
        KNP = K
        IS = 1
  200    CONTINUE
     if (IS == 0) go to 590
! MAKE A(I,JM) A PIVOT.
     IND(KNP,1) = IND(KP,1)
     IND(KP,1) = I
     KP = IP(I,1)
     DO 210 K=KP,IA
        if (IND(K,2) == JM) go to 220
  210    CONTINUE
  220    AM = A(KP)
     A(KP) = A(K)
     A(K) = AM
     IND(K,2) = IND(KP,2)
     IND(KP,2) = JM
     JM = IND(K,2)
     IW(II,4) = I
     W(I) = 2.0D0
  230 CONTINUE
  II = LAST1
  go to 260
  240 IN = M1
  DO 250 IJ=II,LAST1
     IW(IJ,4) = IW(IN,3)
     IN = IN + 1
  250 CONTINUE
  260 LAST2 = LAST1 - 1
  if (M1 == LAST1) go to 570
  DO 270 I=M1,LAST2
     IW(I,3) = IW(I,4)
  270 CONTINUE
  M1 = II
  if (M1 == LAST1) go to 570
!
! CLEAR W
  DO 280 I=1,N
     W(I) = 0.0D0
  280 CONTINUE
!
! PERFORM ELIMINATION
  IR = IW(LAST1,3)
  DO 560 II=M1,LAST1
     IPP = IW(II,3)
     KP = IP(IPP,1)
     KR = IP(IR,1)
     JP = IND(KP,2)
     if (II == LAST1) JP = JM
! SEARCH NON-PIVOT ROW FOR ELEMENT TO BE ELIMINATED.
!  AND BRING IT TO FRONT OF ITS ROW
     KRL = KR + IW(IR,1) - 1
     DO 290 KNP=KR,KRL
        if (JP == IND(KNP,2)) go to 300
  290    CONTINUE
     if (II-LAST1) 560, 590, 560
! BRING ELEMENT TO BE ELIMINATED TO FRONT OF ITS ROW.
  300    AM = A(KNP)
     A(KNP) = A(KR)
     A(KR) = AM
     IND(KNP,2) = IND(KR,2)
     IND(KR,2) = JP
     if (II == LAST1) go to 310
     if (ABS(A(KP)) < U*ABS(AM)) go to 310
     if (ABS(AM) < U*ABS(A(KP))) go to 340
     if (IW(IPP,1) <= IW(IR,1)) go to 340
! PERFORM INTERCHANGE
  310    IW(LAST1,3) = IPP
     IW(II,3) = IR
     IR = IPP
     IPP = IW(II,3)
     K = KR
     KR = KP
     KP = K
     KJ = IP(JP,2)
     DO 320 K=KJ,IA
        if (IND(K,1) == IPP) go to 330
  320    CONTINUE
  330    IND(K,1) = IND(KJ,1)
     IND(KJ,1) = IPP
  340    if (A(KP) == 0.0D0) go to 590
     if (II == LAST1) go to 560
     AM = -A(KR)/A(KP)
! COMPRESS ROW FILE UNLESS IT IS CERTAIN THAT THERE IS ROOM FOR NEW ROW.
     if (LROW+IW(IR,1)+IW(IPP,1)+LENL <= IA) go to 350
     if (NCP >= MCP .OR. LENU+IW(IR,1)+IW(IPP,1)+LENL > IA) go to &
      610
     call LA05ED(A, IND(1,2), IP, N, IW, IA, .TRUE.)
     KP = IP(IPP,1)
     KR = IP(IR,1)
  350    KRL = KR + IW(IR,1) - 1
     KQ = KP + 1
     KPL = KP + IW(IPP,1) - 1
! PLACE PIVOT ROW (EXCLUDING PIVOT ITSELF) IN W.
     if (KQ > KPL) go to 370
     DO 360 K=KQ,KPL
        J = IND(K,2)
        W(J) = A(K)
  360    CONTINUE
  370    IP(IR,1) = LROW + 1
!
! TRANSFER MODIFIED ELEMENTS.
     IND(KR,2) = 0
     KR = KR + 1
     if (KR > KRL) go to 430
     DO 420 KS=KR,KRL
        J = IND(KS,2)
        AU = A(KS) + AM*W(J)
        IND(KS,2) = 0
! if ELEMENT IS VERY SMALL REMOVE IT FROM U.
        if (ABS(AU) <= SMALL) go to 380
        G = MAX(G,ABS(AU))
        LROW = LROW + 1
        A(LROW) = AU
        IND(LROW,2) = J
        go to 410
  380       LENU = LENU - 1
! REMOVE ELEMENT FROM COL FILE.
        K = IP(J,2)
        KL = K + IW(J,2) - 1
        IW(J,2) = KL - K
        DO 390 KK=K,KL
           if (IND(KK,1) == IR) go to 400
  390       CONTINUE
  400       IND(KK,1) = IND(KL,1)
        IND(KL,1) = 0
  410       W(J) = 0.0D0
  420    CONTINUE
!
! SCAN PIVOT ROW FOR FILLS.
  430    if (KQ > KPL) go to 520
     DO 510 KS=KQ,KPL
        J = IND(KS,2)
        AU = AM*W(J)
        if (ABS(AU) <= SMALL) go to 500
        LROW = LROW + 1
        A(LROW) = AU
        IND(LROW,2) = J
        LENU = LENU + 1
!
! CREATE FILL IN COLUMN FILE.
        NZ = IW(J,2)
        K = IP(J,2)
        KL = K + NZ - 1
! if POSSIBLE PLACE NEW ELEMENT AT END OF PRESENT ENTRY.
        if (KL /= LCOL) go to 440
        if (LCOL+LENL >= IA) go to 460
        LCOL = LCOL + 1
        go to 450
  440       if (IND(KL+1,1) /= 0) go to 460
  450       IND(KL+1,1) = IR
        go to 490
! NEW ENTRY HAS TO BE CREATED.
  460       if (LCOL+LENL+NZ+1 < IA) go to 470
! COMPRESS COLUMN FILE if THERE IS NOT ROOM FOR NEW ENTRY.
        if (NCP >= MCP .OR. LENU+LENL+NZ+1 >= IA) go to 610
        call LA05ED(A, IND, IP(1,2), N, IW(1,2), IA, .FALSE.)
        K = IP(J,2)
        KL = K + NZ - 1
! TRANSFER OLD ENTRY INTO NEW.
  470       IP(J,2) = LCOL + 1
        DO 480 KK=K,KL
           LCOL = LCOL + 1
           IND(LCOL,1) = IND(KK,1)
           IND(KK,1) = 0
  480       CONTINUE
! ADD NEW ELEMENT.
        LCOL = LCOL + 1
        IND(LCOL,1) = IR
  490       G = MAX(G,ABS(AU))
        IW(J,2) = NZ + 1
  500       W(J) = 0.0D0
  510    CONTINUE
  520    IW(IR,1) = LROW + 1 - IP(IR,1)
!
! STORE MULTIPLIER
     if (LENL+LCOL+1 <= IA) go to 530
! COMPRESS COL FILE if NECESSARY.
     if (NCP >= MCP) go to 610
     call LA05ED(A, IND, IP(1,2), N, IW(1,2), IA, .FALSE.)
  530    K = IA - LENL
     LENL = LENL + 1
     A(K) = AM
     IND(K,1) = IPP
     IND(K,2) = IR
! CREATE BLANK IN PIVOTAL COLUMN.
     KP = IP(JP,2)
     NZ = IW(JP,2) - 1
     KL = KP + NZ
     DO 540 K=KP,KL
        if (IND(K,1) == IR) go to 550
  540    CONTINUE
  550    IND(K,1) = IND(KL,1)
     IW(JP,2) = NZ
     IND(KL,1) = 0
     LENU = LENU - 1
  560 CONTINUE
!
! CONSTRUCT COLUMN PERMUTATION AND STORE IT IN IW(.,4)
  570 DO 580 II=M,LAST
     I = IW(II,3)
     K = IP(I,1)
     J = IND(K,2)
     IW(II,4) = J
  580 CONTINUE
  return
!
!     THE FOLLOWING INSTRUCTIONS IMPLEMENT THE FAILURE EXITS.
!
  590 if (LP > 0) THEN
     WRITE (XERN1, '(I8)') MM
     call XERMSG ('SLATEC', 'LA05CD', 'SINGULAR MATRIX AFTER ' // &
        'REPLACEMENT OF COLUMN.  INDEX = ' // XERN1, -6, 1)
  end if
  G = -6.0D0
  return
!
  610 if (LP > 0) call XERMSG ('SLATEC', 'LA05CD', &
     'LENGTHS OF ARRAYS A(*) AND IND(*,2) ARE TOO SMALL.', -7, 1)
  G = -7.0D0
  return
!
  620 if (LP > 0) call XERMSG ('SLATEC', 'LA05CD', &
     'EARLIER ENTRY GAVE ERROR RETURN.', -8, 2)
  G = -8.0D0
  return
end
