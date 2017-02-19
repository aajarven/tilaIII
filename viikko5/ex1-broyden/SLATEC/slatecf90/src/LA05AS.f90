subroutine LA05AS (A, IND, NZ, IA, N, IP, IW, W, G, U)
!
!! LA05AS is subsidiary to SPLP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (LA05AS-S, LA05AD-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM
!     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE
!     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING
!     THE FINAL LETTER =S= IN THE NAMES USED HERE.
!     REVISIONS MADE BY R J HANSON, SNLA, AUGUST, 1979.
!     REVISED SEP. 13, 1979.
!
!     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES
!     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL
!     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN
!     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES
!     SPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED.
!
! IP(I,1),IP(I,2) POINT TO THE START OF ROW/COL I.
! IW(I,1),IW(I,2) HOLD THE NUMBER OF NON-ZEROS IN ROW/COL I.
! DURING THE MAIN BODY OF THIS SUBROUTINE THE VECTORS IW(.,3),IW(.,5),
!     IW(.,7) ARE USED TO HOLD DOUBLY LINKED LISTS OF ROWS THAT HAVE
!     NOT BEEN PIVOTAL AND HAVE EQUAL NUMBERS OF NON-ZEROS.
! IW(.,4),IW(.,6),IW(.,8) HOLD SIMILAR LISTS FOR THE COLUMNS.
! IW(I,3),IW(I,4) HOLD FIRST ROW/COLUMN TO HAVE I NON-ZEROS
!     OR ZERO if THERE ARE NONE.
! IW(I,5), IW(I,6) HOLD ROW/COL NUMBER OF ROW/COL PRIOR TO ROW/COL I
!     IN ITS LIST, OR ZERO if NONE.
! IW(I,7), IW(I,8) HOLD ROW/COL NUMBER OF ROW/COL AFTER ROW/COL I
!     IN ITS LIST, OR ZERO if NONE.
! FOR ROWS/COLS THAT HAVE BEEN PIVOTAL IW(I,5),IW(I,6) HOLD NEGATION OF
!     POSITION OF ROW/COL I IN THE PIVOTAL ORDERING.
!
!***SEE ALSO  SPLP
!***ROUTINES CALLED  LA05ES, MC20AS, R1MACH, XERMSG, XSETUN
!***COMMON BLOCKS    LA05DS
!***REVISION HISTORY  (YYMMDD)
!   811215  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890605  Corrected references to XERRWV.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900402  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!***END PROLOGUE  LA05AS
  INTEGER IP(N,2)
  INTEGER IND(IA,2), IW(N,8)
  REAL A(*), AMAX, AU, AM, G, U, SMALL, W(*)
  LOGICAL FIRST
  CHARACTER*8 XERN0, XERN1, XERN2
!
  COMMON /LA05DS/ SMALL, LP, LENL, LENU, NCP, LROW, LCOL
! EPS IS THE RELATIVE ACCURACY OF FLOATING-POINT COMPUTATION
  SAVE EPS, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  LA05AS
  if (FIRST) THEN
     EPS = 2.0E0 * R1MACH(4)
  end if
  FIRST = .FALSE.
!
!     SET THE OUTPUT UNIT NUMBER FOR THE ERROR PROCESSOR.
!     THE USAGE OF THIS ERROR PROCESSOR IS DOCUMENTED IN THE
!     SANDIA LABS. TECH. REPT. SAND78-1189, BY R E JONES.
  call XSETUN(LP)
  if (U > 1.0E0) U = 1.0E0
  if (U < EPS) U = EPS
  if (N < 1) go to 670
  G = 0.
  DO 50 I=1,N
     W(I) = 0.
     DO 40 J=1,5
        IW(I,J) = 0
   40    CONTINUE
   50 CONTINUE
!
! FLUSH OUT SMALL ENTRIES, COUNT ELEMENTS IN ROWS AND COLUMNS
  L = 1
  LENU = NZ
  DO 80 IDUMMY=1,NZ
     if (L > LENU) go to 90
     DO 60 K=L,LENU
        if (ABS(A(K)) <= SMALL) go to 70
        I = IND(K,1)
        J = IND(K,2)
        G = MAX(ABS(A(K)),G)
        if (I < 1 .OR. I > N) go to 680
        if (J < 1 .OR. J > N) go to 680
        IW(I,1) = IW(I,1) + 1
        IW(J,2) = IW(J,2) + 1
   60    CONTINUE
     go to 90
   70    L = K
     A(L) = A(LENU)
     IND(L,1) = IND(LENU,1)
     IND(L,2) = IND(LENU,2)
     LENU = LENU - 1
   80 CONTINUE
!
   90 LENL = 0
  LROW = LENU
  LCOL = LROW
! MCP IS THE MAXIMUM NUMBER OF COMPRESSES PERMITTED BEFORE AN
!     ERROR RETURN RESULTS.
  MCP = MAX(N/10,20)
  NCP = 0
! CHECK FOR NULL ROW OR COLUMN AND INITIALIZE IP(I,2) TO POINT
!     JUST BEYOND WHERE THE LAST COMPONENT OF COLUMN I OF A WILL
!     BE STORED.
  K = 1
  DO 110 IR=1,N
     K = K + IW(IR,2)
     IP(IR,2) = K
     DO 100 L=1,2
        if (IW(IR,L) <= 0) go to 700
  100    CONTINUE
  110 CONTINUE
! REORDER BY ROWS
! CHECK FOR DOUBLE ENTRIES WHILE USING THE NEWLY CONSTRUCTED
!     ROW FILE TO CONSTRUCT THE COLUMN FILE. NOTE THAT BY PUTTING
!    THE ENTRIES IN BACKWARDS AND DECREASING IP(J,2) EACH TIME IT
!     IS USED WE AUTOMATICALLY LEAVE IT POINTING TO THE FIRST ELEMENT.
  call MC20AS(N, LENU, A, IND(1,2), IP, IND(1,1), 0)
  KL = LENU
  DO 130 II=1,N
     IR = N + 1 - II
     KP = IP(IR,1)
     DO 120 K=KP,KL
        J = IND(K,2)
        if (IW(J,5) == IR) go to 660
        IW(J,5) = IR
        KR = IP(J,2) - 1
        IP(J,2) = KR
        IND(KR,1) = IR
  120    CONTINUE
     KL = KP - 1
  130 CONTINUE
!
! SET UP LINKED LISTS OF ROWS AND COLS WITH EQUAL NUMBERS OF NON-ZEROS.
  DO 150 L=1,2
     DO 140 I=1,N
        NZ = IW(I,L)
        IN = IW(NZ,L+2)
        IW(NZ,L+2) = I
        IW(I,L+6) = IN
        IW(I,L+4) = 0
        if (IN /= 0) IW(IN,L+4) = I
  140    CONTINUE
  150 CONTINUE
!
!
! START OF MAIN ELIMINATION LOOP.
  DO 590 IPV=1,N
! FIND PIVOT. JCOST IS MARKOWITZ COST OF CHEAPEST PIVOT FOUND SO FAR,
!     WHICH IS IN ROW IPP AND COLUMN JP.
     JCOST = N*N
! LOOP ON LENGTH OF COLUMN TO BE SEARCHED
     DO 240 NZ=1,N
        if (JCOST <= (NZ-1)**2) go to 250
        J = IW(NZ,4)
! SEARCH COLUMNS WITH NZ NON-ZEROS.
        DO 190 IDUMMY=1,N
           if (J <= 0) go to 200
           KP = IP(J,2)
           KL = KP + IW(J,2) - 1
           DO 180 K=KP,KL
              I = IND(K,1)
              KCOST = (NZ-1)*(IW(I,1)-1)
              if (KCOST >= JCOST) go to 180
              if (NZ == 1) go to 170
! FIND LARGEST ELEMENT IN ROW OF POTENTIAL PIVOT.
              AMAX = 0.
              K1 = IP(I,1)
              K2 = IW(I,1) + K1 - 1
              DO 160 KK=K1,K2
                 AMAX = MAX(AMAX,ABS(A(KK)))
                 if (IND(KK,2) == J) KJ = KK
  160             CONTINUE
! PERFORM STABILITY TEST.
              if (ABS(A(KJ)) < AMAX*U) go to 180
  170             JCOST = KCOST
              IPP = I
              JP = J
              if (JCOST <= (NZ-1)**2) go to 250
  180          CONTINUE
           J = IW(J,8)
  190       CONTINUE
! SEARCH ROWS WITH NZ NON-ZEROS.
  200       I = IW(NZ,3)
        DO 230 IDUMMY=1,N
           if (I <= 0) go to 240
           AMAX = 0.
           KP = IP(I,1)
           KL = KP + IW(I,1) - 1
! FIND LARGEST ELEMENT IN THE ROW
           DO 210 K=KP,KL
              AMAX = MAX(ABS(A(K)),AMAX)
  210          CONTINUE
           AU = AMAX*U
           DO 220 K=KP,KL
! PERFORM STABILITY TEST.
              if (ABS(A(K)) < AU) go to 220
              J = IND(K,2)
              KCOST = (NZ-1)*(IW(J,2)-1)
              if (KCOST >= JCOST) go to 220
              JCOST = KCOST
              IPP = I
              JP = J
              if (JCOST <= (NZ-1)**2) go to 250
  220          CONTINUE
           I = IW(I,7)
  230       CONTINUE
  240    CONTINUE
!
! PIVOT FOUND.
! REMOVE ROWS AND COLUMNS INVOLVED IN ELIMINATION FROM ORDERING VECTORS.
  250    KP = IP(JP,2)
     KL = IW(JP,2) + KP - 1
     DO 290 L=1,2
        DO 280 K=KP,KL
           I = IND(K,L)
           IL = IW(I,L+4)
           IN = IW(I,L+6)
           if (IL == 0) go to 260
           IW(IL,L+6) = IN
           go to 270
  260          NZ = IW(I,L)
           IW(NZ,L+2) = IN
  270          if (IN > 0) IW(IN,L+4) = IL
  280       CONTINUE
        KP = IP(IPP,1)
        KL = KP + IW(IPP,1) - 1
  290    CONTINUE
! STORE PIVOT
     IW(IPP,5) = -IPV
     IW(JP,6) = -IPV
! ELIMINATE PIVOTAL ROW FROM COLUMN FILE AND FIND PIVOT IN ROW FILE.
     DO 320 K=KP,KL
        J = IND(K,2)
        KPC = IP(J,2)
        IW(J,2) = IW(J,2) - 1
        KLC = KPC + IW(J,2)
        DO 300 KC=KPC,KLC
           if (IPP == IND(KC,1)) go to 310
  300       CONTINUE
  310       IND(KC,1) = IND(KLC,1)
        IND(KLC,1) = 0
        if (J == JP) KR = K
  320    CONTINUE
! BRING PIVOT TO FRONT OF PIVOTAL ROW.
     AU = A(KR)
     A(KR) = A(KP)
     A(KP) = AU
     IND(KR,2) = IND(KP,2)
     IND(KP,2) = JP
!
! PERFORM ELIMINATION ITSELF, LOOPING ON NON-ZEROS IN PIVOT COLUMN.
     NZC = IW(JP,2)
     if (NZC == 0) go to 550
     DO 540 NC=1,NZC
        KC = IP(JP,2) + NC - 1
        IR = IND(KC,1)
! SEARCH NON-PIVOT ROW FOR ELEMENT TO BE ELIMINATED.
        KR = IP(IR,1)
        KRL = KR + IW(IR,1) - 1
        DO 330 KNP=KR,KRL
           if (JP == IND(KNP,2)) go to 340
  330       CONTINUE
! BRING ELEMENT TO BE ELIMINATED TO FRONT OF ITS ROW.
  340       AM = A(KNP)
        A(KNP) = A(KR)
        A(KR) = AM
        IND(KNP,2) = IND(KR,2)
        IND(KR,2) = JP
        AM = -A(KR)/A(KP)
! COMPRESS ROW FILE UNLESS IT IS CERTAIN THAT THERE IS ROOM FOR NEW ROW.
        if (LROW+IW(IR,1)+IW(IPP,1)+LENL <= IA) go to 350
        if (NCP >= MCP .OR. LENU+IW(IR,1)+IW(IPP,1)+LENL > IA) GO &
         TO 710
        call LA05ES(A, IND(1,2), IP, N, IW, IA, .TRUE.)
        KP = IP(IPP,1)
        KR = IP(IR,1)
  350       KRL = KR + IW(IR,1) - 1
        KQ = KP + 1
        KPL = KP + IW(IPP,1) - 1
! PLACE PIVOT ROW (EXCLUDING PIVOT ITSELF) IN W.
        if (KQ > KPL) go to 370
        DO 360 K=KQ,KPL
           J = IND(K,2)
           W(J) = A(K)
  360       CONTINUE
  370       IP(IR,1) = LROW + 1
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
  380          LENU = LENU - 1
! REMOVE ELEMENT FROM COL FILE.
           K = IP(J,2)
           KL = K + IW(J,2) - 1
           IW(J,2) = KL - K
           DO 390 KK=K,KL
              if (IND(KK,1) == IR) go to 400
  390          CONTINUE
  400          IND(KK,1) = IND(KL,1)
           IND(KL,1) = 0
  410          W(J) = 0.
  420       CONTINUE
!
! SCAN PIVOT ROW FOR FILLS.
  430       if (KQ > KPL) go to 520
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
           if (NZ  ==  0) go to 460
! if POSSIBLE PLACE NEW ELEMENT AT END OF PRESENT ENTRY.
           if (KL /= LCOL) go to 440
           if (LCOL+LENL >= IA) go to 460
           LCOL = LCOL + 1
           go to 450
  440          if (IND(KL+1,1) /= 0) go to 460
  450          IND(KL+1,1) = IR
           go to 490
! NEW ENTRY HAS TO BE CREATED.
  460          if (LCOL+LENL+NZ+1 < IA) go to 470
! COMPRESS COLUMN FILE if THERE IS NOT ROOM FOR NEW ENTRY.
           if (NCP >= MCP .OR. LENU+LENL+NZ+1 >= IA) go to 710
           call LA05ES(A, IND, IP(1,2), N, IW(1,2), IA, .FALSE.)
           K = IP(J,2)
           KL = K + NZ - 1
! TRANSFER OLD ENTRY INTO NEW.
  470          IP(J,2) = LCOL + 1
           if (KL  <  K) go to 485
           DO 480 KK=K,KL
              LCOL = LCOL + 1
              IND(LCOL,1) = IND(KK,1)
              IND(KK,1) = 0
  480          CONTINUE
  485          CONTINUE
! ADD NEW ELEMENT.
           LCOL = LCOL + 1
           IND(LCOL,1) = IR
  490          G = MAX(G,ABS(AU))
           IW(J,2) = NZ + 1
  500          W(J) = 0.
  510       CONTINUE
  520       IW(IR,1) = LROW + 1 - IP(IR,1)
!
! STORE MULTIPLIER
        if (LENL+LCOL+1 <= IA) go to 530
! COMPRESS COL FILE if NECESSARY.
        if (NCP >= MCP) go to 710
        call LA05ES(A, IND, IP(1,2), N, IW(1,2), IA, .FALSE.)
  530       K = IA - LENL
        LENL = LENL + 1
        A(K) = AM
        IND(K,1) = IPP
        IND(K,2) = IR
        LENU = LENU - 1
  540    CONTINUE
!
! INSERT ROWS AND COLUMNS INVOLVED IN ELIMINATION IN LINKED LISTS
!     OF EQUAL NUMBERS OF NON-ZEROS.
  550    K1 = IP(JP,2)
     K2 = IW(JP,2) + K1 - 1
     IW(JP,2) = 0
     DO 580 L=1,2
        if (K2 < K1) go to 570
        DO 560 K=K1,K2
           IR = IND(K,L)
           if (L == 1) IND(K,L) = 0
           NZ = IW(IR,L)
           if (NZ <= 0) go to 720
           IN = IW(NZ,L+2)
           IW(IR,L+6) = IN
           IW(IR,L+4) = 0
           IW(NZ,L+2) = IR
           if (IN /= 0) IW(IN,L+4) = IR
  560       CONTINUE
  570       K1 = IP(IPP,1) + 1
        K2 = IW(IPP,1) + K1 - 2
  580    CONTINUE
  590 CONTINUE
!
! RESET COLUMN FILE TO REFER TO U AND STORE ROW/COL NUMBERS IN
!     PIVOTAL ORDER IN IW(.,3),IW(.,4)
  DO 600 I=1,N
     J = -IW(I,5)
     IW(J,3) = I
     J = -IW(I,6)
     IW(J,4) = I
     IW(I,2) = 0
  600 CONTINUE
  DO 620 I=1,N
     KP = IP(I,1)
     KL = IW(I,1) + KP - 1
     DO 610 K=KP,KL
        J = IND(K,2)
        IW(J,2) = IW(J,2) + 1
  610    CONTINUE
  620 CONTINUE
  K = 1
  DO 630 I=1,N
     K = K + IW(I,2)
     IP(I,2) = K
  630 CONTINUE
  LCOL = K - 1
  DO 650 II=1,N
     I = IW(II,3)
     KP = IP(I,1)
     KL = IW(I,1) + KP - 1
     DO 640 K=KP,KL
        J = IND(K,2)
        KN = IP(J,2) - 1
        IP(J,2) = KN
        IND(KN,1) = I
  640    CONTINUE
  650 CONTINUE
  return
!
!     THE FOLLOWING INSTRUCTIONS IMPLEMENT THE FAILURE EXITS.
!
  660 if (LP > 0) THEN
     WRITE (XERN1, '(I8)') IR
     WRITE (XERN2, '(I8)') J
     call XERMSG ('SLATEC', 'LA05AS', 'MORE THAN ONE MATRIX ' // &
        'ENTRY.  HERE ROW = ' // XERN1 // ' AND COL = ' // XERN2, &
        -4, 1)
  end if
  G = -4.
  return
!
  670 if (LP > 0) call XERMSG ('SLATEC', 'LA05AS', &
     'THE ORDER OF THE SYSTEM, N, IS NOT POSITIVE.', -1, 1)
  G = -1.0E0
  return
!
  680 if (LP > 0) THEN
     WRITE (XERN0, '(I8)') K
     WRITE (XERN1, '(I8)') I
     WRITE (XERN2, '(I8)') J
     call XERMSG ('SLATEC', 'LA05AS', 'ELEMENT K = ' // XERN0 // &
        ' IS OUT OF BOUNDS.$$HERE ROW = ' // XERN1 // &
        ' AND COL = ' // XERN2, -3, 1)
  end if
  G = -3.
  return
!
  700 if (LP > 0) THEN
     WRITE (XERN1, '(I8)') L
     call XERMSG ('SLATEC', 'LA05AS', 'ROW OR COLUMN HAS NO ' // &
        'ELEMENTS.  HERE INDEX = ' // XERN1, -2, 1)
  end if
  G = -2.
  return
!
  710 if (LP > 0) call XERMSG ('SLATEC', 'LA05AS', &
     'LENGTHS OF ARRAYS A(*) AND IND(*,2) ARE TOO SMALL.', -7, 1)
  G = -7.
  return
!
  720 IPV = IPV + 1
  IW(IPV,1) = IR
  DO 730 I=1,N
     II = -IW(I,L+4)
     if (II > 0) IW(II,1) = I
  730 CONTINUE
!
  if (LP > 0) THEN
     XERN1 = 'ROWS'
     if (L == 2) XERN1 = 'COLUMNS'
     call XERMSG ('SLATEC', 'LA05AS', 'DEPENDANT ' // XERN1, -5, 1)
!
  740    WRITE (XERN1, '(I8)') IW(I,1)
     XERN2 = ' '
     if (I+1 <= IPV) WRITE (XERN2, '(I8)') IW(I+1,1)
     call XERMSG ('SLATEC', 'LA05AS', &
        'DEPENDENT VECTOR INDICES ARE ' // XERN1 // ' AND ' // &
        XERN2, -5, 1)
     I = I + 2
     if (I <= IPV) go to 740
  end if
  G = -5.
  return
end
