subroutine DES (F, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID, YPOUT, &
     YP, YY, WT, P, PHI, ALPHA, BETA, PSI, V, W, SIG, G, GI, H, EPS, &
     X, XOLD, HOLD, TOLD, DELSGN, TSTOP, TWOU, FOURU, START, PHASE1, &
     NORND, STIFF, INTOUT, NS, KORD, KOLD, INIT, KSTEPS, KLE4, &
     IQUIT, KPREV, IVC, IV, KGI, RPAR, IPAR)
!
!! DES is subsidiary to DEABM.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (DES-S, DDES-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   DEABM  merely allocates storage for  DES  to relieve the user of the
!   inconvenience of a long call list.  Consequently  DES  is used as
!   described in the comments for  DEABM .
!
!***SEE ALSO  DEABM
!***ROUTINES CALLED  R1MACH, SINTRP, STEPS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800501  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls, replace GOTOs with
!           IF-THEN-ELSEs.  (RWC)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DES
!
  LOGICAL STIFF,CRASH,START,PHASE1,NORND,INTOUT
!
  DIMENSION Y(*),YY(*),WT(*),PHI(NEQ,16),P(*),YP(*), &
    YPOUT(*),PSI(12),ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13), &
    GI(11),IV(10),INFO(15),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
  CHARACTER*8 XERN1
  CHARACTER*16 XERN3, XERN4
!
  EXTERNAL F
!
!.......................................................................
!
!  THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
!  NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE COUNTER
!  IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE EXCESSIVE
!  WORK.
!
  SAVE MAXNUM
  DATA MAXNUM/500/
!
!.......................................................................
!
!***FIRST EXECUTABLE STATEMENT  DES
  if (INFO(1)  ==  0) THEN
!
! ON THE FIRST call , PERFORM INITIALIZATION --
!        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
!        FUNCTION ROUTINE  R1MACH. THE USER MUST MAKE SURE THAT THE
!        VALUES SET IN R1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
!
     U=R1MACH(4)
!                       -- SET ASSOCIATED MACHINE DEPENDENT PARAMETERS
     TWOU=2.*U
     FOURU=4.*U
!                       -- SET TERMINATION FLAG
     IQUIT=0
!                       -- SET INITIALIZATION INDICATOR
     INIT=0
!                       -- SET COUNTER FOR ATTEMPTED STEPS
     KSTEPS=0
!                       -- SET INDICATOR FOR INTERMEDIATE-OUTPUT
     INTOUT= .FALSE.
!                       -- SET INDICATOR FOR STIFFNESS DETECTION
     STIFF= .FALSE.
!                       -- SET STEP COUNTER FOR STIFFNESS DETECTION
     KLE4=0
!                       -- SET INDICATORS FOR STEPS CODE
     START= .TRUE.
     PHASE1= .TRUE.
     NORND= .TRUE.
!                       -- RESET INFO(1) FOR SUBSEQUENT CALLS
     INFO(1)=1
  end if
!
!.......................................................................
!
!      CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY
!
  if (INFO(1)  /=  0  .AND.  INFO(1)  /=  1) THEN
     WRITE (XERN1, '(I8)') INFO(1)
     call XERMSG ('SLATEC', 'DES', &
        'IN DEABM, INFO(1) MUST BE ' // &
        'SET TO 0 FOR THE START OF A NEW PROBLEM, AND MUST BE ' // &
        'SET TO 1 FOLLOWING AN INTERRUPTED TASK.  YOU ARE ' // &
        'ATTEMPTING TO CONTINUE THE INTEGRATION ILLEGALLY BY ' // &
        'CALLING THE CODE WITH INFO(1) = ' // XERN1, 3, 1)
     IDID=-33
  end if
!
  if (INFO(2)  /=  0  .AND.  INFO(2)  /=  1) THEN
     WRITE (XERN1, '(I8)') INFO(2)
     call XERMSG ('SLATEC', 'DES', &
        'IN DEABM, INFO(2) MUST BE 0 OR 1 ' // &
        'INDICATING SCALAR AND VECTOR ERROR TOLERANCES, ' // &
        'RESPECTIVELY.  YOU HAVE CALLED THE CODE WITH INFO(2) = ' // &
        XERN1, 4, 1)
     IDID=-33
  end if
!
  if (INFO(3)  /=  0  .AND.  INFO(3)  /=  1) THEN
     WRITE (XERN1, '(I8)') INFO(3)
     call XERMSG ('SLATEC', 'DES', &
        'IN DEABM, INFO(3) MUST BE 0 OR 1 ' // &
        'INDICATING THE INTERVAL OR INTERMEDIATE-OUTPUT MODE OF ' // &
        'INTEGRATION, RESPECTIVELY.  YOU HAVE CALLED THE CODE ' // &
        'WITH  INFO(3) = ' // XERN1, 5, 1)
     IDID=-33
  end if
!
  if (INFO(4)  /=  0  .AND.  INFO(4)  /=  1) THEN
     WRITE (XERN1, '(I8)') INFO(4)
     call XERMSG ('SLATEC', 'DES', &
        'IN DEABM, INFO(4) MUST BE 0 OR 1 ' // &
        'INDICATING WHETHER OR NOT THE INTEGRATION INTERVAL IS ' // &
        'TO BE RESTRICTED BY A POINT TSTOP.  YOU HAVE CALLED ' // &
        'THE CODE WITH INFO(4) = ' // XERN1, 14, 1)
     IDID=-33
  end if
!
  if (NEQ  <  1) THEN
     WRITE (XERN1, '(I8)') NEQ
     call XERMSG ('SLATEC', 'DES', &
        'IN DEABM, THE NUMBER OF EQUATIONS ' // &
        'NEQ MUST BE A POSITIVE INTEGER.  YOU HAVE CALLED THE ' // &
        'CODE WITH  NEQ = ' // XERN1, 6, 1)
     IDID=-33
  end if
!
  NRTOLP = 0
  NATOLP = 0
  DO 90 K=1,NEQ
     if (NRTOLP  ==  0 .AND. RTOL(K)  <  0.) THEN
        WRITE (XERN1, '(I8)') K
        WRITE (XERN3, '(1PE15.6)') RTOL(K)
        call XERMSG ('SLATEC', 'DES', &
           'IN DEABM, THE RELATIVE ERROR ' // &
           'TOLERANCES RTOL MUST BE NON-NEGATIVE.  YOU HAVE ' // &
           'CALLED THE CODE WITH  RTOL(' // XERN1 // ') = ' // &
           XERN3 // '.  IN THE CASE OF VECTOR ERROR TOLERANCES, ' // &
           'NO FURTHER CHECKING OF RTOL COMPONENTS IS DONE.', 7, 1)
        IDID = -33
        NRTOLP = 1
     ENDIF
!
     if (NATOLP  ==  0 .AND. ATOL(K)  <  0.) THEN
        WRITE (XERN1, '(I8)') K
        WRITE (XERN3, '(1PE15.6)') ATOL(K)
        call XERMSG ('SLATEC', 'DES', &
           'IN DEABM, THE ABSOLUTE ERROR ' // &
           'TOLERANCES ATOL MUST BE NON-NEGATIVE.  YOU HAVE ' // &
           'CALLED THE CODE WITH  ATOL(' // XERN1 // ') = ' // &
           XERN3 // '.  IN THE CASE OF VECTOR ERROR TOLERANCES, ' // &
           'NO FURTHER CHECKING OF ATOL COMPONENTS IS DONE.', 8, 1)
        IDID = -33
        NATOLP = 1
     ENDIF
!
     if (INFO(2)  ==  0) go to 100
     if (NATOLP > 0 .AND. NRTOLP > 0) go to 100
   90 CONTINUE
!
  100 if (INFO(4)  ==  1) THEN
     if (SIGN(1.,TOUT-T)  /=  SIGN(1.,TSTOP-T) &
        .OR. ABS(TOUT-T)  >  ABS(TSTOP-T)) THEN
        WRITE (XERN3, '(1PE15.6)') TOUT
        WRITE (XERN4, '(1PE15.6)') TSTOP
        call XERMSG ('SLATEC', 'DES', &
           'IN DEABM, YOU HAVE CALLED THE ' // &
           'CODE WITH  TOUT = ' // XERN3 // ' BUT YOU HAVE ' // &
           'ALSO TOLD THE CODE (INFO(4) = 1) NOT TO INTEGRATE ' // &
           'PAST THE POINT TSTOP = ' // XERN4 // ' THESE ' // &
           'INSTRUCTIONS CONFLICT.', 14, 1)
        IDID=-33
     ENDIF
  end if
!
!     CHECK SOME CONTINUATION POSSIBILITIES
!
  if (INIT  /=  0) THEN
     if (T  ==  TOUT) THEN
        WRITE (XERN3, '(1PE15.6)') T
        call XERMSG ('SLATEC', 'DES', &
           'IN DEABM, YOU HAVE CALLED THE ' // &
           'CODE WITH  T = TOUT = ' // XERN3 // '$$THIS IS NOT ' // &
           'ALLOWED ON CONTINUATION CALLS.', 9, 1)
        IDID=-33
     ENDIF
!
     if (T  /=  TOLD) THEN
        WRITE (XERN3, '(1PE15.6)') TOLD
        WRITE (XERN4, '(1PE15.6)') T
        call XERMSG ('SLATEC', 'DES', &
           'IN DEABM, YOU HAVE CHANGED THE ' // &
           'VALUE OF T FROM ' // XERN3 // ' TO ' // XERN4 // &
           '  THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 10, 1)
        IDID=-33
     ENDIF
!
     if (INIT  /=  1) THEN
        if (DELSGN*(TOUT-T)  <  0.) THEN
           WRITE (XERN3, '(1PE15.6)') TOUT
           call XERMSG ('SLATEC', 'DES', &
              'IN DEABM, BY CALLING THE ' // &
              'CODE WITH TOUT = ' // XERN3 // ' YOU ARE ' // &
              'ATTEMPTING TO CHANGE THE DIRECTION OF ' // &
              'INTEGRATION.$$THIS IS NOT ALLOWED WITHOUT ' // &
              'RESTARTING.', 11, 1)
           IDID=-33
        ENDIF
     ENDIF
  end if
!
!     INVALID INPUT DETECTED
!
  if (IDID  ==  (-33)) THEN
     if (IQUIT  /=  (-33)) THEN
        IQUIT = -33
        INFO(1) = -1
     ELSE
        call XERMSG ('SLATEC', 'DES', &
           'IN DEABM, INVALID INPUT WAS ' // &
           'DETECTED ON SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE ' // &
           'TO PROCEED BECAUSE YOU HAVE NOT CORRECTED THE ' // &
           'PROBLEM, SO EXECUTION IS BEING TERMINATED.', 12, 2)
     ENDIF
     return
  end if
!
!.......................................................................
!
!     RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED AS
!     ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS CASE,
!     THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE SMALLEST VALUE
!     FOURU WHICH IS LIKELY TO BE REASONABLE FOR THIS METHOD AND MACHINE
!
  DO 180 K=1,NEQ
    if (RTOL(K)+ATOL(K)  >  0.) go to 170
    RTOL(K)=FOURU
    IDID=-2
  170   if (INFO(2)  ==  0) go to 190
  180   CONTINUE
!
  190 if (IDID  /=  (-2)) go to 200
!                       RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A
!                                                SMALL POSITIVE VALUE
  INFO(1)=-1
  return
!
!     BRANCH ON STATUS OF INITIALIZATION INDICATOR
!            INIT=0 MEANS INITIAL DERIVATIVES AND NOMINAL STEP SIZE
!                   AND DIRECTION NOT YET SET
!            INIT=1 MEANS NOMINAL STEP SIZE AND DIRECTION NOT YET SET
!            INIT=2 MEANS NO FURTHER INITIALIZATION REQUIRED
!
  200 if (INIT  ==  0) go to 210
  if (INIT  ==  1) go to 220
  go to 240
!
!.......................................................................
!
!     MORE INITIALIZATION --
!                         -- EVALUATE INITIAL DERIVATIVES
!
  210 INIT=1
  A=T
  call F(A,Y,YP,RPAR,IPAR)
  if (T  /=  TOUT) go to 220
  IDID=2
  DO 215 L = 1,NEQ
  215    YPOUT(L) = YP(L)
  TOLD=T
  return
!
!                         -- SET INDEPENDENT AND DEPENDENT VARIABLES
!                                              X AND YY(*) FOR STEPS
!                         -- SET SIGN OF INTEGRATION DIRECTION
!                         -- INITIALIZE THE STEP SIZE
!
  220 INIT = 2
  X = T
  DO 230 L = 1,NEQ
  230   YY(L) = Y(L)
  DELSGN = SIGN(1.0,TOUT-T)
  H = SIGN(MAX(FOURU*ABS(X),ABS(TOUT-X)),TOUT-X)
!
!.......................................................................
!
!   ON EACH call SET INFORMATION WHICH DETERMINES THE ALLOWED INTERVAL
!   OF INTEGRATION BEFORE RETURNING WITH AN ANSWER AT TOUT
!
  240 DEL = TOUT - T
  ABSDEL = ABS(DEL)
!
!.......................................................................
!
!   if ALREADY PAST OUTPUT POINT, INTERPOLATE AND RETURN
!
  250 if ( ABS(X-T)  <  ABSDEL) go to 260
  call SINTRP(X,YY,TOUT,Y,YPOUT,NEQ,KOLD,PHI,IVC,IV,KGI,GI, &
                                          ALPHA,G,W,XOLD,P)
  IDID = 3
  if (X  /=  TOUT) go to 255
  IDID = 2
  INTOUT = .FALSE.
  255 T = TOUT
  TOLD = T
  return
!
!   if CANNOT GO PAST TSTOP AND SUFFICIENTLY CLOSE,
!   EXTRAPOLATE AND RETURN
!
  260 if (INFO(4)  /=  1) go to 280
  if (ABS(TSTOP-X)  >=  FOURU*ABS(X)) go to 280
  DT = TOUT - X
  DO 270 L = 1,NEQ
  270   Y(L) = YY(L) + DT*YP(L)
  call F(TOUT,Y,YPOUT,RPAR,IPAR)
  IDID = 3
  T = TOUT
  TOLD = T
  return
!
  280 if (INFO(3)  ==  0  .OR.  .NOT.INTOUT) go to 300
!
!   INTERMEDIATE-OUTPUT MODE
!
  IDID = 1
  DO 290 L = 1,NEQ
    Y(L)=YY(L)
  290   YPOUT(L) = YP(L)
  T = X
  TOLD = T
  INTOUT = .FALSE.
  return
!
!.......................................................................
!
!     MONITOR NUMBER OF STEPS ATTEMPTED
!
  300 if (KSTEPS  <=  MAXNUM) go to 330
!
!                       A SIGNIFICANT AMOUNT OF WORK HAS BEEN EXPENDED
  IDID=-1
  KSTEPS=0
  if (.NOT. STIFF) go to 310
!
!                       PROBLEM APPEARS TO BE STIFF
  IDID=-4
  STIFF= .FALSE.
  KLE4=0
!
  310 DO 320 L = 1,NEQ
    Y(L) = YY(L)
  320   YPOUT(L) = YP(L)
  T = X
  TOLD = T
  INFO(1) = -1
  INTOUT = .FALSE.
  return
!
!.......................................................................
!
!   LIMIT STEP SIZE, SET WEIGHT VECTOR AND TAKE A STEP
!
  330 HA = ABS(H)
  if (INFO(4)  /=  1) go to 340
  HA = MIN(HA,ABS(TSTOP-X))
  340 H = SIGN(HA,H)
  EPS = 1.0
  LTOL = 1
  DO 350 L = 1,NEQ
    if (INFO(2)  ==  1) LTOL = L
    WT(L) = RTOL(LTOL)*ABS(YY(L)) + ATOL(LTOL)
    if (WT(L)  <=  0.0) go to 360
  350   CONTINUE
  go to 380
!
!                       RELATIVE ERROR CRITERION INAPPROPRIATE
  360 IDID = -3
  DO 370 L = 1,NEQ
    Y(L) = YY(L)
  370   YPOUT(L) = YP(L)
  T = X
  TOLD = T
  INFO(1) = -1
  INTOUT = .FALSE.
  return
!
  380 call STEPS(F,NEQ,YY,X,H,EPS,WT,START,HOLD,KORD,KOLD,CRASH,PHI,P, &
             YP,PSI,ALPHA,BETA,SIG,V,W,G,PHASE1,NS,NORND,KSTEPS, &
             TWOU,FOURU,XOLD,KPREV,IVC,IV,KGI,GI,RPAR,IPAR)
!
!.......................................................................
!
  if ( .NOT.CRASH) go to 420
!
!                       TOLERANCES TOO SMALL
  IDID = -2
  RTOL(1) = EPS*RTOL(1)
  ATOL(1) = EPS*ATOL(1)
  if (INFO(2)  ==  0) go to 400
  DO 390 L = 2,NEQ
    RTOL(L) = EPS*RTOL(L)
  390   ATOL(L) = EPS*ATOL(L)
  400 DO 410 L = 1,NEQ
    Y(L) = YY(L)
  410   YPOUT(L) = YP(L)
  T = X
  TOLD = T
  INFO(1) = -1
  INTOUT = .FALSE.
  return
!
!   (STIFFNESS TEST) COUNT NUMBER OF CONSECUTIVE STEPS TAKEN WITH THE
!   ORDER OF THE METHOD BEING LESS OR EQUAL TO FOUR
!
  420 KLE4 = KLE4 + 1
  if ( KOLD  >  4) KLE4 = 0
  if ( KLE4  >=  50) STIFF = .TRUE.
  INTOUT = .TRUE.
  go to 250
end
