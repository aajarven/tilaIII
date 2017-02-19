subroutine LSOD (F, NEQ, T, Y, TOUT, RTOL, ATOL, IDID, YPOUT, YH, &
     YH1, EWT, SAVF, ACOR, WM, IWM, JAC, INTOUT, TSTOP, TOLFAC, &
     DELSGN, RPAR, IPAR)
!
!! LSOD is subsidiary to DEBDF.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (LSOD-S, DLSOD-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   DEBDF  merely allocates storage for  LSOD  to relieve the user of
!   the inconvenience of a long call list.  Consequently  LSOD  is used
!   as described in the comments for  DEBDF .
!
!***SEE ALSO  DEBDF
!***ROUTINES CALLED  HSTART, INTYD, R1MACH, STOD, VNWRMS, XERMSG
!***COMMON BLOCKS    DEBDF1
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!***END PROLOGUE  LSOD
!
!
  LOGICAL INTOUT
!
  DIMENSION Y(*),YPOUT(*),YH(NEQ,6),YH1(*),EWT(*),SAVF(*), &
            ACOR(*),WM(*),IWM(*),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
  CHARACTER*8 XERN1
  CHARACTER*16 XERN3, XERN4
!
  COMMON /DEBDF1/ TOLD, ROWNS(210), &
     EL0, H, HMIN, HMXI, HU, X, U, &
     IQUIT, INIT, LYH, LEWT, LACOR, LSAVF, LWM, KSTEPS, &
     IBEGIN, ITOL, IINTEG, ITSTOP, IJAC, IBAND, IOWNS(6), &
     IER, JSTART, KFLAG, LDUM, METH, MITER, MAXORD, N, NQ, NST, &
     NFE, NJE, NQU
!
  EXTERNAL F, JAC
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
!***FIRST EXECUTABLE STATEMENT  LSOD
  if (IBEGIN  ==  0) THEN
!
!        ON THE FIRST call , PERFORM INITIALIZATION --
!        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
!        FUNCTION ROUTINE R1MACH. THE USER MUST MAKE SURE THAT THE
!        VALUES SET IN R1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
!
     U = R1MACH(4)
!                          -- SET ASSOCIATED MACHINE DEPENDENT PARAMETER
     WM(1) = SQRT(U)
!                          -- SET TERMINATION FLAG
     IQUIT = 0
!                          -- SET INITIALIZATION INDICATOR
     INIT = 0
!                          -- SET COUNTER FOR ATTEMPTED STEPS
     KSTEPS = 0
!                          -- SET INDICATOR FOR INTERMEDIATE-OUTPUT
     INTOUT = .FALSE.
!                          -- SET START INDICATOR FOR STOD CODE
     JSTART = 0
!                          -- SET BDF METHOD INDICATOR
     METH = 2
!                          -- SET MAXIMUM ORDER FOR BDF METHOD
     MAXORD = 5
!                          -- SET ITERATION MATRIX INDICATOR
!
     if (IJAC  ==  0 .AND. IBAND  ==  0) MITER = 2
     if (IJAC  ==  1 .AND. IBAND  ==  0) MITER = 1
     if (IJAC  ==  0 .AND. IBAND  ==  1) MITER = 5
     if (IJAC  ==  1 .AND. IBAND  ==  1) MITER = 4
!
!                          -- SET OTHER NECESSARY ITEMS IN COMMON BLOCK
     N = NEQ
     NST = 0
     NJE = 0
     HMXI = 0.
     NQ = 1
     H = 1.
!                          -- RESET IBEGIN FOR SUBSEQUENT CALLS
     IBEGIN=1
  end if
!
!.......................................................................
!
!      CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY
!
  if (NEQ  <  1) THEN
     WRITE (XERN1, '(I8)') NEQ
     call XERMSG ('SLATEC', 'LSOD', &
        'IN DEBDF, THE NUMBER OF EQUATIONS MUST BE A POSITIVE ' // &
        'INTEGER.$$YOU HAVE CALLED THE CODE WITH NEQ = ' // XERN1, &
        6, 1)
     IDID=-33
  end if
!
  NRTOLP = 0
  NATOLP = 0
  DO 60 K = 1,NEQ
     if (NRTOLP  <=  0) THEN
        if (RTOL(K)  <  0.) THEN
           WRITE (XERN1, '(I8)') K
           WRITE (XERN3, '(1PE15.6)') RTOL(K)
           call XERMSG ('SLATEC', 'LSOD', &
              'IN DEBDF, THE RELATIVE ERROR TOLERANCES MUST ' // &
              'BE NON-NEGATIVE.$$YOU HAVE CALLED THE CODE WITH ' // &
              'RTOL(' // XERN1 // ') = ' // XERN3 // '$$IN THE ' // &
              'CASE OF VECTOR ERROR TOLERANCES, NO FURTHER ' // &
              'CHECKING OF RTOL COMPONENTS IS DONE.', 7, 1)
           IDID = -33
           if (NATOLP  >  0) go to 70
           NRTOLP = 1
        ELSEIF (NATOLP  >  0) THEN
           go to 50
        ENDIF
     ENDIF
!
     if (ATOL(K)  <  0.) THEN
        WRITE (XERN1, '(I8)') K
        WRITE (XERN3, '(1PE15.6)') ATOL(K)
        call XERMSG ('SLATEC', 'LSOD', &
           'IN DEBDF, THE ABSOLUTE ERROR ' // &
           'TOLERANCES MUST BE NON-NEGATIVE.$$YOU HAVE CALLED ' // &
           'THE CODE WITH ATOL(' // XERN1 // ') = ' // XERN3 // &
           '$$IN THE CASE OF VECTOR ERROR TOLERANCES, NO FURTHER ' &
           // 'CHECKING OF ATOL COMPONENTS IS DONE.', 8, 1)
        IDID=-33
        if (NRTOLP  >  0) go to 70
        NATOLP=1
     ENDIF
   50    if (ITOL  ==  0) go to 70
   60 CONTINUE
!
   70 if (ITSTOP  ==  1) THEN
     if (SIGN(1.,TOUT-T)  /=  SIGN(1.,TSTOP-T) .OR. &
        ABS(TOUT-T)  >  ABS(TSTOP-T)) THEN
        WRITE (XERN3, '(1PE15.6)') TOUT
        WRITE (XERN4, '(1PE15.6)') TSTOP
        call XERMSG ('SLATEC', 'LSOD', &
           'IN DEBDF, YOU HAVE CALLED THE ' // &
           'CODE WITH TOUT = ' // XERN3 // '$$BUT YOU HAVE ' // &
           'ALSO TOLD THE CODE NOT TO INTEGRATE PAST THE POINT ' // &
           'TSTOP = ' // XERN4 // ' BY SETTING INFO(4) = 1.  ' // &
           'THESE INSTRUCTIONS CONFLICT.', 14, 1)
        IDID=-33
     ENDIF
  end if
!
!        CHECK SOME CONTINUATION POSSIBILITIES
!
  if (INIT  /=  0) THEN
     if (T  ==  TOUT) THEN
        WRITE (XERN3, '(1PE15.6)') T
        call XERMSG ('SLATEC', 'LSOD', &
           'IN DEBDF, YOU HAVE CALLED THE CODE WITH T = TOUT = ' // &
           XERN3 // '  THIS IS NOT ALLOWED ON CONTINUATION CALLS.', &
           9, 1)
        IDID=-33
     ENDIF
!
     if (T  /=  TOLD) THEN
        WRITE (XERN3, '(1PE15.6)') TOLD
        WRITE (XERN4, '(1PE15.6)') T
        call XERMSG ('SLATEC', 'LSOD', &
           'IN DEBDF, YOU HAVE CHANGED THE VALUE OF T FROM ' // &
           XERN3 // ' TO ' // XERN4 // &
           '  THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 10, 1)
        IDID=-33
     ENDIF
!
     if (INIT  /=  1) THEN
        if (DELSGN*(TOUT-T)  <  0.) THEN
           WRITE (XERN3, '(1PE15.6)') TOUT
           call XERMSG ('SLATEC', 'LSOD', &
              'IN DEBDF, BY CALLING THE CODE WITH TOUT = ' // &
              XERN3 // ' YOU ARE ATTEMPTING TO CHANGE THE ' // &
              'DIRECTION OF INTEGRATION.$$' // &
              'THIS IS NOT ALLOWED WITHOUT RESTARTING.', 11, 1)
           IDID=-33
        ENDIF
     ENDIF
  end if
!
  if (IDID  ==  (-33)) THEN
     if (IQUIT  /=  (-33)) THEN
!                       INVALID INPUT DETECTED
        IQUIT=-33
        IBEGIN=-1
     ELSE
        call XERMSG ('SLATEC', 'LSOD', &
           'IN DEBDF, INVALID INPUT WAS ' // &
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
!     100*U WHICH IS LIKELY TO BE REASONABLE FOR THIS METHOD AND MACHINE
!
  DO 170 K=1,NEQ
    if (RTOL(K)+ATOL(K)  >  0.) go to 160
    RTOL(K)=100.*U
    IDID=-2
  160   if (ITOL  ==  0) go to 180
  170   CONTINUE
!
  180 if (IDID  /=  (-2)) go to 190
!                       RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A
!                                                SMALL POSITIVE VALUE
  IBEGIN=-1
  return
!
!     BRANCH ON STATUS OF INITIALIZATION INDICATOR
!            INIT=0 MEANS INITIAL DERIVATIVES AND NOMINAL STEP SIZE
!                   AND DIRECTION NOT YET SET
!            INIT=1 MEANS NOMINAL STEP SIZE AND DIRECTION NOT YET SET
!            INIT=2 MEANS NO FURTHER INITIALIZATION REQUIRED
!
  190 if (INIT  ==  0) go to 200
  if (INIT  ==  1) go to 220
  go to 240
!
!.......................................................................
!
!     MORE INITIALIZATION --
!                         -- EVALUATE INITIAL DERIVATIVES
!
  200 INIT=1
  call F(T,Y,YH(1,2),RPAR,IPAR)
  NFE=1
  if (T  /=  TOUT) go to 220
  IDID=2
  DO 210 L = 1,NEQ
  210    YPOUT(L) = YH(L,2)
  TOLD=T
  return
!
!                         -- COMPUTE INITIAL STEP SIZE
!                         -- SAVE SIGN OF INTEGRATION DIRECTION
!                         -- SET INDEPENDENT AND DEPENDENT VARIABLES
!                                              X AND YH(*) FOR STOD
!
  220 LTOL = 1
  DO 225 L=1,NEQ
    if (ITOL  ==  1) LTOL = L
    TOL = RTOL(LTOL)*ABS(Y(L)) + ATOL(LTOL)
    if (TOL  ==  0.) go to 380
  225   EWT(L) = TOL
!
  BIG = SQRT(R1MACH(2))
  call HSTART (F,NEQ,T,TOUT,Y,YH(1,2),EWT,1,U,BIG, &
               YH(1,3),YH(1,4),YH(1,5),YH(1,6),RPAR,IPAR,H)
!
  DELSGN = SIGN(1.0,TOUT-T)
  X = T
  DO 230 L = 1,NEQ
    YH(L,1) = Y(L)
  230   YH(L,2) = H*YH(L,2)
  INIT = 2
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
  250 if (ABS(X-T)  <  ABSDEL) go to 270
  call INTYD(TOUT,0,YH,NEQ,Y,INTFLG)
  call INTYD(TOUT,1,YH,NEQ,YPOUT,INTFLG)
  IDID = 3
  if (X  /=  TOUT) go to 260
  IDID = 2
  INTOUT = .FALSE.
  260 T = TOUT
  TOLD = T
  return
!
!   if CANNOT GO PAST TSTOP AND SUFFICIENTLY CLOSE,
!   EXTRAPOLATE AND RETURN
!
  270 if (ITSTOP  /=  1) go to 290
  if (ABS(TSTOP-X)  >=  100.*U*ABS(X)) go to 290
  DT = TOUT - X
  DO 280 L = 1,NEQ
  280   Y(L) = YH(L,1) + (DT/H)*YH(L,2)
  call F(TOUT,Y,YPOUT,RPAR,IPAR)
  NFE = NFE + 1
  IDID = 3
  T = TOUT
  TOLD = T
  return
!
  290 if (IINTEG  ==  0  .OR.  .NOT.INTOUT) go to 300
!
!   INTERMEDIATE-OUTPUT MODE
!
  IDID = 1
  go to 500
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
  IBEGIN = -1
  go to 500
!
!.......................................................................
!
!   LIMIT STEP SIZE AND SET WEIGHT VECTOR
!
  330 HMIN = 100.*U*ABS(X)
  HA = MAX(ABS(H),HMIN)
  if (ITSTOP  /=  1) go to 340
  HA = MIN(HA,ABS(TSTOP-X))
  340 H = SIGN(HA,H)
  LTOL = 1
  DO 350 L = 1,NEQ
    if (ITOL  ==  1) LTOL = L
    EWT(L) = RTOL(LTOL)*ABS(YH(L,1)) + ATOL(LTOL)
    if (EWT(L)  <=  0.0) go to 380
  350   CONTINUE
  TOLFAC = U*VNWRMS(NEQ,YH,EWT)
  if (TOLFAC  <=  1.) go to 400
!
!                       TOLERANCES TOO SMALL
  IDID = -2
  TOLFAC = 2.*TOLFAC
  RTOL(1) = TOLFAC*RTOL(1)
  ATOL(1) = TOLFAC*ATOL(1)
  if (ITOL  ==  0) go to 370
  DO 360 L = 2,NEQ
    RTOL(L) = TOLFAC*RTOL(L)
  360   ATOL(L) = TOLFAC*ATOL(L)
  370 IBEGIN = -1
  go to 500
!
!                       RELATIVE ERROR CRITERION INAPPROPRIATE
  380 IDID = -3
  IBEGIN = -1
  go to 500
!
!.......................................................................
!
!     TAKE A STEP
!
  400 call STOD(NEQ,Y,YH,NEQ,YH1,EWT,SAVF,ACOR,WM,IWM,F,JAC,RPAR,IPAR)
!
  JSTART = -2
  INTOUT = .TRUE.
  if (KFLAG  ==  0) go to 250
!
!.......................................................................
!
  if (KFLAG  ==  -1) go to 450
!
!                       REPEATED CORRECTOR CONVERGENCE FAILURES
  IDID = -6
  IBEGIN = -1
  go to 500
!
!                       REPEATED ERROR TEST FAILURES
  450 IDID = -7
  IBEGIN = -1
!
!.......................................................................
!
!                       STORE VALUES BEFORE RETURNING TO DEBDF
  500 DO 555 L = 1,NEQ
    Y(L) = YH(L,1)
  555   YPOUT(L) = YH(L,2)/H
  T = X
  TOLD = T
  INTOUT = .FALSE.
  return
end
