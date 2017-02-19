subroutine DLSOD (DF, NEQ, T, Y, TOUT, RTOL, ATOL, IDID, YPOUT, &
     YH, YH1, EWT, SAVF, ACOR, WM, IWM, DJAC, INTOUT, TSTOP, TOLFAC, &
     DELSGN, RPAR, IPAR)
!
!! DLSOD is subsidiary to DDEBDF.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (LSOD-S, DLSOD-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   DDEBDF  merely allocates storage for  DLSOD  to relieve the user of
!   the inconvenience of a long call list.  Consequently  DLSOD  is used
!   as described in the comments for  DDEBDF .
!
!***SEE ALSO  DDEBDF
!***ROUTINES CALLED  D1MACH, DHSTRT, DINTYD, DSTOD, DVNRMS, XERMSG
!***COMMON BLOCKS    DDEBD1
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!***END PROLOGUE  DLSOD
!
  INTEGER IBAND, IBEGIN, IDID, IER, IINTEG, IJAC, INIT, INTFLG, &
        IOWNS, IPAR, IQUIT, ITOL, ITSTOP, IWM, JSTART, K, KFLAG, &
        KSTEPS, L, LACOR, LDUM, LEWT, LSAVF, LTOL, LWM, LYH, MAXNUM, &
        MAXORD, METH, MITER, N, NATOLP, NEQ, NFE, NJE, NQ, NQU, &
        NRTOLP, NST
  DOUBLE PRECISION ABSDEL, ACOR, ATOL, BIG, D1MACH, DEL, &
        DELSGN, DT, DVNRMS, EL0, EWT, &
        H, HA, HMIN, HMXI, HU, ROWNS, RPAR, RTOL, SAVF, T, TOL, &
        TOLD, TOLFAC, TOUT, TSTOP, U, WM, X, Y, YH, YH1, YPOUT
  LOGICAL INTOUT
  CHARACTER*8 XERN1
  CHARACTER*16 XERN3, XERN4
!
  DIMENSION Y(*),YPOUT(*),YH(NEQ,6),YH1(*),EWT(*),SAVF(*), &
            ACOR(*),WM(*),IWM(*),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
!
!
  COMMON /DDEBD1/ TOLD,ROWNS(210),EL0,H,HMIN,HMXI,HU,X,U,IQUIT,INIT, &
                  LYH,LEWT,LACOR,LSAVF,LWM,KSTEPS,IBEGIN,ITOL, &
                  IINTEG,ITSTOP,IJAC,IBAND,IOWNS(6),IER,JSTART, &
                  KFLAG,LDUM,METH,MITER,MAXORD,N,NQ,NST,NFE,NJE,NQU
!
  EXTERNAL DF, DJAC
!
!     ..................................................................
!
!       THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
!       NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE
!       COUNTER IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE
!       EXCESSIVE WORK.
  SAVE MAXNUM
!
  DATA MAXNUM /500/
!
!     ..................................................................
!
!***FIRST EXECUTABLE STATEMENT  DLSOD
  if (IBEGIN  ==  0) THEN
!
!        ON THE FIRST call , PERFORM INITIALIZATION --
!        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
!        FUNCTION ROUTINE D1MACH. THE USER MUST MAKE SURE THAT THE
!        VALUES SET IN D1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
!
     U = D1MACH(4)
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
!                          -- SET START INDICATOR FOR DSTOD CODE
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
     HMXI = 0.0D0
     NQ = 1
     H = 1.0D0
!                          -- RESET IBEGIN FOR SUBSEQUENT CALLS
     IBEGIN = 1
  end if
!
!     ..................................................................
!
!      CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY
!
  if (NEQ  <  1) THEN
     WRITE (XERN1, '(I8)') NEQ
     call XERMSG ('SLATEC', 'DLSOD', &
        'IN DDEBDF, THE NUMBER OF EQUATIONS MUST BE A ' // &
        'POSITIVE INTEGER.$$YOU HAVE CALLED THE CODE WITH NEQ = ' // &
        XERN1, 6, 1)
     IDID=-33
  end if
!
  NRTOLP = 0
  NATOLP = 0
  DO 60 K = 1, NEQ
     if (NRTOLP  <=  0) THEN
        if (RTOL(K)  <  0.) THEN
           WRITE (XERN1, '(I8)') K
           WRITE (XERN3, '(1PE15.6)') RTOL(K)
           call XERMSG ('SLATEC', 'DLSOD', &
              'IN DDEBDF, THE RELATIVE ERROR TOLERANCES MUST ' // &
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
        call XERMSG ('SLATEC', 'DLSOD', &
           'IN DDEBDF, THE ABSOLUTE ERROR ' // &
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
     if (SIGN(1.0D0,TOUT-T)  /=  SIGN(1.0D0,TSTOP-T) .OR. &
        ABS(TOUT-T)  >  ABS(TSTOP-T)) THEN
        WRITE (XERN3, '(1PE15.6)') TOUT
        WRITE (XERN4, '(1PE15.6)') TSTOP
        call XERMSG ('SLATEC', 'DLSOD', &
           'IN DDEBDF, YOU HAVE CALLED THE ' // &
           'CODE WITH TOUT = ' // XERN3 // '$$BUT YOU HAVE ' // &
           'ALSO TOLD THE CODE NOT TO INTEGRATE PAST THE POINT ' // &
           'TSTOP = ' // XERN4 // ' BY SETTING INFO(4) = 1.$$' // &
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
        call XERMSG ('SLATEC', 'DLSOD', &
           'IN DDEBDF, YOU HAVE CALLED THE CODE WITH T = TOUT = ' // &
           XERN3 // '$$THIS IS NOT ALLOWED ON CONTINUATION CALLS.', &
           9, 1)
        IDID=-33
     ENDIF
!
     if (T  /=  TOLD) THEN
        WRITE (XERN3, '(1PE15.6)') TOLD
        WRITE (XERN4, '(1PE15.6)') T
        call XERMSG ('SLATEC', 'DLSOD', &
           'IN DDEBDF, YOU HAVE CHANGED THE VALUE OF T FROM ' // &
           XERN3 // ' TO ' // XERN4 // &
           '  THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 10, 1)
        IDID=-33
     ENDIF
!
     if (INIT  /=  1) THEN
        if (DELSGN*(TOUT-T)  <  0.0D0) THEN
           WRITE (XERN3, '(1PE15.6)') TOUT
           call XERMSG ('SLATEC', 'DLSOD', &
              'IN DDEBDF, BY CALLING THE CODE WITH TOUT = ' // &
              XERN3 // ' YOU ARE ATTEMPTING TO CHANGE THE ' // &
              'DIRECTION OF INTEGRATION.$$THIS IS NOT ALLOWED ' // &
              'WITHOUT RESTARTING.', 11, 1)
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
        call XERMSG ('SLATEC', 'DLSOD', &
           'IN DDEBDF, INVALID INPUT WAS DETECTED ON ' // &
           'SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE TO PROCEED ' // &
           'BECAUSE YOU HAVE NOT CORRECTED THE PROBLEM, ' // &
           'SO EXECUTION IS BEING TERMINATED.', 12, 2)
     ENDIF
     return
  end if
!
!        ...............................................................
!
!             RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED
!             AS ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS
!             CASE, THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE
!             SMALLEST VALUE 100*U WHICH IS LIKELY TO BE REASONABLE FOR
!             THIS METHOD AND MACHINE
!
  DO 180 K = 1, NEQ
     if (RTOL(K) + ATOL(K)  >  0.0D0) go to 170
        RTOL(K) = 100.0D0*U
        IDID = -2
  170    CONTINUE
!     ...EXIT
     if (ITOL  ==  0) go to 190
  180 CONTINUE
  190 CONTINUE
!
  if (IDID  /=  (-2)) go to 200
!        RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A
!                                 SMALL POSITIVE VALUE
     IBEGIN = -1
  go to 460
  200 CONTINUE
!        BEGIN BLOCK PERMITTING ...EXITS TO 450
!           BEGIN BLOCK PERMITTING ...EXITS TO 430
!              BEGIN BLOCK PERMITTING ...EXITS TO 260
!                 BEGIN BLOCK PERMITTING ...EXITS TO 230
!
!                    BRANCH ON STATUS OF INITIALIZATION INDICATOR
!                           INIT=0 MEANS INITIAL DERIVATIVES AND
!                           NOMINAL STEP SIZE
!                                  AND DIRECTION NOT YET SET
!                           INIT=1 MEANS NOMINAL STEP SIZE AND
!                           DIRECTION NOT YET SET INIT=2 MEANS NO
!                           FURTHER INITIALIZATION REQUIRED
!
                 if (INIT  ==  0) go to 210
!                 ......EXIT
                    if (INIT  ==  1) go to 230
!              .........EXIT
                    go to 260
  210                CONTINUE
!
!                    ................................................
!
!                         MORE INITIALIZATION --
!                                             -- EVALUATE INITIAL
!                                             DERIVATIVES
!
                 INIT = 1
                 call DF(T,Y,YH(1,2),RPAR,IPAR)
                 NFE = 1
!                 ...EXIT
                 if (T  /=  TOUT) go to 230
                 IDID = 2
                 DO 220 L = 1, NEQ
                    YPOUT(L) = YH(L,2)
  220                CONTINUE
                 TOLD = T
!        ............EXIT
                 go to 450
  230             CONTINUE
!
!                 -- COMPUTE INITIAL STEP SIZE
!                 -- SAVE SIGN OF INTEGRATION DIRECTION
!                 -- SET INDEPENDENT AND DEPENDENT VARIABLES
!                                      X AND YH(*) FOR DSTOD
!
              LTOL = 1
              DO 240 L = 1, NEQ
                 if (ITOL  ==  1) LTOL = L
                 TOL = RTOL(LTOL)*ABS(Y(L)) + ATOL(LTOL)
                 if (TOL  ==  0.0D0) go to 390
                 EWT(L) = TOL
  240             CONTINUE
!
              BIG = SQRT(D1MACH(2))
              call DHSTRT(DF,NEQ,T,TOUT,Y,YH(1,2),EWT,1,U,BIG, &
                          YH(1,3),YH(1,4),YH(1,5),YH(1,6),RPAR, &
                          IPAR,H)
!
              DELSGN = SIGN(1.0D0,TOUT-T)
              X = T
              DO 250 L = 1, NEQ
                 YH(L,1) = Y(L)
                 YH(L,2) = H*YH(L,2)
  250             CONTINUE
              INIT = 2
  260          CONTINUE
!
!              ......................................................
!
!                 ON EACH call SET INFORMATION WHICH DETERMINES THE
!                 ALLOWED INTERVAL OF INTEGRATION BEFORE RETURNING
!                 WITH AN ANSWER AT TOUT
!
           DEL = TOUT - T
           ABSDEL = ABS(DEL)
!
!              ......................................................
!
!                 if ALREADY PAST OUTPUT POINT, INTERPOLATE AND
!                 return
!
  270          CONTINUE
!                 BEGIN BLOCK PERMITTING ...EXITS TO 400
!                    BEGIN BLOCK PERMITTING ...EXITS TO 380
                    if (ABS(X-T)  <  ABSDEL) go to 290
                       call DINTYD(TOUT,0,YH,NEQ,Y,INTFLG)
                       call DINTYD(TOUT,1,YH,NEQ,YPOUT,INTFLG)
                       IDID = 3
                       if (X  /=  TOUT) go to 280
                          IDID = 2
                          INTOUT = .FALSE.
  280                      CONTINUE
                       T = TOUT
                       TOLD = T
!        ..................EXIT
                       go to 450
  290                   CONTINUE
!
!                       if CANNOT GO PAST TSTOP AND SUFFICIENTLY
!                       CLOSE, EXTRAPOLATE AND RETURN
!
                    if (ITSTOP  /=  1) go to 310
                    if (ABS(TSTOP-X)  >=  100.0D0*U*ABS(X)) &
                       go to 310
                       DT = TOUT - X
                       DO 300 L = 1, NEQ
                          Y(L) = YH(L,1) + (DT/H)*YH(L,2)
  300                      CONTINUE
                       call DF(TOUT,Y,YPOUT,RPAR,IPAR)
                       NFE = NFE + 1
                       IDID = 3
                       T = TOUT
                       TOLD = T
!        ..................EXIT
                       go to 450
  310                   CONTINUE
!
                    if (IINTEG  ==  0 .OR. .NOT.INTOUT) go to 320
!
!                          INTERMEDIATE-OUTPUT MODE
!
                       IDID = 1
                    go to 370
  320                   CONTINUE
!
!                       .............................................
!
!                            MONITOR NUMBER OF STEPS ATTEMPTED
!
                    if (KSTEPS  <=  MAXNUM) go to 330
!
!                          A SIGNIFICANT AMOUNT OF WORK HAS BEEN
!                          EXPENDED
                       IDID = -1
                       KSTEPS = 0
                       IBEGIN = -1
                    go to 370
  330                   CONTINUE
!
!                          ..........................................
!
!                             LIMIT STEP SIZE AND SET WEIGHT VECTOR
!
                       HMIN = 100.0D0*U*ABS(X)
                       HA = MAX(ABS(H),HMIN)
                       if (ITSTOP  ==  1) &
                          HA = MIN(HA,ABS(TSTOP-X))
                       H = SIGN(HA,H)
                       LTOL = 1
                       DO 340 L = 1, NEQ
                          if (ITOL  ==  1) LTOL = L
                          EWT(L) = RTOL(LTOL)*ABS(YH(L,1)) &
                                   + ATOL(LTOL)
!                    .........EXIT
                          if (EWT(L)  <=  0.0D0) go to 380
  340                      CONTINUE
                       TOLFAC = U*DVNRMS(NEQ,YH,EWT)
!                 .........EXIT
                       if (TOLFAC  <=  1.0D0) go to 400
!
!                          TOLERANCES TOO SMALL
                       IDID = -2
                       TOLFAC = 2.0D0*TOLFAC
                       RTOL(1) = TOLFAC*RTOL(1)
                       ATOL(1) = TOLFAC*ATOL(1)
                       if (ITOL  ==  0) go to 360
                          DO 350 L = 2, NEQ
                             RTOL(L) = TOLFAC*RTOL(L)
                             ATOL(L) = TOLFAC*ATOL(L)
  350                         CONTINUE
  360                      CONTINUE
                       IBEGIN = -1
  370                   CONTINUE
!           ............EXIT
                    go to 430
  380                CONTINUE
!
!                    RELATIVE ERROR CRITERION INAPPROPRIATE
  390                CONTINUE
                 IDID = -3
                 IBEGIN = -1
!           .........EXIT
                 go to 430
  400             CONTINUE
!
!                 ...................................................
!
!                      TAKE A STEP
!
              call DSTOD(NEQ,Y,YH,NEQ,YH1,EWT,SAVF,ACOR,WM,IWM, &
                         DF,DJAC,RPAR,IPAR)
!
              JSTART = -2
              INTOUT = .TRUE.
           if (KFLAG  ==  0) go to 270
!
!              ......................................................
!
           if (KFLAG  ==  -1) go to 410
!
!                 REPEATED CORRECTOR CONVERGENCE FAILURES
              IDID = -6
              IBEGIN = -1
           go to 420
  410          CONTINUE
!
!                 REPEATED ERROR TEST FAILURES
              IDID = -7
              IBEGIN = -1
  420          CONTINUE
  430       CONTINUE
!
!           .........................................................
!
!                                  STORE VALUES BEFORE RETURNING TO
!                                  DDEBDF
        DO 440 L = 1, NEQ
           Y(L) = YH(L,1)
           YPOUT(L) = YH(L,2)/H
  440       CONTINUE
        T = X
        TOLD = T
        INTOUT = .FALSE.
  450    CONTINUE
  460 CONTINUE
  return
end
