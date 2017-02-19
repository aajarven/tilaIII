subroutine DRKFS (DF, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID, H, &
     TOLFAC, YP, F1, F2, F3, F4, F5, YS, TOLD, DTSIGN, U26, RER, &
     INIT, KSTEPS, KOP, IQUIT, STIFF, NONSTF, NTSTEP, NSTIFS, RPAR, &
     IPAR)
!
!! DRKFS integrates a system of ODE's for DDERKF.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DDERKF
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (DERKFS-S, DRKFS-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     Fehlberg Fourth-Fifth Order Runge-Kutta Method
! **********************************************************************
!
!     DRKFS integrates a system of first order ordinary differential
!     equations as described in the comments for DDERKF .
!
!     The arrays YP,F1,F2,F3,F4,F5,and YS  (of length at least NEQ)
!     appear in the call list for variable dimensioning purposes.
!
!     The variables H,TOLFAC,TOLD,DTSIGN,U26,RER,INIT,KSTEPS,KOP,IQUIT,
!     STIFF,NONSTF,NTSTEP, and NSTIFS are used internally by the code
!     and appear in the call list to eliminate local retention of
!     variables between calls. Accordingly, these variables and the
!     array YP should not be altered.
!     Items of possible interest are
!         H  - An appropriate step size to be used for the next step
!         TOLFAC - Factor of change in the tolerances
!         YP - Derivative of solution vector at T
!         KSTEPS - Counter on the number of steps attempted
!
! **********************************************************************
!
!***SEE ALSO  DDERKF
!***ROUTINES CALLED  D1MACH, DFEHL, DHSTRT, DHVNRM, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891024  Changed references from DVNORM to DHVNRM.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls, change GOTOs to
!           IF-THEN-ELSEs.  (RWC)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DRKFS
!
  INTEGER IDID, INFO, INIT, IPAR, IQUIT, K, KOP, KSTEPS, KTOL, &
        MXKOP, MXSTEP, NATOLP, NEQ, NRTOLP, NSTIFS, NTSTEP
  DOUBLE PRECISION A, ATOL, BIG, D1MACH, &
        DT, DTSIGN, DHVNRM, DY, EE, EEOET, ES, ESTIFF, &
        ESTTOL, ET, F1, F2, F3, F4, F5, H, HMIN, REMIN, RER, RPAR, &
        RTOL, S, T, TOL, TOLD, TOLFAC, TOUT, U, U26, UTE, Y, YAVG, &
        YP, YS
  LOGICAL HFAILD,OUTPUT,STIFF,NONSTF
  CHARACTER*8 XERN1
  CHARACTER*16 XERN3, XERN4
!
  DIMENSION Y(*),YP(*),F1(*),F2(*),F3(*),F4(*),F5(*), &
            YS(*),INFO(15),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
!
  EXTERNAL DF
!
!     ..................................................................
!
!       A FIFTH ORDER METHOD WILL GENERALLY NOT BE CAPABLE OF DELIVERING
!       ACCURACIES NEAR LIMITING PRECISION ON COMPUTERS WITH LONG
!       WORDLENGTHS. TO PROTECT AGAINST LIMITING PRECISION DIFFICULTIES
!       ARISING FROM UNREASONABLE ACCURACY REQUESTS, AN APPROPRIATE
!       TOLERANCE THRESHOLD REMIN IS ASSIGNED FOR THIS METHOD. THIS
!       VALUE SHOULD NOT BE CHANGED ACROSS DIFFERENT MACHINES.
!
  SAVE REMIN, MXSTEP, MXKOP
  DATA REMIN /1.0D-12/
!
!     ..................................................................
!
!       THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
!       NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MXSTEP, THE
!       COUNTER IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE
!       EXCESSIVE WORK.
!
  DATA MXSTEP /500/
!
!     ..................................................................
!
!       INEFFICIENCY CAUSED BY TOO FREQUENT OUTPUT IS MONITORED BY
!       COUNTING THE NUMBER OF STEP SIZES WHICH ARE SEVERELY SHORTENED
!       DUE SOLELY TO THE CHOICE OF OUTPUT POINTS. WHEN THE NUMBER OF
!       ABUSES EXCEED MXKOP, THE COUNTER IS RESET TO ZERO AND THE USER
!       IS INFORMED ABOUT POSSIBLE MISUSE OF THE CODE.
!
  DATA MXKOP /100/
!
!     ..................................................................
!
!***FIRST EXECUTABLE STATEMENT  DRKFS
  if (INFO(1)  ==  0) THEN
!
! ON THE FIRST call , PERFORM INITIALIZATION --
!        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
!        FUNCTION ROUTINE  D1MACH. THE USER MUST MAKE SURE THAT THE
!        VALUES SET IN D1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
!
     U = D1MACH(4)
!                       -- SET ASSOCIATED MACHINE DEPENDENT PARAMETERS
     U26 = 26.0D0*U
     RER = 2.0D0*U + REMIN
!                       -- SET TERMINATION FLAG
     IQUIT = 0
!                       -- SET INITIALIZATION INDICATOR
     INIT = 0
!                       -- SET COUNTER FOR IMPACT OF OUTPUT POINTS
     KOP = 0
!                       -- SET COUNTER FOR ATTEMPTED STEPS
     KSTEPS = 0
!                       -- SET INDICATORS FOR STIFFNESS DETECTION
     STIFF = .FALSE.
     NONSTF = .FALSE.
!                       -- SET STEP COUNTERS FOR STIFFNESS DETECTION
     NTSTEP = 0
     NSTIFS = 0
!                       -- RESET INFO(1) FOR SUBSEQUENT CALLS
     INFO(1) = 1
  end if
!
!.......................................................................
!
!        CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY
!
  if (INFO(1)  /=  0 .AND. INFO(1)  /=  1) THEN
     WRITE (XERN1, '(I8)') INFO(1)
     call XERMSG ('SLATEC', 'DRKFS', &
        'IN DDERKF, INFO(1) MUST BE SET TO 0 ' // &
        'FOR THE START OF A NEW PROBLEM, AND MUST BE SET TO 1 ' // &
        'FOLLOWING AN INTERRUPTED TASK.  YOU ARE ATTEMPTING TO ' // &
        'CONTINUE THE INTEGRATION ILLEGALLY BY CALLING THE CODE ' // &
        'WITH  INFO(1) = ' // XERN1, 3, 1)
     IDID = -33
  end if
!
  if (INFO(2)  /=  0 .AND. INFO(2)  /=  1) THEN
     WRITE (XERN1, '(I8)') INFO(2)
     call XERMSG ('SLATEC', 'DRKFS', &
        'IN DDERKF, INFO(2) MUST BE 0 OR 1 ' // &
        'INDICATING SCALAR AND VECTOR ERROR TOLERANCES, ' // &
        'RESPECTIVELY.  YOU HAVE CALLED THE CODE WITH INFO(2) = ' // &
        XERN1, 4, 1)
     IDID = -33
  end if
!
  if (INFO(3)  /=  0 .AND. INFO(3)  /=  1) THEN
     WRITE (XERN1, '(I8)') INFO(3)
     call XERMSG ('SLATEC', 'DRKFS', &
        'IN DDERKF, INFO(3) MUST BE 0 OR 1 ' // &
        'INDICATING THE INTERVAL OR INTERMEDIATE-OUTPUT MODE OF ' // &
        'INTEGRATION, RESPECTIVELY.  YOU HAVE CALLED THE CODE ' // &
        'WITH  INFO(3) = ' // XERN1, 5, 1)
     IDID = -33
  end if
!
  if (NEQ  <  1) THEN
     WRITE (XERN1, '(I8)') NEQ
     call XERMSG ('SLATEC', 'DRKFS', &
        'IN DDERKF, THE NUMBER OF EQUATIONS ' // &
        'NEQ MUST BE A POSITIVE INTEGER.  YOU HAVE CALLED THE ' // &
        'CODE WITH NEQ = ' // XERN1, 6, 1)
     IDID = -33
  end if
!
  NRTOLP = 0
  NATOLP = 0
  DO 10 K=1,NEQ
     if (NRTOLP  ==  0 .AND. RTOL(K)  <  0.D0) THEN
        WRITE (XERN1, '(I8)') K
        WRITE (XERN3, '(1PE15.6)') RTOL(K)
        call XERMSG ('SLATEC', 'DRKFS', &
           'IN DDERKF, THE RELATIVE ERROR ' // &
           'TOLERANCES RTOL MUST BE NON-NEGATIVE.  YOU HAVE ' // &
           'CALLED THE CODE WITH  RTOL(' // XERN1 // ') = ' // &
           XERN3 // '.  IN THE CASE OF VECTOR ERROR TOLERANCES, ' // &
           'NO FURTHER CHECKING OF RTOL COMPONENTS IS DONE.', 7, 1)
        IDID = -33
        NRTOLP = 1
     ENDIF
!
     if (NATOLP  ==  0 .AND. ATOL(K)  <  0.D0) THEN
        WRITE (XERN1, '(I8)') K
        WRITE (XERN3, '(1PE15.6)') ATOL(K)
        call XERMSG ('SLATEC', 'DRKFS', &
           'IN DDERKF, THE ABSOLUTE ERROR ' // &
           'TOLERANCES ATOL MUST BE NON-NEGATIVE.  YOU HAVE ' // &
           'CALLED THE CODE WITH  ATOL(' // XERN1 // ') = ' // &
           XERN3 // '.  IN THE CASE OF VECTOR ERROR TOLERANCES, ' // &
           'NO FURTHER CHECKING OF ATOL COMPONENTS IS DONE.', 8, 1)
        IDID = -33
        NATOLP = 1
     ENDIF
!
     if (INFO(2)  ==  0) go to 20
     if (NATOLP > 0 .AND. NRTOLP > 0) go to 20
   10 CONTINUE
!
!
!     CHECK SOME CONTINUATION POSSIBILITIES
!
   20 if (INIT  /=  0) THEN
     if (T  ==  TOUT) THEN
        WRITE (XERN3, '(1PE15.6)') T
        call XERMSG ('SLATEC', 'DRKFS', &
           'IN DDERKF, YOU HAVE CALLED THE ' // &
           'CODE WITH  T = TOUT = ' // XERN3 // '$$THIS IS NOT ' // &
           'ALLOWED ON CONTINUATION CALLS.', 9, 1)
        IDID=-33
     ENDIF
!
     if (T  /=  TOLD) THEN
        WRITE (XERN3, '(1PE15.6)') TOLD
        WRITE (XERN4, '(1PE15.6)') T
        call XERMSG ('SLATEC', 'DRKFS', &
           'IN DDERKF, YOU HAVE CHANGED THE ' // &
           'VALUE OF T FROM ' // XERN3 // ' TO ' // XERN4 // &
           '$$THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 10, 1)
        IDID=-33
     ENDIF
!
     if (INIT  /=  1) THEN
        if (DTSIGN*(TOUT-T)  <  0.D0) THEN
           WRITE (XERN3, '(1PE15.6)') TOUT
           call XERMSG ('SLATEC', 'DRKFS', &
              'IN DDERKF, BY CALLING THE CODE WITH TOUT = ' // &
              XERN3 // ' YOU ARE ATTEMPTING TO CHANGE THE ' // &
              'DIRECTION OF INTEGRATION.$$THIS IS NOT ALLOWED ' // &
              'WITHOUT RESTARTING.', 11, 1)
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
        GOTO 540
     ELSE
        call XERMSG ('SLATEC', 'DRKFS', &
           'IN DDERKF, INVALID INPUT WAS ' // &
           'DETECTED ON SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE ' // &
           'TO PROCEED BECAUSE YOU HAVE NOT CORRECTED THE ' // &
           'PROBLEM, SO EXECUTION IS BEING TERMINATED.', 12, 2)
        return
     ENDIF
  end if
!
!           ............................................................
!
!                RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND
!                INTERPRETED AS ASKING FOR THE MOST ACCURATE SOLUTION
!                POSSIBLE. IN THIS CASE, THE RELATIVE ERROR TOLERANCE
!                RTOL IS RESET TO THE SMALLEST VALUE RER WHICH IS LIKELY
!                TO BE REASONABLE FOR THIS METHOD AND MACHINE.
!
        DO 190 K = 1, NEQ
           if (RTOL(K) + ATOL(K)  >  0.0D0) go to 180
              RTOL(K) = RER
              IDID = -2
  180          CONTINUE
!           ...EXIT
           if (INFO(2)  ==  0) go to 200
  190       CONTINUE
  200       CONTINUE
!
        if (IDID  /=  (-2)) go to 210
!
!              RTOL=ATOL=0 ON INPUT, SO RTOL WAS CHANGED TO A
!                                       SMALL POSITIVE VALUE
           TOLFAC = 1.0D0
        go to 530
  210       CONTINUE
!
!                       BRANCH ON STATUS OF INITIALIZATION INDICATOR
!                              INIT=0 MEANS INITIAL DERIVATIVES AND
!                              STARTING STEP SIZE
!                                     NOT YET COMPUTED
!                              INIT=1 MEANS STARTING STEP SIZE NOT YET
!                              COMPUTED INIT=2 MEANS NO FURTHER
!                              INITIALIZATION REQUIRED
!
                    if (INIT  ==  0) go to 220
!                    ......EXIT
                       if (INIT  ==  1) go to 240
!                 .........EXIT
                       go to 260
  220                   CONTINUE
!
!                       ................................................
!
!                            MORE INITIALIZATION --
!                                                -- EVALUATE INITIAL
!                                                DERIVATIVES
!
                    INIT = 1
                    A = T
                    call DF(A,Y,YP,RPAR,IPAR)
                    if (T  /=  TOUT) go to 230
!
!                          INTERVAL MODE
                       IDID = 2
                       T = TOUT
                       TOLD = T
!     .....................EXIT
                       go to 560
  230                   CONTINUE
  240                CONTINUE
!
!                    -- SET SIGN OF INTEGRATION DIRECTION  AND
!                    -- ESTIMATE STARTING STEP SIZE
!
                 INIT = 2
                 DTSIGN = SIGN(1.0D0,TOUT-T)
                 U = D1MACH(4)
                 BIG = SQRT(D1MACH(2))
                 UTE = U**0.375D0
                 DY = UTE*DHVNRM(Y,NEQ)
                 if (DY  ==  0.0D0) DY = UTE
                 KTOL = 1
                 DO 250 K = 1, NEQ
                    if (INFO(2)  ==  1) KTOL = K
                    TOL = RTOL(KTOL)*ABS(Y(K)) + ATOL(KTOL)
                    if (TOL  ==  0.0D0) TOL = DY*RTOL(KTOL)
                    F1(K) = TOL
  250                CONTINUE
!
                 call DHSTRT(DF,NEQ,T,TOUT,Y,YP,F1,4,U,BIG,F2,F3,F4, &
                             F5,RPAR,IPAR,H)
  260             CONTINUE
!
!                 ......................................................
!
!                      SET STEP SIZE FOR INTEGRATION IN THE DIRECTION
!                      FROM T TO TOUT AND SET OUTPUT POINT INDICATOR
!
              DT = TOUT - T
              H = SIGN(H,DT)
              OUTPUT = .FALSE.
!
!                 TEST TO SEE if DDERKF IS BEING SEVERELY IMPACTED BY
!                 TOO MANY OUTPUT POINTS
!
              if (ABS(H)  >=  2.0D0*ABS(DT)) KOP = KOP + 1
              if (KOP  <=  MXKOP) go to 270
!
!                    UNNECESSARY FREQUENCY OF OUTPUT IS RESTRICTING
!                                              THE STEP SIZE CHOICE
                 IDID = -5
                 KOP = 0
              go to 510
  270             CONTINUE
!
                 if (ABS(DT)  >  U26*ABS(T)) go to 290
!
!                       if TOO CLOSE TO OUTPUT POINT,EXTRAPOLATE AND
!                       return
!
                    DO 280 K = 1, NEQ
                       Y(K) = Y(K) + DT*YP(K)
  280                   CONTINUE
                    A = TOUT
                    call DF(A,Y,YP,RPAR,IPAR)
                    KSTEPS = KSTEPS + 1
                 go to 500
  290                CONTINUE
!                       BEGIN BLOCK PERMITTING ...EXITS TO 490
!
!                          *********************************************
!                          *********************************************
!                               STEP BY STEP INTEGRATION
!
  300                      CONTINUE
!                             BEGIN BLOCK PERMITTING ...EXITS TO 480
                             HFAILD = .FALSE.
!
!                                TO PROTECT AGAINST IMPOSSIBLE ACCURACY
!                                REQUESTS, COMPUTE A TOLERANCE FACTOR
!                                BASED ON THE REQUESTED ERROR TOLERANCE
!                                AND A LEVEL OF ACCURACY ACHIEVABLE AT
!                                LIMITING PRECISION
!
                             TOLFAC = 0.0D0
                             KTOL = 1
                             DO 330 K = 1, NEQ
                                if (INFO(2)  ==  1) KTOL = K
                                ET = RTOL(KTOL)*ABS(Y(K)) &
                                     + ATOL(KTOL)
                                if (ET  >  0.0D0) go to 310
                                   TOLFAC = MAX(TOLFAC, &
                                                  RER/RTOL(KTOL))
                                go to 320
  310                               CONTINUE
                                   TOLFAC = MAX(TOLFAC, &
                                                  ABS(Y(K)) &
                                                  *(RER/ET))
  320                               CONTINUE
  330                            CONTINUE
                             if (TOLFAC  <=  1.0D0) go to 340
!
!                          REQUESTED ERROR UNATTAINABLE DUE TO LIMITED
!                                                  PRECISION AVAILABLE
                                TOLFAC = 2.0D0*TOLFAC
                                IDID = -2
!              .....................EXIT
                                go to 520
  340                            CONTINUE
!
!                                SET SMALLEST ALLOWABLE STEP SIZE
!
                             HMIN = U26*ABS(T)
!
!                                ADJUST STEP SIZE if NECESSARY TO HIT
!                                THE OUTPUT POINT -- LOOK AHEAD TWO
!                                STEPS TO AVOID DRASTIC CHANGES IN THE
!                                STEP SIZE AND THUS LESSEN THE IMPACT OF
!                                OUTPUT POINTS ON THE CODE.  STRETCH THE
!                                STEP SIZE BY, AT MOST, AN AMOUNT EQUAL
!                                TO THE SAFETY FACTOR OF 9/10.
!
                             DT = TOUT - T
                             if (ABS(DT)  >=  2.0D0*ABS(H)) &
                                go to 370
                                if (ABS(DT)  >  ABS(H)/0.9D0) &
                                   go to 350
!
!                                      THE NEXT STEP, if SUCCESSFUL,
!                                      WILL COMPLETE THE INTEGRATION TO
!                                      THE OUTPUT POINT
!
                                   OUTPUT = .TRUE.
                                   H = DT
                                go to 360
  350                               CONTINUE
!
                                   H = 0.5D0*DT
  360                               CONTINUE
  370                            CONTINUE
!
!
!                                ***************************************
!                                     CORE INTEGRATOR FOR TAKING A
!                                     SINGLE STEP
!                                ***************************************
!                                     TO AVOID PROBLEMS WITH ZERO
!                                     CROSSINGS, RELATIVE ERROR IS
!                                     MEASURED USING THE AVERAGE OF THE
!                                     MAGNITUDES OF THE SOLUTION AT THE
!                                     BEGINNING AND END OF A STEP.
!                                     THE ERROR ESTIMATE FORMULA HAS
!                                     BEEN GROUPED TO CONTROL LOSS OF
!                                     SIGNIFICANCE.
!                                     LOCAL ERROR ESTIMATES FOR A FIRST
!                                     ORDER METHOD USING THE SAME
!                                     STEP SIZE AS THE FEHLBERG METHOD
!                                     ARE CALCULATED AS PART OF THE
!                                     TEST FOR STIFFNESS.
!                                     TO DISTINGUISH THE VARIOUS
!                                     ARGUMENTS, H IS NOT PERMITTED
!                                     TO BECOME SMALLER THAN 26 UNITS OF
!                                     ROUNDOFF IN T.  PRACTICAL LIMITS
!                                     ON THE CHANGE IN THE STEP SIZE ARE
!                                     ENFORCED TO SMOOTH THE STEP SIZE
!                                     SELECTION PROCESS AND TO AVOID
!                                     EXCESSIVE CHATTERING ON PROBLEMS
!                                     HAVING DISCONTINUITIES.  TO
!                                     PREVENT UNNECESSARY FAILURES, THE
!                                     CODE USES 9/10 THE STEP SIZE
!                                     IT ESTIMATES WILL SUCCEED.
!                                     AFTER A STEP FAILURE, THE STEP
!                                     SIZE IS NOT ALLOWED TO INCREASE
!                                     FOR THE NEXT ATTEMPTED STEP. THIS
!                                     MAKES THE CODE MORE EFFICIENT ON
!                                     PROBLEMS HAVING DISCONTINUITIES
!                                     AND MORE EFFECTIVE IN GENERAL
!                                     SINCE LOCAL EXTRAPOLATION IS BEING
!                                     USED AND EXTRA CAUTION SEEMS
!                                     WARRANTED.
!                                .......................................
!
!                                     MONITOR NUMBER OF STEPS ATTEMPTED
!
  380                            CONTINUE
                                if (KSTEPS  <=  MXSTEP) go to 390
!
!                                      A SIGNIFICANT AMOUNT OF WORK HAS
!                                      BEEN EXPENDED
                                   IDID = -1
                                   KSTEPS = 0
!              ........................EXIT
                                   if (.NOT.STIFF) go to 520
!
!                                      PROBLEM APPEARS TO BE STIFF
                                   IDID = -4
                                   STIFF = .FALSE.
                                   NONSTF = .FALSE.
                                   NTSTEP = 0
                                   NSTIFS = 0
!              ........................EXIT
                                   go to 520
  390                               CONTINUE
!
!                                   ADVANCE AN APPROXIMATE SOLUTION OVER
!                                   ONE STEP OF LENGTH H
!
                                call DFEHL(DF,NEQ,T,Y,H,YP,F1,F2,F3, &
                                           F4,F5,YS,RPAR,IPAR)
                                KSTEPS = KSTEPS + 1
!
!                                   ....................................
!
!                                        COMPUTE AND TEST ALLOWABLE
!                                        TOLERANCES VERSUS LOCAL ERROR
!                                        ESTIMATES.  NOTE THAT RELATIVE
!                                        ERROR IS MEASURED WITH RESPECT
!                                        TO THE AVERAGE OF THE
!                                        MAGNITUDES OF THE SOLUTION AT
!                                        THE BEGINNING AND END OF THE
!                                        STEP.  LOCAL ERROR ESTIMATES
!                                        FOR A SPECIAL FIRST ORDER
!                                        METHOD ARE CALCULATED ONLY WHEN
!                                        THE STIFFNESS DETECTION IS
!                                        TURNED ON.
!
                                EEOET = 0.0D0
                                ESTIFF = 0.0D0
                                KTOL = 1
                                DO 420 K = 1, NEQ
                                   YAVG = 0.5D0 &
                                          *(ABS(Y(K)) &
                                            + ABS(YS(K)))
                                   if (INFO(2)  ==  1) KTOL = K
                                   ET = RTOL(KTOL)*YAVG + ATOL(KTOL)
                                   if (ET  >  0.0D0) go to 400
!
!           PURE RELATIVE ERROR INAPPROPRIATE WHEN SOLUTION
!                                                  VANISHES
                                      IDID = -3
!              ...........................EXIT
                                      go to 520
  400                                  CONTINUE
!
                                   EE = ABS((-2090.0D0*YP(K) &
                                              +(21970.0D0*F3(K) &
                                                -15048.0D0*F4(K))) &
                                             +(22528.0D0*F2(K) &
                                               -27360.0D0*F5(K)))
                                   if (STIFF .OR. NONSTF) go to 410
                                      ES = ABS(H &
                                                *(0.055455D0*YP(K) &
                                                  -0.035493D0*F1(K) &
                                                  -0.036571D0*F2(K) &
                                                  +0.023107D0*F3(K) &
                                                  -0.009515D0*F4(K) &
                                                  +0.003017D0*F5(K)) &
                                                  )
                                      ESTIFF = MAX(ESTIFF,ES/ET)
  410                                  CONTINUE
                                   EEOET = MAX(EEOET,EE/ET)
  420                               CONTINUE
!
                                ESTTOL = ABS(H)*EEOET/752400.0D0
!
!                                ...EXIT
                                if (ESTTOL  <=  1.0D0) go to 440
!
!                                   ....................................
!
!                                        UNSUCCESSFUL STEP
!
                                if (ABS(H)  >  HMIN) go to 430
!
!                             REQUESTED ERROR UNATTAINABLE AT SMALLEST
!                                                  ALLOWABLE STEP SIZE
                                   TOLFAC = 1.69D0*ESTTOL
                                   IDID = -2
!              ........................EXIT
                                   go to 520
  430                               CONTINUE
!
!                                   REDUCE THE STEP SIZE , TRY AGAIN
!                                   THE DECREASE IS LIMITED TO A FACTOR
!                                   OF 1/10
!
                                HFAILD = .TRUE.
                                OUTPUT = .FALSE.
                                S = 0.1D0
                                if (ESTTOL  <  59049.0D0) &
                                   S = 0.9D0/ESTTOL**0.2D0
                                H = SIGN(MAX(S*ABS(H),HMIN),H)
                             go to 380
  440                            CONTINUE
!
!                                .......................................
!
!                                SUCCESSFUL STEP
!                                                  STORE SOLUTION AT T+H
!                                                  AND EVALUATE
!                                                  DERIVATIVES THERE
!
                             T = T + H
                             DO 450 K = 1, NEQ
                                Y(K) = YS(K)
  450                            CONTINUE
                             A = T
                             call DF(A,Y,YP,RPAR,IPAR)
!
!                                CHOOSE NEXT STEP SIZE
!                                THE INCREASE IS LIMITED TO A FACTOR OF
!                                5 if STEP FAILURE HAS JUST OCCURRED,
!                                NEXT
!                                   STEP SIZE IS NOT ALLOWED TO INCREASE
!
                             S = 5.0D0
                             if (ESTTOL  >  1.889568D-4) &
                                S = 0.9D0/ESTTOL**0.2D0
                             if (HFAILD) S = MIN(S,1.0D0)
                             H = SIGN(MAX(S*ABS(H),HMIN),H)
!
!                                .......................................
!
!                                     CHECK FOR STIFFNESS (IF NOT
!                                     ALREADY DETECTED)
!
!                                     IN A SEQUENCE OF 50 SUCCESSFUL
!                                     STEPS BY THE FEHLBERG METHOD, 25
!                                     SUCCESSFUL STEPS BY THE FIRST
!                                     ORDER METHOD INDICATES STIFFNESS
!                                     AND TURNS THE TEST OFF. if 26
!                                     FAILURES BY THE FIRST ORDER METHOD
!                                     OCCUR, THE TEST IS TURNED OFF
!                                     UNTIL THIS SEQUENCE OF 50 STEPS BY
!                                     THE FEHLBERG METHOD IS COMPLETED.
!
!                             ...EXIT
                             if (STIFF) go to 480
                             NTSTEP = MOD(NTSTEP+1,50)
                             if (NTSTEP  ==  1) NONSTF = .FALSE.
!                             ...EXIT
                             if (NONSTF) go to 480
                             if (ESTIFF  >  1.0D0) go to 460
!
!                                   SUCCESSFUL STEP WITH FIRST ORDER
!                                   METHOD
                                NSTIFS = NSTIFS + 1
!                                   TURN TEST OFF AFTER 25 INDICATIONS
!                                   OF STIFFNESS
                                if (NSTIFS  ==  25) STIFF = .TRUE.
                             go to 470
  460                            CONTINUE
!
!                                UNSUCCESSFUL STEP WITH FIRST ORDER
!                                METHOD
                             if (NTSTEP - NSTIFS  <=  25) go to 470
!               TURN STIFFNESS DETECTION OFF FOR THIS BLOCK OF
!                                                  FIFTY STEPS
                                NONSTF = .TRUE.
!                                   RESET STIFF STEP COUNTER
                                NSTIFS = 0
  470                            CONTINUE
  480                         CONTINUE
!
!                             ******************************************
!                                  END OF CORE INTEGRATOR
!                             ******************************************
!
!
!                                  SHOULD WE TAKE ANOTHER STEP
!
!                       ......EXIT
                          if (OUTPUT) go to 490
                       if (INFO(3)  ==  0) go to 300
!
!                          *********************************************
!                          *********************************************
!
!                               INTEGRATION SUCCESSFULLY COMPLETED
!
!                                           ONE-STEP MODE
                       IDID = 1
                       TOLD = T
!     .....................EXIT
                       go to 560
  490                   CONTINUE
  500                CONTINUE
!
!                    INTERVAL MODE
                 IDID = 2
                 T = TOUT
                 TOLD = T
!     ...............EXIT
                 go to 560
  510             CONTINUE
  520          CONTINUE
  530       CONTINUE
  540    CONTINUE
!
!        INTEGRATION TASK INTERRUPTED
!
     INFO(1) = -1
     TOLD = T
!     ...EXIT
     if (IDID  /=  (-2)) go to 560
!
!        THE ERROR TOLERANCES ARE INCREASED TO VALUES
!                WHICH ARE APPROPRIATE FOR CONTINUING
     RTOL(1) = TOLFAC*RTOL(1)
     ATOL(1) = TOLFAC*ATOL(1)
!     ...EXIT
     if (INFO(2)  ==  0) go to 560
     DO 550 K = 2, NEQ
        RTOL(K) = TOLFAC*RTOL(K)
        ATOL(K) = TOLFAC*ATOL(K)
  550    CONTINUE
  560 CONTINUE
  return
end
