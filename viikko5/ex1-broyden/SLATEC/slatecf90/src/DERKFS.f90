subroutine DERKFS (F, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID, H, &
     TOLFAC, YP, F1, F2, F3, F4, F5, YS, TOLD, DTSIGN, U26, RER, &
     INIT, KSTEPS, KOP, IQUIT, STIFF, NONSTF, NTSTEP, NSTIFS, RPAR, &
     IPAR)
!
!! DERKFS is subsidiary to DERKF.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (DERKFS-S, DRKFS-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     Fehlberg Fourth-Fifth order Runge-Kutta Method
! **********************************************************************
!
!     DERKFS integrates a system of first order ordinary differential
!     equations as described in the comments for DERKF .
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
!***SEE ALSO  DERKF
!***ROUTINES CALLED  DEFEHL, HSTART, HVNRM, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800501  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891024  Changed references from VNORM to HVNRM.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls, replace GOTOs with
!           IF-THEN-ELSEs.  (RWC)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DERKFS
!
  LOGICAL HFAILD,OUTPUT,STIFF,NONSTF
  CHARACTER*8 XERN1
  CHARACTER*16 XERN3, XERN4
!
  DIMENSION Y(*),YP(*),F1(*),F2(*),F3(*),F4(*),F5(*), &
            YS(*),INFO(15),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
!
  EXTERNAL F
!
!.......................................................................
!
!  A FIFTH ORDER METHOD WILL GENERALLY NOT BE CAPABLE OF DELIVERING
!  ACCURACIES NEAR LIMITING PRECISION ON COMPUTERS WITH LONG
!  WORDLENGTHS. TO PROTECT AGAINST LIMITING PRECISION DIFFICULTIES
!  ARISING FROM UNREASONABLE ACCURACY REQUESTS, AN APPROPRIATE
!  TOLERANCE THRESHOLD REMIN IS ASSIGNED FOR THIS METHOD. THIS VALUE
!  SHOULD NOT BE CHANGED ACROSS DIFFERENT MACHINES.
!
  SAVE REMIN, MXSTEP, MXKOP
  DATA REMIN/1.E-12/
!
!.......................................................................
!
!  THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
!  NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MXSTEP, THE COUNTER
!  IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE EXCESSIVE
!  WORK.
!
  DATA MXSTEP/500/
!
!.......................................................................
!
!  INEFFICIENCY CAUSED BY TOO FREQUENT OUTPUT IS MONITORED BY COUNTING
!  THE NUMBER OF STEP SIZES WHICH ARE SEVERELY SHORTENED DUE SOLELY TO
!  THE CHOICE OF OUTPUT POINTS. WHEN THE NUMBER OF ABUSES EXCEED MXKOP,
!  THE COUNTER IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE
!  MISUSE OF THE CODE.
!
  DATA MXKOP/100/
!
!.......................................................................
!
!***FIRST EXECUTABLE STATEMENT  DERKFS
  if (INFO(1)  ==  0) THEN
!
! ON THE FIRST call , PERFORM INITIALIZATION --
!        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
!        FUNCTION ROUTINE  R1MACH. THE USER MUST MAKE SURE THAT THE
!        VALUES SET IN R1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
!
     U = R1MACH(4)
!                       -- SET ASSOCIATED MACHINE DEPENDENT PARAMETERS
     U26 = 26.*U
     RER = 2.*U+REMIN
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
     call XERMSG ('SLATEC', 'DERKFS', &
        'IN DERKF, INFO(1) MUST BE SET TO 0 ' // &
        'FOR THE START OF A NEW PROBLEM, AND MUST BE SET TO 1 ' // &
        'FOLLOWING AN INTERRUPTED TASK.  YOU ARE ATTEMPTING TO ' // &
        'CONTINUE THE INTEGRATION ILLEGALLY BY CALLING THE CODE ' // &
        'WITH  INFO(1) = ' // XERN1, 3, 1)
     IDID = -33
  end if
!
  if (INFO(2)  /=  0 .AND. INFO(2)  /=  1) THEN
     WRITE (XERN1, '(I8)') INFO(2)
     call XERMSG ('SLATEC', 'DERKFS', &
        'IN DERKF, INFO(2) MUST BE 0 OR 1 INDICATING SCALAR ' // &
        'AND VECTOR ERROR TOLERANCES, RESPECTIVELY.  YOU HAVE ' // &
        'CALLED THE CODE WITH INFO(2) = ' // XERN1, 4, 1)
     IDID = -33
  end if
!
  if (INFO(3)  /=  0 .AND. INFO(3)  /=  1) THEN
     WRITE (XERN1, '(I8)') INFO(3)
     call XERMSG ('SLATEC', 'DERKFS', &
        'IN DERKF, INFO(3) MUST BE 0 OR 1 INDICATING THE ' // &
        'OR INTERMEDIATE-OUTPUT MODE OF INTEGRATION, ' // &
        'RESPECTIVELY.  YOU HAVE CALLED THE CODE ' // &
        'WITH  INFO(3) = ' // XERN1, 5, 1)
     IDID = -33
  end if
!
  if (NEQ  <  1) THEN
     WRITE (XERN1, '(I8)') NEQ
     call XERMSG ('SLATEC', 'DERKFS', &
        'IN DERKF, THE NUMBER OF EQUATIONS NEQ MUST BE A ' // &
        'POSITIVE INTEGER.  YOU HAVE CALLED THE ' // &
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
        call XERMSG ('SLATEC', 'DERKFS', &
           'IN DERKF, THE RELATIVE ERROR ' // &
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
        call XERMSG ('SLATEC', 'DERKFS', &
           'IN DERKF, THE ABSOLUTE ERROR ' // &
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
        call XERMSG ('SLATEC', 'DERKFS', &
           'IN DERKF, YOU HAVE CALLED THE ' // &
           'CODE WITH  T = TOUT = ' // XERN3 // '$$THIS IS NOT ' // &
           'ALLOWED ON CONTINUATION CALLS.', 9, 1)
        IDID=-33
     ENDIF
!
     if (T  /=  TOLD) THEN
        WRITE (XERN3, '(1PE15.6)') TOLD
        WRITE (XERN4, '(1PE15.6)') T
        call XERMSG ('SLATEC', 'DERKFS', &
           'IN DERKF, YOU HAVE CHANGED THE ' // &
           'VALUE OF T FROM ' // XERN3 // ' TO ' // XERN4 // &
           '$$THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 10, 1)
        IDID=-33
     ENDIF
!
     if (INIT  /=  1) THEN
        if (DTSIGN*(TOUT-T)  <  0.D0) THEN
           WRITE (XERN3, '(1PE15.6)') TOUT
           call XERMSG ('SLATEC', 'DERKFS', &
              'IN DERKF, BY CALLING THE CODE ' // &
              'WITH TOUT = ' // XERN3 // ' YOU ARE ATTEMPTING ' // &
              'TO CHANGE THE DIRECTION OF INTEGRATION.$$THIS IS ' // &
              'NOT ALLOWED WITHOUT RESTARTING.', 11, 1)
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
        GOTO 909
     ELSE
        call XERMSG ('SLATEC', 'DERKFS', &
           'IN DERKF, INVALID INPUT WAS ' // &
           'DETECTED ON SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE ' // &
           'TO PROCEED BECAUSE YOU HAVE NOT CORRECTED THE ' // &
           'PROBLEM, SO EXECUTION IS BEING TERMINATED.', 12, 2)
        return
     ENDIF
  end if
!
!.......................................................................
!
!     RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED AS
!     ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS CASE,
!     THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE SMALLEST VALUE
!     RER WHICH IS LIKELY TO BE REASONABLE FOR THIS METHOD AND MACHINE.
!
  DO 50 K=1,NEQ
    if (RTOL(K)+ATOL(K)  >  0.) go to 45
    RTOL(K)=RER
    IDID=-2
   45   if (INFO(2)  ==  0) go to 55
   50   CONTINUE
!
   55 if (IDID  /=  (-2)) go to 60
!
!                       RTOL=ATOL=0 ON INPUT, SO RTOL WAS CHANGED TO A
!                                                SMALL POSITIVE VALUE
  TOLFAC=1.
  go to 909
!
!     BRANCH ON STATUS OF INITIALIZATION INDICATOR
!            INIT=0 MEANS INITIAL DERIVATIVES AND STARTING STEP SIZE
!                   NOT YET COMPUTED
!            INIT=1 MEANS STARTING STEP SIZE NOT YET COMPUTED
!            INIT=2 MEANS NO FURTHER INITIALIZATION REQUIRED
!
   60 if (INIT  ==  0) go to 65
  if (INIT  ==  1) go to 70
  go to 80
!
!.......................................................................
!
!     MORE INITIALIZATION --
!                         -- EVALUATE INITIAL DERIVATIVES
!
   65 INIT=1
  A=T
  call F(A,Y,YP,RPAR,IPAR)
  if (T  ==  TOUT) go to 666
!
!                         -- SET SIGN OF INTEGRATION DIRECTION  AND
!                         -- ESTIMATE STARTING STEP SIZE
!
   70 INIT=2
  DTSIGN=SIGN(1.,TOUT-T)
  U=R1MACH(4)
  BIG=SQRT(R1MACH(2))
  UTE=U**0.375
  DY=UTE*HVNRM(Y,NEQ)
  if (DY  ==  0.) DY=UTE
  KTOL=1
  DO 75 K=1,NEQ
    if (INFO(2)  ==  1)  KTOL=K
    TOL=RTOL(KTOL)*ABS(Y(K))+ATOL(KTOL)
    if (TOL  ==  0.) TOL=DY*RTOL(KTOL)
   75   F1(K)=TOL
!
  call HSTART (F,NEQ,T,TOUT,Y,YP,F1,4,U,BIG,F2,F3,F4,F5,RPAR,IPAR,H)
!
!.......................................................................
!
!     SET STEP SIZE FOR INTEGRATION IN THE DIRECTION FROM T TO TOUT
!     AND SET OUTPUT POINT INDICATOR
!
   80 DT=TOUT-T
  H=SIGN(H,DT)
  OUTPUT= .FALSE.
!
!     TEST TO SEE if DERKF IS BEING SEVERELY IMPACTED BY TOO MANY
!     OUTPUT POINTS
!
  if (ABS(H)  >=  2.*ABS(DT)) KOP=KOP+1
  if (KOP  <=  MXKOP) go to 85
!
!                       UNNECESSARY FREQUENCY OF OUTPUT IS RESTRICTING
!                                                 THE STEP SIZE CHOICE
  IDID=-5
  KOP=0
  go to 909
!
   85 if (ABS(DT)  >  U26*ABS(T)) go to 100
!
!     if TOO CLOSE TO OUTPUT POINT,EXTRAPOLATE AND RETURN
!
  DO 90 K=1,NEQ
   90   Y(K)=Y(K)+DT*YP(K)
  A=TOUT
  call F(A,Y,YP,RPAR,IPAR)
  KSTEPS=KSTEPS+1
  go to 666
!
! **********************************************************************
! **********************************************************************
!     STEP BY STEP INTEGRATION
!
  100 HFAILD= .FALSE.
!
!     TO PROTECT AGAINST IMPOSSIBLE ACCURACY REQUESTS, COMPUTE A
!     TOLERANCE FACTOR BASED ON THE REQUESTED ERROR TOLERANCE AND A
!     LEVEL OF ACCURACY ACHIEVABLE AT LIMITING PRECISION
!
  TOLFAC=0.
  KTOL=1
  DO 125 K=1,NEQ
    if (INFO(2)  ==  1) KTOL=K
    ET=RTOL(KTOL)*ABS(Y(K))+ATOL(KTOL)
    if (ET  >  0.) go to 120
    TOLFAC=MAX(TOLFAC,RER/RTOL(KTOL))
    go to 125
  120   TOLFAC=MAX(TOLFAC,ABS(Y(K))*(RER/ET))
  125   CONTINUE
  if (TOLFAC  <=  1.) go to 150
!
!                       REQUESTED ERROR UNATTAINABLE DUE TO LIMITED
!                                               PRECISION AVAILABLE
  TOLFAC=2.*TOLFAC
  IDID=-2
  go to 909
!
!     SET SMALLEST ALLOWABLE STEP SIZE
!
  150 HMIN=U26*ABS(T)
!
!     ADJUST STEP SIZE if NECESSARY TO HIT THE OUTPUT POINT --
!     LOOK AHEAD TWO STEPS TO AVOID DRASTIC CHANGES IN THE STEP SIZE AND
!     THUS LESSEN THE IMPACT OF OUTPUT POINTS ON THE CODE.
!     STRETCH THE STEP SIZE BY, AT MOST, AN AMOUNT EQUAL TO THE
!     SAFETY FACTOR OF 9/10.
!
  DT=TOUT-T
  if (ABS(DT)  >=  2.*ABS(H)) go to 200
  if (ABS(DT)  >  ABS(H)/0.9) go to 175
!
!     THE NEXT STEP, if SUCCESSFUL, WILL COMPLETE THE INTEGRATION TO
!     THE OUTPUT POINT
!
  OUTPUT= .TRUE.
  H=DT
  go to 200
!
  175 H=0.5*DT
!
!
! **********************************************************************
!     CORE INTEGRATOR FOR TAKING A SINGLE STEP
! **********************************************************************
!     TO AVOID PROBLEMS WITH ZERO CROSSINGS, RELATIVE ERROR IS MEASURED
!     USING THE AVERAGE OF THE MAGNITUDES OF THE SOLUTION AT THE
!     BEGINNING AND END OF A STEP.
!     THE ERROR ESTIMATE FORMULA HAS BEEN GROUPED TO CONTROL LOSS OF
!     SIGNIFICANCE.
!     LOCAL ERROR ESTIMATES FOR A FIRST ORDER METHOD USING THE SAME
!     STEP SIZE AS THE FEHLBERG METHOD ARE CALCULATED AS PART OF THE
!     TEST FOR STIFFNESS.
!     TO DISTINGUISH THE VARIOUS ARGUMENTS, H IS NOT PERMITTED
!     TO BECOME SMALLER THAN 26 UNITS OF ROUNDOFF IN T.
!     PRACTICAL LIMITS ON THE CHANGE IN THE STEP SIZE ARE ENFORCED TO
!     SMOOTH THE STEP SIZE SELECTION PROCESS AND TO AVOID EXCESSIVE
!     CHATTERING ON PROBLEMS HAVING DISCONTINUITIES.
!     TO PREVENT UNNECESSARY FAILURES, THE CODE USES 9/10 THE STEP SIZE
!     IT ESTIMATES WILL SUCCEED.
!     AFTER A STEP FAILURE, THE STEP SIZE IS NOT ALLOWED TO INCREASE FOR
!     THE NEXT ATTEMPTED STEP. THIS MAKES THE CODE MORE EFFICIENT ON
!     PROBLEMS HAVING DISCONTINUITIES AND MORE EFFECTIVE IN GENERAL
!     SINCE LOCAL EXTRAPOLATION IS BEING USED AND EXTRA CAUTION SEEMS
!     WARRANTED.
!.......................................................................
!
!     MONITOR NUMBER OF STEPS ATTEMPTED
!
  200 if (KSTEPS  <=  MXSTEP) go to 222
!
!                       A SIGNIFICANT AMOUNT OF WORK HAS BEEN EXPENDED
  IDID=-1
  KSTEPS=0
  if (.NOT. STIFF) go to 909
!
!                       PROBLEM APPEARS TO BE STIFF
  IDID=-4
  STIFF= .FALSE.
  NONSTF= .FALSE.
  NTSTEP=0
  NSTIFS=0
  go to 909
!
!     ADVANCE AN APPROXIMATE SOLUTION OVER ONE STEP OF LENGTH H
!
  222 call DEFEHL(F,NEQ,T,Y,H,YP,F1,F2,F3,F4,F5,YS,RPAR,IPAR)
  KSTEPS=KSTEPS+1
!
!.......................................................................
!
!     COMPUTE AND TEST ALLOWABLE TOLERANCES VERSUS LOCAL ERROR
!     ESTIMATES.  NOTE THAT RELATIVE ERROR IS MEASURED WITH RESPECT TO
!     THE AVERAGE OF THE MAGNITUDES OF THE SOLUTION AT THE BEGINNING
!     AND END OF THE STEP.
!     LOCAL ERROR ESTIMATES FOR A SPECIAL FIRST ORDER METHOD ARE
!     CALCULATED ONLY WHEN THE STIFFNESS DETECTION IS TURNED ON.
!
  EEOET=0.
  ESTIFF=0.
  KTOL=1
  DO 350 K=1,NEQ
    YAVG=0.5*(ABS(Y(K))+ABS(YS(K)))
    if (INFO(2)  ==  1) KTOL=K
    ET=RTOL(KTOL)*YAVG+ATOL(KTOL)
    if (ET  >  0.) go to 325
!
!                       PURE RELATIVE ERROR INAPPROPRIATE WHEN SOLUTION
!                                                              VANISHES
    IDID=-3
    go to 909
!
  325   EE=ABS((-2090.*YP(K)+(21970.*F3(K)-15048.*F4(K)))+ &
                          (22528.*F2(K)-27360.*F5(K)))
    if (STIFF .OR. NONSTF) go to 350
    ES=ABS(H*(0.055455*YP(K)-0.035493*F1(K)-0.036571*F2(K)+ &
              0.023107*F3(K)-0.009515*F4(K)+0.003017*F5(K)))
    ESTIFF=MAX(ESTIFF,ES/ET)
  350   EEOET=MAX(EEOET,EE/ET)
!
  ESTTOL=ABS(H)*EEOET/752400.
!
  if (ESTTOL  <=  1.) go to 500
!
!.......................................................................
!
!     UNSUCCESSFUL STEP
!
  if (ABS(H)  >  HMIN) go to 400
!
!                       REQUESTED ERROR UNATTAINABLE AT SMALLEST
!                                            ALLOWABLE STEP SIZE
  TOLFAC=1.69*ESTTOL
  IDID=-2
  go to 909
!
!                       REDUCE THE STEP SIZE , TRY AGAIN
!                       THE DECREASE IS LIMITED TO A FACTOR OF 1/10
!
  400 HFAILD= .TRUE.
  OUTPUT= .FALSE.
  S=0.1
  if (ESTTOL  <  59049.) S=0.9/ESTTOL**0.2
  H=SIGN(MAX(S*ABS(H),HMIN),H)
  go to 200
!
!.......................................................................
!
!     SUCCESSFUL STEP
!                       STORE SOLUTION AT T+H
!                       AND EVALUATE DERIVATIVES THERE
!
  500 T=T+H
  DO 525 K=1,NEQ
  525   Y(K)=YS(K)
  A=T
  call F(A,Y,YP,RPAR,IPAR)
!
!                       CHOOSE NEXT STEP SIZE
!                       THE INCREASE IS LIMITED TO A FACTOR OF 5
!                       if STEP FAILURE HAS JUST OCCURRED, NEXT
!                          STEP SIZE IS NOT ALLOWED TO INCREASE
!
  S=5.
  if (ESTTOL  >  1.889568E-4) S=0.9/ESTTOL**0.2
  if (HFAILD) S=MIN(S,1.)
  H=SIGN(MAX(S*ABS(H),HMIN),H)
!
!.......................................................................
!
!     CHECK FOR STIFFNESS (IF NOT ALREADY DETECTED)
!
!     IN A SEQUENCE OF 50 SUCCESSFUL STEPS BY THE FEHLBERG METHOD, 25
!     SUCCESSFUL STEPS BY THE FIRST ORDER METHOD INDICATES STIFFNESS
!     AND TURNS THE TEST OFF. if 26 FAILURES BY THE FIRST ORDER METHOD
!     OCCUR, THE TEST IS TURNED OFF UNTIL THIS SEQUENCE OF 50 STEPS
!     BY THE FEHLBERG METHOD IS COMPLETED.
!
  if (STIFF) go to 600
  NTSTEP=MOD(NTSTEP+1,50)
  if (NTSTEP  ==  1) NONSTF= .FALSE.
  if (NONSTF) go to 600
  if (ESTIFF  >  1.) go to 550
!
!                       SUCCESSFUL STEP WITH FIRST ORDER METHOD
  NSTIFS=NSTIFS+1
!                       TURN TEST OFF AFTER 25 INDICATIONS OF STIFFNESS
  if (NSTIFS  ==  25) STIFF= .TRUE.
  go to 600
!
!                       UNSUCCESSFUL STEP WITH FIRST ORDER METHOD
  550 if (NTSTEP-NSTIFS  <=  25) go to 600
!                       TURN STIFFNESS DETECTION OFF FOR THIS BLOCK OF
!                                                          FIFTY STEPS
  NONSTF= .TRUE.
!                       RESET STIFF STEP COUNTER
  NSTIFS=0
!
! **********************************************************************
!     END OF CORE INTEGRATOR
! **********************************************************************
!
!
!     SHOULD WE TAKE ANOTHER STEP
!
  600 if (OUTPUT) go to 666
  if (INFO(3)  ==  0) go to 100
!
! **********************************************************************
! **********************************************************************
!
!     INTEGRATION SUCCESSFULLY COMPLETED
!
!                 ONE-STEP MODE
  IDID=1
  TOLD=T
  return
!
!                 INTERVAL MODE
  666 IDID=2
  T=TOUT
  TOLD=T
  return
!
!     INTEGRATION TASK INTERRUPTED
!
  909 INFO(1)=-1
  TOLD=T
  if (IDID  /=  (-2)) RETURN
!
!                       THE ERROR TOLERANCES ARE INCREASED TO VALUES
!                               WHICH ARE APPROPRIATE FOR CONTINUING
  RTOL(1)=TOLFAC*RTOL(1)
  ATOL(1)=TOLFAC*ATOL(1)
  if (INFO(2)  ==  0) RETURN
  DO 939 K=2,NEQ
    RTOL(K)=TOLFAC*RTOL(K)
  939   ATOL(K)=TOLFAC*ATOL(K)
  return
end
