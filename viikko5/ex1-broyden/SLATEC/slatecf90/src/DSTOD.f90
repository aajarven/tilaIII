subroutine DSTOD (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR, WM, IWM, &
     DF, DJAC, RPAR, IPAR)
!
!! DSTOD is subsidiary to DDEBDF.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (STOD-S, DSTOD-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   DSTOD integrates a system of first order odes over one step in the
!   integrator package DDEBDF.
! ----------------------------------------------------------------------
! DSTOD  performs one step of the integration of an initial value
! problem for a system of ordinary differential equations.
! Note.. DSTOD  is independent of the value of the iteration method
! indicator MITER, when this is  /=  0, and hence is independent
! of the type of chord method used, or the Jacobian structure.
! Communication with DSTOD  is done with the following variables..
!
! Y      = An array of length  >=  N used as the Y argument in
!          all calls to DF and DJAC.
! NEQ    = Integer array containing problem size in NEQ(1), and
!          passed as the NEQ argument in all calls to DF and DJAC.
! YH     = An NYH by LMAX array containing the dependent variables
!          and their approximate scaled derivatives, where
!          LMAX = MAXORD + 1.  YH(I,J+1) contains the approximate
!          J-th derivative of Y(I), scaled by H**J/FACTORIAL(J)
!          (J = 0,1,...,NQ).  On entry for the first step, the first
!          two columns of YH must be set from the initial values.
! NYH    = A constant integer  >=  N, the first dimension of YH.
! YH1    = A one-dimensional array occupying the same space as YH.
! EWT    = An array of N elements with which the estimated local
!          errors in YH are compared.
! SAVF   = An array of working storage, of length N.
! ACOR   = A work array of length N, used for the accumulated
!          corrections.  On a successful return, ACOR(I) contains
!          the estimated one-step local error in Y(I).
! WM,IWM = DOUBLE PRECISION and INTEGER work arrays associated with
!          matrix operations in chord iteration (MITER  /=  0).
! DPJAC   = Name of routine to evaluate and preprocess Jacobian matrix
!          if a chord method is being used.
! DSLVS   = Name of routine to solve linear system in chord iteration.
! H      = The step size to be attempted on the next step.
!          H is altered by the error control algorithm during the
!          problem.  H can be either positive or negative, but its
!          sign must remain constant throughout the problem.
! HMIN   = The minimum absolute value of the step size H to be used.
! HMXI   = Inverse of the maximum absolute value of H to be used.
!          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
!          HMIN and HMXI may be changed at any time, but will not
!          take effect until the next change of H is considered.
! TN     = The independent variable. TN is updated on each step taken.
! JSTART = An integer used for input only, with the following
!          values and meanings..
!               0  Perform the first step.
!            > 0  Take a new step continuing from the last.
!              -1  Take the next step with a new value of H, MAXORD,
!                    N, METH, MITER, and/or matrix parameters.
!              -2  Take the next step with a new value of H,
!                    but with other inputs unchanged.
!          On return, JSTART is set to 1 to facilitate continuation.
! KFLAG  = a completion code with the following meanings..
!               0  The step was successful.
!              -1  The requested error could not be achieved.
!              -2  Corrector convergence could not be achieved.
!          A return with KFLAG = -1 or -2 means either
!          ABS(H) = HMIN or 10 consecutive failures occurred.
!          On a return with KFLAG negative, the values of TN and
!          the YH array are as of the beginning of the last
!          step, and H is the last step size attempted.
! MAXORD = The maximum order of integration method to be allowed.
! METH/MITER = The method flags.  See description in driver.
! N      = The number of first-order differential equations.
! ----------------------------------------------------------------------
!
!***SEE ALSO  DDEBDF
!***ROUTINES CALLED  DCFOD, DPJAC, DSLVS, DVNRMS
!***COMMON BLOCKS    DDEBD1
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!   920422  Changed DIMENSION statement.  (WRB)
!***END PROLOGUE  DSTOD
!
  INTEGER I, I1, IALTH, IER, IOD, IOWND, IPAR, IPUP, IREDO, IRET, &
        IWM, J, JB, JSTART, KFLAG, KSTEPS, L, LMAX, M, MAXORD, &
        MEO, METH, MITER, N, NCF, NEQ, NEWQ, NFE, NJE, NQ, NQNYH, &
        NQU, NST, NSTEPJ, NYH
  DOUBLE PRECISION ACOR, CONIT, CRATE, DCON, DDN, &
        DEL, DELP, DSM, DUP, DVNRMS, EL, EL0, ELCO, &
        EWT, EXDN, EXSM, EXUP, H, HMIN, HMXI, HOLD, HU, R, RC, &
        RH, RHDN, RHSM, RHUP, RMAX, ROWND, RPAR, SAVF, TESCO, &
        TN, TOLD, UROUND, WM, Y, YH, YH1
  EXTERNAL DF, DJAC
!
  DIMENSION Y(*),YH(NYH,*),YH1(*),EWT(*),SAVF(*),ACOR(*),WM(*), &
            IWM(*),RPAR(*),IPAR(*)
  COMMON /DDEBD1/ ROWND,CONIT,CRATE,EL(13),ELCO(13,12),HOLD,RC,RMAX, &
                  TESCO(3,12),EL0,H,HMIN,HMXI,HU,TN,UROUND,IOWND(7), &
                  KSTEPS,IOD(6),IALTH,IPUP,LMAX,MEO,NQNYH,NSTEPJ, &
                  IER,JSTART,KFLAG,L,METH,MITER,MAXORD,N,NQ,NST,NFE, &
                  NJE,NQU
!
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 690
!        BEGIN BLOCK PERMITTING ...EXITS TO 60
!***FIRST EXECUTABLE STATEMENT  DSTOD
        KFLAG = 0
        TOLD = TN
        NCF = 0
        if (JSTART  >  0) go to 160
        if (JSTART  ==  -1) go to 10
           if (JSTART  ==  -2) go to 90
!              ---------------------------------------------------------
!               ON THE FIRST CALL, THE ORDER IS SET TO 1, AND OTHER
!               VARIABLES ARE INITIALIZED.  RMAX IS THE MAXIMUM RATIO BY
!               WHICH H CAN BE INCREASED IN A SINGLE STEP.  IT IS
!               INITIALLY 1.E4 TO COMPENSATE FOR THE SMALL INITIAL H,
!               BUT THEN IS NORMALLY EQUAL TO 10.  if A FAILURE OCCURS
!               (IN CORRECTOR CONVERGENCE OR ERROR TEST), RMAX IS SET AT
!               2 FOR THE NEXT INCREASE.
!              ---------------------------------------------------------
           LMAX = MAXORD + 1
           NQ = 1
           L = 2
           IALTH = 2
           RMAX = 10000.0D0
           RC = 0.0D0
           EL0 = 1.0D0
           CRATE = 0.7D0
           DELP = 0.0D0
           HOLD = H
           MEO = METH
           NSTEPJ = 0
           IRET = 3
        go to 50
   10       CONTINUE
!              BEGIN BLOCK PERMITTING ...EXITS TO 30
!                 ------------------------------------------------------
!                  THE FOLLOWING BLOCK HANDLES PRELIMINARIES NEEDED WHEN
!                  JSTART = -1.  IPUP IS SET TO MITER TO FORCE A MATRIX
!                  UPDATE.  if AN ORDER INCREASE IS ABOUT TO BE
!                  CONSIDERED (IALTH = 1), IALTH IS RESET TO 2 TO
!                  POSTPONE CONSIDERATION ONE MORE STEP.  if THE CALLER
!                  HAS CHANGED METH, DCFOD  IS CALLED TO RESET THE
!                  COEFFICIENTS OF THE METHOD.  if THE CALLER HAS
!                  CHANGED MAXORD TO A VALUE LESS THAN THE CURRENT
!                  ORDER NQ, NQ IS REDUCED TO MAXORD, AND A NEW H CHOSEN
!                  ACCORDINGLY.  if H IS TO BE CHANGED, YH MUST BE
!                  RESCALED.  if H OR METH IS BEING CHANGED, IALTH IS
!                  RESET TO L = NQ + 1 TO PREVENT FURTHER CHANGES IN H
!                  FOR THAT MANY STEPS.
!                 ------------------------------------------------------
              IPUP = MITER
              LMAX = MAXORD + 1
              if (IALTH  ==  1) IALTH = 2
              if (METH  ==  MEO) go to 20
                 call DCFOD(METH,ELCO,TESCO)
                 MEO = METH
!              ......EXIT
                 if (NQ  >  MAXORD) go to 30
                 IALTH = L
                 IRET = 1
!        ............EXIT
                 go to 60
   20             CONTINUE
              if (NQ  <=  MAXORD) go to 90
   30          CONTINUE
           NQ = MAXORD
           L = LMAX
           DO 40 I = 1, L
              EL(I) = ELCO(I,NQ)
   40          CONTINUE
           NQNYH = NQ*NYH
           RC = RC*EL(1)/EL0
           EL0 = EL(1)
           CONIT = 0.5D0/(NQ+2)
           DDN = DVNRMS(N,SAVF,EWT)/TESCO(1,L)
           EXDN = 1.0D0/L
           RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
           RH = MIN(RHDN,1.0D0)
           IREDO = 3
           if (H  ==  HOLD) go to 660
           RH = MIN(RH,ABS(H/HOLD))
           H = HOLD
           go to 100
   50       CONTINUE
!           ------------------------------------------------------------
!            DCFOD  IS CALLED TO GET ALL THE INTEGRATION COEFFICIENTS
!            FOR THE CURRENT METH.  THEN THE EL VECTOR AND RELATED
!            CONSTANTS ARE RESET WHENEVER THE ORDER NQ IS CHANGED, OR AT
!            THE START OF THE PROBLEM.
!           ------------------------------------------------------------
        call DCFOD(METH,ELCO,TESCO)
   60    CONTINUE
   70    CONTINUE
!           BEGIN BLOCK PERMITTING ...EXITS TO 680
           DO 80 I = 1, L
              EL(I) = ELCO(I,NQ)
   80          CONTINUE
           NQNYH = NQ*NYH
           RC = RC*EL(1)/EL0
           EL0 = EL(1)
           CONIT = 0.5D0/(NQ+2)
           go to (90,660,160), IRET
!              ---------------------------------------------------------
!               if H IS BEING CHANGED, THE H RATIO RH IS CHECKED AGAINST
!               RMAX, HMIN, AND HMXI, AND THE YH ARRAY RESCALED.  IALTH
!               IS SET TO L = NQ + 1 TO PREVENT A CHANGE OF H FOR THAT
!               MANY STEPS, UNLESS FORCED BY A CONVERGENCE OR ERROR TEST
!               FAILURE.
!              ---------------------------------------------------------
   90          CONTINUE
           if (H  ==  HOLD) go to 160
           RH = H/HOLD
           H = HOLD
           IREDO = 3
  100          CONTINUE
  110          CONTINUE
              RH = MIN(RH,RMAX)
              RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
              R = 1.0D0
              DO 130 J = 2, L
                 R = R*RH
                 DO 120 I = 1, N
                    YH(I,J) = YH(I,J)*R
  120                CONTINUE
  130             CONTINUE
              H = H*RH
              RC = RC*RH
              IALTH = L
              if (IREDO  /=  0) go to 150
                 RMAX = 10.0D0
                 R = 1.0D0/TESCO(2,NQU)
                 DO 140 I = 1, N
                    ACOR(I) = ACOR(I)*R
  140                CONTINUE
!     ...............EXIT
                 go to 690
  150             CONTINUE
!                 ------------------------------------------------------
!                  THIS SECTION COMPUTES THE PREDICTED VALUES BY
!                  EFFECTIVELY MULTIPLYING THE YH ARRAY BY THE PASCAL
!                  TRIANGLE MATRIX.  RC IS THE RATIO OF NEW TO OLD
!                  VALUES OF THE COEFFICIENT  H*EL(1).  WHEN RC DIFFERS
!                  FROM 1 BY MORE THAN 30 PERCENT, IPUP IS SET TO MITER
!                  TO FORCE DPJAC TO BE CALLED, if A JACOBIAN IS
!                  INVOLVED.  IN ANY CASE, DPJAC IS CALLED AT LEAST
!                  EVERY 20-TH STEP.
!                 ------------------------------------------------------
  160             CONTINUE
  170             CONTINUE
!                    BEGIN BLOCK PERMITTING ...EXITS TO 610
!                       BEGIN BLOCK PERMITTING ...EXITS TO 490
                       if (ABS(RC-1.0D0)  >  0.3D0) IPUP = MITER
                       if (NST  >=  NSTEPJ + 20) IPUP = MITER
                       TN = TN + H
                       I1 = NQNYH + 1
                       DO 190 JB = 1, NQ
                          I1 = I1 - NYH
                          DO 180 I = I1, NQNYH
                             YH1(I) = YH1(I) + YH1(I+NYH)
  180                         CONTINUE
  190                      CONTINUE
                       KSTEPS = KSTEPS + 1
!                          ---------------------------------------------
!                           UP TO 3 CORRECTOR ITERATIONS ARE TAKEN.  A
!                           CONVERGENCE TEST IS MADE ON THE R.M.S. NORM
!                           OF EACH CORRECTION, WEIGHTED BY THE ERROR
!                           WEIGHT VECTOR EWT.  THE SUM OF THE
!                           CORRECTIONS IS ACCUMULATED IN THE VECTOR
!                           ACOR(I).  THE YH ARRAY IS NOT ALTERED IN THE
!                           CORRECTOR LOOP.
!                          ---------------------------------------------
  200                      CONTINUE
                          M = 0
                          DO 210 I = 1, N
                             Y(I) = YH(I,1)
  210                         CONTINUE
                          call DF(TN,Y,SAVF,RPAR,IPAR)
                          NFE = NFE + 1
                          if (IPUP  <=  0) go to 220
!                                ---------------------------------------
!                                 if INDICATED, THE MATRIX P = I -
!                                 H*EL(1)*J IS REEVALUATED AND
!                                 PREPROCESSED BEFORE STARTING THE
!                                 CORRECTOR ITERATION.  IPUP IS SET TO 0
!                                 AS AN INDICATOR THAT THIS HAS BEEN
!                                 DONE.
!                                ---------------------------------------
                             IPUP = 0
                             RC = 1.0D0
                             NSTEPJ = NST
                             CRATE = 0.7D0
                             call DPJAC(NEQ,Y,YH,NYH,EWT,ACOR,SAVF, &
                                        WM,IWM,DF,DJAC,RPAR,IPAR)
!                          ......EXIT
                             if (IER  /=  0) go to 440
  220                         CONTINUE
                          DO 230 I = 1, N
                             ACOR(I) = 0.0D0
  230                         CONTINUE
  240                         CONTINUE
                             if (MITER  /=  0) go to 270
!                                   ------------------------------------
!                                    IN THE CASE OF FUNCTIONAL
!                                    ITERATION, UPDATE Y DIRECTLY FROM
!                                    THE RESULT OF THE LAST FUNCTION
!                                    EVALUATION.
!                                   ------------------------------------
                                DO 250 I = 1, N
                                   SAVF(I) = H*SAVF(I) - YH(I,2)
                                   Y(I) = SAVF(I) - ACOR(I)
  250                               CONTINUE
                                DEL = DVNRMS(N,Y,EWT)
                                DO 260 I = 1, N
                                   Y(I) = YH(I,1) + EL(1)*SAVF(I)
                                   ACOR(I) = SAVF(I)
  260                               CONTINUE
                             go to 300
  270                            CONTINUE
!                                   ------------------------------------
!                                    IN THE CASE OF THE CHORD METHOD,
!                                    COMPUTE THE CORRECTOR ERROR, AND
!                                    SOLVE THE LINEAR SYSTEM WITH THAT
!                                    AS RIGHT-HAND SIDE AND P AS
!                                    COEFFICIENT MATRIX.
!                                   ------------------------------------
                                DO 280 I = 1, N
                                   Y(I) = H*SAVF(I) &
                                          - (YH(I,2) + ACOR(I))
  280                               CONTINUE
                                call DSLVS(WM,IWM,Y,SAVF)
!                             ......EXIT
                                if (IER  /=  0) go to 430
                                DEL = DVNRMS(N,Y,EWT)
                                DO 290 I = 1, N
                                   ACOR(I) = ACOR(I) + Y(I)
                                   Y(I) = YH(I,1) + EL(1)*ACOR(I)
  290                               CONTINUE
  300                            CONTINUE
!                                ---------------------------------------
!                                 TEST FOR CONVERGENCE.  if M > 0, AN
!                                 ESTIMATE OF THE CONVERGENCE RATE
!                                 CONSTANT IS STORED IN CRATE, AND THIS
!                                 IS USED IN THE TEST.
!                                ---------------------------------------
                             if (M  /=  0) &
                                CRATE = MAX(0.2D0*CRATE,DEL/DELP)
                             DCON = DEL*MIN(1.0D0,1.5D0*CRATE) &
                                    /(TESCO(2,NQ)*CONIT)
                             if (DCON  >  1.0D0) go to 420
!                                   ------------------------------------
!                                    THE CORRECTOR HAS CONVERGED.  IPUP
!                                    IS SET TO -1 if MITER  /=  0, TO
!                                    SIGNAL THAT THE JACOBIAN INVOLVED
!                                    MAY NEED UPDATING LATER.  THE LOCAL
!                                    ERROR TEST IS MADE AND CONTROL
!                                    PASSES TO STATEMENT 500 if IT
!                                    FAILS.
!                                   ------------------------------------
                                if (MITER  /=  0) IPUP = -1
                                if (M  ==  0) DSM = DEL/TESCO(2,NQ)
                                if (M  >  0) &
                                   DSM = DVNRMS(N,ACOR,EWT) &
                                         /TESCO(2,NQ)
                                if (DSM  >  1.0D0) go to 380
!                                      BEGIN BLOCK
!                                      PERMITTING ...EXITS TO 360
!                                         ------------------------------
!                                          AFTER A SUCCESSFUL STEP,
!                                          UPDATE THE YH ARRAY.
!                                          CONSIDER CHANGING H if IALTH
!                                          = 1.  OTHERWISE DECREASE
!                                          IALTH BY 1.  if IALTH IS THEN
!                                          1 AND NQ  <  MAXORD, THEN
!                                          ACOR IS SAVED FOR USE IN A
!                                          POSSIBLE ORDER INCREASE ON
!                                          THE NEXT STEP.  if A CHANGE
!                                          IN H IS CONSIDERED, AN
!                                          INCREASE OR DECREASE IN ORDER
!                                          BY ONE IS CONSIDERED ALSO.  A
!                                          CHANGE IN H IS MADE ONLY IF
!                                          IT IS BY A FACTOR OF AT LEAST
!                                          1.1.  if NOT, IALTH IS SET TO
!                                          3 TO PREVENT TESTING FOR THAT
!                                          MANY STEPS.
!                                         ------------------------------
                                      KFLAG = 0
                                      IREDO = 0
                                      NST = NST + 1
                                      HU = H
                                      NQU = NQ
                                      DO 320 J = 1, L
                                         DO 310 I = 1, N
                                            YH(I,J) = YH(I,J) &
                                                      + EL(J) &
                                                        *ACOR(I)
  310                                        CONTINUE
  320                                     CONTINUE
                                      IALTH = IALTH - 1
                                      if (IALTH  /=  0) go to 340
!                                            ---------------------------
!                                             REGARDLESS OF THE SUCCESS
!                                             OR FAILURE OF THE STEP,
!                                             FACTORS RHDN, RHSM, AND
!                                             RHUP ARE COMPUTED, BY
!                                             WHICH H COULD BE
!                                             MULTIPLIED AT ORDER NQ -
!                                             1, ORDER NQ, OR ORDER NQ +
!                                             1, RESPECTIVELY.  IN THE
!                                             CASE OF FAILURE, RHUP =
!                                             0.0 TO AVOID AN ORDER
!                                             INCREASE.  THE LARGEST OF
!                                             THESE IS DETERMINED AND
!                                             THE NEW ORDER CHOSEN
!                                             ACCORDINGLY.  if THE ORDER
!                                             IS TO BE INCREASED, WE
!                                             COMPUTE ONE ADDITIONAL
!                                             SCALED DERIVATIVE.
!                                            ---------------------------
                                         RHUP = 0.0D0
!                       .....................EXIT
                                         if (L  ==  LMAX) go to 490
                                         DO 330 I = 1, N
                                            SAVF(I) = ACOR(I) &
                                                      - YH(I,LMAX)
  330                                        CONTINUE
                                         DUP = DVNRMS(N,SAVF,EWT) &
                                               /TESCO(3,NQ)
                                         EXUP = 1.0D0/(L+1)
                                         RHUP = 1.0D0 &
                                                /(1.4D0*DUP**EXUP &
                                                  + 0.0000014D0)
!                       .....................EXIT
                                         go to 490
  340                                     CONTINUE
!                                      ...EXIT
                                      if (IALTH  >  1) go to 360
!                                      ...EXIT
                                      if (L  ==  LMAX) go to 360
                                      DO 350 I = 1, N
                                         YH(I,LMAX) = ACOR(I)
  350                                     CONTINUE
  360                                  CONTINUE
                                   R = 1.0D0/TESCO(2,NQU)
                                   DO 370 I = 1, N
                                      ACOR(I) = ACOR(I)*R
  370                                  CONTINUE
!     .................................EXIT
                                   go to 690
  380                               CONTINUE
!                                   ------------------------------------
!                                    THE ERROR TEST FAILED.  KFLAG KEEPS
!                                    TRACK OF MULTIPLE FAILURES.
!                                    RESTORE TN AND THE YH ARRAY TO
!                                    THEIR PREVIOUS VALUES, AND PREPARE
!                                    TO TRY THE STEP AGAIN.  COMPUTE THE
!                                    OPTIMUM STEP SIZE FOR THIS OR ONE
!                                    LOWER ORDER.  AFTER 2 OR MORE
!                                    FAILURES, H IS FORCED TO DECREASE
!                                    BY A FACTOR OF 0.2 OR LESS.
!                                   ------------------------------------
                                KFLAG = KFLAG - 1
                                TN = TOLD
                                I1 = NQNYH + 1
                                DO 400 JB = 1, NQ
                                   I1 = I1 - NYH
                                   DO 390 I = I1, NQNYH
                                      YH1(I) = YH1(I) - YH1(I+NYH)
  390                                  CONTINUE
  400                               CONTINUE
                                RMAX = 2.0D0
                                if (ABS(H)  >  HMIN*1.00001D0) &
                                   go to 410
!                                      ---------------------------------
!                                       ALL RETURNS ARE MADE THROUGH
!                                       THIS SECTION.  H IS SAVED IN
!                                       HOLD TO ALLOW THE CALLER TO
!                                       CHANGE H ON THE NEXT STEP.
!                                      ---------------------------------
                                   KFLAG = -1
!     .................................EXIT
                                   go to 690
  410                               CONTINUE
!                    ...............EXIT
                                if (KFLAG  <=  -3) go to 610
                                IREDO = 2
                                RHUP = 0.0D0
!                       ............EXIT
                                go to 490
  420                            CONTINUE
                             M = M + 1
!                             ...EXIT
                             if (M  ==  3) go to 430
!                             ...EXIT
                             if (M  >=  2 .AND. DEL  >  2.0D0*DELP) &
                                go to 430
                             DELP = DEL
                             call DF(TN,Y,SAVF,RPAR,IPAR)
                             NFE = NFE + 1
                          go to 240
  430                         CONTINUE
!                             ------------------------------------------
!                              THE CORRECTOR ITERATION FAILED TO
!                              CONVERGE IN 3 TRIES.  if MITER  /=  0 AND
!                              THE JACOBIAN IS OUT OF DATE, DPJAC IS
!                              CALLED FOR THE NEXT TRY.  OTHERWISE THE
!                              YH ARRAY IS RETRACTED TO ITS VALUES
!                              BEFORE PREDICTION, AND H IS REDUCED, IF
!                              POSSIBLE.  if H CANNOT BE REDUCED OR 10
!                              FAILURES HAVE OCCURRED, EXIT WITH KFLAG =
!                              -2.
!                             ------------------------------------------
!                          ...EXIT
                          if (IPUP  ==  0) go to 440
                          IPUP = MITER
                       go to 200
  440                      CONTINUE
                       TN = TOLD
                       NCF = NCF + 1
                       RMAX = 2.0D0
                       I1 = NQNYH + 1
                       DO 460 JB = 1, NQ
                          I1 = I1 - NYH
                          DO 450 I = I1, NQNYH
                             YH1(I) = YH1(I) - YH1(I+NYH)
  450                         CONTINUE
  460                      CONTINUE
                       if (ABS(H)  >  HMIN*1.00001D0) go to 470
                          KFLAG = -2
!     ........................EXIT
                          go to 690
  470                      CONTINUE
                       if (NCF  /=  10) go to 480
                          KFLAG = -2
!     ........................EXIT
                          go to 690
  480                      CONTINUE
                       RH = 0.25D0
                       IPUP = MITER
                       IREDO = 1
!                 .........EXIT
                       go to 650
  490                   CONTINUE
                    EXSM = 1.0D0/L
                    RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
                    RHDN = 0.0D0
                    if (NQ  ==  1) go to 500
                       DDN = DVNRMS(N,YH(1,L),EWT)/TESCO(1,NQ)
                       EXDN = 1.0D0/NQ
                       RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
  500                   CONTINUE
                    if (RHSM  >=  RHUP) go to 550
                       if (RHUP  <=  RHDN) go to 540
                          NEWQ = L
                          RH = RHUP
                          if (RH  >=  1.1D0) go to 520
                             IALTH = 3
                             R = 1.0D0/TESCO(2,NQU)
                             DO 510 I = 1, N
                                ACOR(I) = ACOR(I)*R
  510                            CONTINUE
!     ...........................EXIT
                             go to 690
  520                         CONTINUE
                          R = EL(L)/L
                          DO 530 I = 1, N
                             YH(I,NEWQ+1) = ACOR(I)*R
  530                         CONTINUE
                          NQ = NEWQ
                          L = NQ + 1
                          IRET = 2
!           ..................EXIT
                          go to 680
  540                      CONTINUE
                    go to 580
  550                   CONTINUE
                    if (RHSM  <  RHDN) go to 580
                       NEWQ = NQ
                       RH = RHSM
                       if (KFLAG  ==  0 .AND. RH  <  1.1D0) &
                          go to 560
                          if (KFLAG  <=  -2) RH = MIN(RH,0.2D0)
!                             ------------------------------------------
!                              if THERE IS A CHANGE OF ORDER, RESET NQ,
!                              L, AND THE COEFFICIENTS.  IN ANY CASE H
!                              IS RESET ACCORDING TO RH AND THE YH ARRAY
!                              IS RESCALED.  THEN EXIT FROM 680 if THE
!                              STEP WAS OK, OR REDO THE STEP OTHERWISE.
!                             ------------------------------------------
!                 ............EXIT
                          if (NEWQ  ==  NQ) go to 650
                          NQ = NEWQ
                          L = NQ + 1
                          IRET = 2
!           ..................EXIT
                          go to 680
  560                      CONTINUE
                       IALTH = 3
                       R = 1.0D0/TESCO(2,NQU)
                       DO 570 I = 1, N
                          ACOR(I) = ACOR(I)*R
  570                      CONTINUE
!     .....................EXIT
                       go to 690
  580                   CONTINUE
                    NEWQ = NQ - 1
                    RH = RHDN
                    if (KFLAG  <  0 .AND. RH  >  1.0D0) RH = 1.0D0
                    if (KFLAG  ==  0 .AND. RH  <  1.1D0) go to 590
                       if (KFLAG  <=  -2) RH = MIN(RH,0.2D0)
!                          ---------------------------------------------
!                           if THERE IS A CHANGE OF ORDER, RESET NQ, L,
!                           AND THE COEFFICIENTS.  IN ANY CASE H IS
!                           RESET ACCORDING TO RH AND THE YH ARRAY IS
!                           RESCALED.  THEN EXIT FROM 680 if THE STEP
!                           WAS OK, OR REDO THE STEP OTHERWISE.
!                          ---------------------------------------------
!                 .........EXIT
                       if (NEWQ  ==  NQ) go to 650
                       NQ = NEWQ
                       L = NQ + 1
                       IRET = 2
!           ...............EXIT
                       go to 680
  590                   CONTINUE
                    IALTH = 3
                    R = 1.0D0/TESCO(2,NQU)
                    DO 600 I = 1, N
                       ACOR(I) = ACOR(I)*R
  600                   CONTINUE
!     ..................EXIT
                    go to 690
  610                CONTINUE
!                    ---------------------------------------------------
!                     CONTROL REACHES THIS SECTION if 3 OR MORE FAILURES
!                     HAVE OCCURRED.  if 10 FAILURES HAVE OCCURRED, EXIT
!                     WITH KFLAG = -1.  IT IS ASSUMED THAT THE
!                     DERIVATIVES THAT HAVE ACCUMULATED IN THE YH ARRAY
!                     HAVE ERRORS OF THE WRONG ORDER.  HENCE THE FIRST
!                     DERIVATIVE IS RECOMPUTED, AND THE ORDER IS SET TO
!                     1.  THEN H IS REDUCED BY A FACTOR OF 10, AND THE
!                     STEP IS RETRIED, UNTIL IT SUCCEEDS OR H REACHES
!                     HMIN.
!                    ---------------------------------------------------
                 if (KFLAG  /=  -10) go to 620
!                       ------------------------------------------------
!                        ALL RETURNS ARE MADE THROUGH THIS SECTION.  H
!                        IS SAVED IN HOLD TO ALLOW THE CALLER TO CHANGE
!                        H ON THE NEXT STEP.
!                       ------------------------------------------------
                    KFLAG = -1
!     ..................EXIT
                    go to 690
  620                CONTINUE
                 RH = 0.1D0
                 RH = MAX(HMIN/ABS(H),RH)
                 H = H*RH
                 DO 630 I = 1, N
                    Y(I) = YH(I,1)
  630                CONTINUE
                 call DF(TN,Y,SAVF,RPAR,IPAR)
                 NFE = NFE + 1
                 DO 640 I = 1, N
                    YH(I,2) = H*SAVF(I)
  640                CONTINUE
                 IPUP = MITER
                 IALTH = 5
!              ......EXIT
                 if (NQ  /=  1) go to 670
              go to 170
  650             CONTINUE
  660             CONTINUE
              RH = MAX(RH,HMIN/ABS(H))
           go to 110
  670          CONTINUE
           NQ = 1
           L = 2
           IRET = 3
  680       CONTINUE
     go to 70
  690 CONTINUE
  HOLD = H
  JSTART = 1
  return
!     ----------------------- END OF SUBROUTINE DSTOD
!     -----------------------
end
