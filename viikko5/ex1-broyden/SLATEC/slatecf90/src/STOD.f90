subroutine STOD (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR, WM, IWM, &
     F, JAC, RPAR, IPAR)
!
!! STOD integrates a system of first order ODE's over one step for DEBDF.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DEBDF
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (STOD-S, DSTOD-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   STOD integrates a system of first order odes over one step in the
!   integrator package DEBDF.
! ----------------------------------------------------------------------
! STOD  performs one step of the integration of an initial value
! problem for a system of ordinary differential equations.
! Note.. STOD  is independent of the value of the iteration method
! indicator MITER, when this is  /=  0, and hence is independent
! of the type of chord method used, or the Jacobian structure.
! Communication with STOD  is done with the following variables..
!
! Y      = An array of length  >=  n used as the Y argument in
!          all calls to F and JAC.
! NEQ    = Integer array containing problem size in NEQ(1), and
!          passed as the NEQ argument in all calls to F and JAC.
! YH     = An NYH by LMAX array containing the dependent variables
!          and their approximate scaled derivatives, where
!          LMAX = MAXORD + 1.  YH(I,J+1) contains the approximate
!          J-th derivative of Y(I), scaled by H**J/Factorial(j)
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
! WM,IWM = Real and integer work arrays associated with matrix
!          operations in chord iteration (MITER  /=  0).
! PJAC   = Name of routine to evaluate and preprocess Jacobian matrix
!          if a chord method is being used.
! SLVS   = Name of routine to solve linear system in chord iteration.
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
!***SEE ALSO  DEBDF
!***ROUTINES CALLED  CFOD, PJAC, SLVS, VNWRMS
!***COMMON BLOCKS    DEBDF1
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!   920422  Changed DIMENSION statement.  (WRB)
!***END PROLOGUE  STOD
  EXTERNAL F, JAC
!
!LLL. OPTIMIZE
  INTEGER NEQ, NYH, IWM, I, I1, IALTH, IER, IOWND, IREDO, IRET, &
     IPUP, J, JB, JSTART, KFLAG, L, LMAX, M, MAXORD, MEO, METH, &
     MITER, N, NCF, NEWQ, NFE, NJE, NQ, NQNYH, NQU, NST, NSTEPJ
  REAL Y, YH, YH1, EWT, SAVF, ACOR, WM, &
     ROWND, CONIT, CRATE, EL, ELCO, HOLD, RC, RMAX, TESCO, &
     EL0, H, HMIN, HMXI, HU, TN, UROUND, &
     DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP, &
     R, RH, RHDN, RHSM, RHUP, TOLD, VNWRMS
  DIMENSION         Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*), &
     ACOR(*), WM(*), IWM(*), RPAR(*), IPAR(*)
  COMMON /DEBDF1/ ROWND, CONIT, CRATE, EL(13), ELCO(13,12), &
     HOLD, RC, RMAX, TESCO(3,12), &
     EL0, H, HMIN, HMXI, HU, TN, UROUND, IOWND(7), KSTEPS, IOD(6), &
     IALTH, IPUP, LMAX, MEO, NQNYH, NSTEPJ, &
     IER, JSTART, KFLAG, L, METH, MITER, MAXORD, N, NQ, NST, NFE, &
     NJE, NQU
!
!
!***FIRST EXECUTABLE STATEMENT  STOD
  KFLAG = 0
  TOLD = TN
  NCF = 0
  if (JSTART  >  0) go to 200
  if (JSTART  ==  -1) go to 100
  if (JSTART  ==  -2) go to 160
!-----------------------------------------------------------------------
! ON THE FIRST CALL, THE ORDER IS SET TO 1, AND OTHER VARIABLES ARE
! INITIALIZED.  RMAX IS THE MAXIMUM RATIO BY WHICH H CAN BE INCREASED
! IN A SINGLE STEP.  IT IS INITIALLY 1.E4 TO COMPENSATE FOR THE SMALL
! INITIAL H, BUT THEN IS NORMALLY EQUAL TO 10.  if A FAILURE
! OCCURS (IN CORRECTOR CONVERGENCE OR ERROR TEST), RMAX IS SET AT 2
! FOR THE NEXT INCREASE.
!-----------------------------------------------------------------------
  LMAX = MAXORD + 1
  NQ = 1
  L = 2
  IALTH = 2
  RMAX = 10000.0E0
  RC = 0.0E0
  EL0 = 1.0E0
  CRATE = 0.7E0
  DELP = 0.0E0
  HOLD = H
  MEO = METH
  NSTEPJ = 0
  IRET = 3
  go to 140
!-----------------------------------------------------------------------
! THE FOLLOWING BLOCK HANDLES PRELIMINARIES NEEDED WHEN JSTART = -1.
! IPUP IS SET TO MITER TO FORCE A MATRIX UPDATE.
! if AN ORDER INCREASE IS ABOUT TO BE CONSIDERED (IALTH = 1),
! IALTH IS RESET TO 2 TO POSTPONE CONSIDERATION ONE MORE STEP.
! if THE CALLER HAS CHANGED METH, CFOD  IS CALLED TO RESET
! THE COEFFICIENTS OF THE METHOD.
! if THE CALLER HAS CHANGED MAXORD TO A VALUE LESS THAN THE CURRENT
! ORDER NQ, NQ IS REDUCED TO MAXORD, AND A NEW H CHOSEN ACCORDINGLY.
! if H IS TO BE CHANGED, YH MUST BE RESCALED.
! if H OR METH IS BEING CHANGED, IALTH IS RESET TO L = NQ + 1
! TO PREVENT FURTHER CHANGES IN H FOR THAT MANY STEPS.
!-----------------------------------------------------------------------
 100  IPUP = MITER
  LMAX = MAXORD + 1
  if (IALTH  ==  1) IALTH = 2
  if (METH  ==  MEO) go to 110
  call CFOD  (METH, ELCO, TESCO)
  MEO = METH
  if (NQ  >  MAXORD) go to 120
  IALTH = L
  IRET = 1
  go to 150
 110  if (NQ  <=  MAXORD) go to 160
 120  NQ = MAXORD
  L = LMAX
  DO 125 I = 1,L
 125    EL(I) = ELCO(I,NQ)
  NQNYH = NQ*NYH
  RC = RC*EL(1)/EL0
  EL0 = EL(1)
  CONIT = 0.5E0/(NQ+2)
  DDN = VNWRMS (N, SAVF, EWT)/TESCO(1,L)
  EXDN = 1.0E0/L
  RHDN = 1.0E0/(1.3E0*DDN**EXDN + 0.0000013E0)
  RH = MIN(RHDN,1.0E0)
  IREDO = 3
  if (H  ==  HOLD) go to 170
  RH = MIN(RH,ABS(H/HOLD))
  H = HOLD
  go to 175
!-----------------------------------------------------------------------
! CFOD  IS CALLED TO GET ALL THE INTEGRATION COEFFICIENTS FOR THE
! CURRENT METH.  THEN THE EL VECTOR AND RELATED CONSTANTS ARE RESET
! WHENEVER THE ORDER NQ IS CHANGED, OR AT THE START OF THE PROBLEM.
!-----------------------------------------------------------------------
 140  call CFOD  (METH, ELCO, TESCO)
 150  DO 155 I = 1,L
 155    EL(I) = ELCO(I,NQ)
  NQNYH = NQ*NYH
  RC = RC*EL(1)/EL0
  EL0 = EL(1)
  CONIT = 0.5E0/(NQ+2)
  go to (160, 170, 200), IRET
!-----------------------------------------------------------------------
! if H IS BEING CHANGED, THE H RATIO RH IS CHECKED AGAINST
! RMAX, HMIN, AND HMXI, AND THE YH ARRAY RESCALED.  IALTH IS SET TO
! L = NQ + 1 TO PREVENT A CHANGE OF H FOR THAT MANY STEPS, UNLESS
! FORCED BY A CONVERGENCE OR ERROR TEST FAILURE.
!-----------------------------------------------------------------------
 160  if (H  ==  HOLD) go to 200
  RH = H/HOLD
  H = HOLD
  IREDO = 3
  go to 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,RMAX)
  RH = RH/MAX(1.0E0,ABS(H)*HMXI*RH)
  R = 1.0E0
  DO 180 J = 2,L
    R = R*RH
    DO 180 I = 1,N
 180      YH(I,J) = YH(I,J)*R
  H = H*RH
  RC = RC*RH
  IALTH = L
  if (IREDO  ==  0) go to 680
!-----------------------------------------------------------------------
! THIS SECTION COMPUTES THE PREDICTED VALUES BY EFFECTIVELY
! MULTIPLYING THE YH ARRAY BY THE PASCAL TRIANGLE MATRIX.
! RC IS THE RATIO OF NEW TO OLD VALUES OF THE COEFFICIENT  H*EL(1).
! WHEN RC DIFFERS FROM 1 BY MORE THAN 30 PERCENT, IPUP IS SET TO MITER
! TO FORCE PJAC TO BE CALLED, if A JACOBIAN IS INVOLVED.
! IN ANY CASE, PJAC IS CALLED AT LEAST EVERY 20-TH STEP.
!-----------------------------------------------------------------------
 200  if (ABS(RC-1.0E0)  >  0.3E0) IPUP = MITER
  if (NST  >=  NSTEPJ+20) IPUP = MITER
  TN = TN + H
  I1 = NQNYH + 1
  DO 215 JB = 1,NQ
    I1 = I1 - NYH
    DO 210 I = I1,NQNYH
 210      YH1(I) = YH1(I) + YH1(I+NYH)
 215    CONTINUE
  KSTEPS = KSTEPS + 1
!-----------------------------------------------------------------------
! UP TO 3 CORRECTOR ITERATIONS ARE TAKEN.  A CONVERGENCE TEST IS
! MADE ON THE R.M.S. NORM OF EACH CORRECTION, WEIGHTED BY THE ERROR
! WEIGHT VECTOR EWT.  THE SUM OF THE CORRECTIONS IS ACCUMULATED IN THE
! VECTOR ACOR(I).  THE YH ARRAY IS NOT ALTERED IN THE CORRECTOR LOOP.
!-----------------------------------------------------------------------
 220  M = 0
  DO 230 I = 1,N
 230    Y(I) = YH(I,1)
  call F (TN, Y, SAVF, RPAR, IPAR)
  NFE = NFE + 1
  if (IPUP  <=  0) go to 250
!-----------------------------------------------------------------------
! if INDICATED, THE MATRIX P = I - H*EL(1)*J IS REEVALUATED AND
! PREPROCESSED BEFORE STARTING THE CORRECTOR ITERATION.  IPUP IS SET
! TO 0 AS AN INDICATOR THAT THIS HAS BEEN DONE.
!-----------------------------------------------------------------------
  IPUP = 0
  RC = 1.0E0
  NSTEPJ = NST
  CRATE = 0.7E0
  call PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC, &
            RPAR, IPAR)
  if (IER  /=  0) go to 430
 250  DO 260 I = 1,N
 260    ACOR(I) = 0.0E0
 270  if (MITER  /=  0) go to 350
!-----------------------------------------------------------------------
! IN THE CASE OF FUNCTIONAL ITERATION, UPDATE Y DIRECTLY FROM
! THE RESULT OF THE LAST FUNCTION EVALUATION.
!-----------------------------------------------------------------------
  DO 290 I = 1,N
    SAVF(I) = H*SAVF(I) - YH(I,2)
 290    Y(I) = SAVF(I) - ACOR(I)
  DEL = VNWRMS (N, Y, EWT)
  DO 300 I = 1,N
    Y(I) = YH(I,1) + EL(1)*SAVF(I)
 300    ACOR(I) = SAVF(I)
  go to 400
!-----------------------------------------------------------------------
! IN THE CASE OF THE CHORD METHOD, COMPUTE THE CORRECTOR ERROR,
! AND SOLVE THE LINEAR SYSTEM WITH THAT AS RIGHT-HAND SIDE AND
! P AS COEFFICIENT MATRIX.
!-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
  call SLVS (WM, IWM, Y, SAVF)
  if (IER  /=  0) go to 410
  DEL = VNWRMS (N, Y, EWT)
  DO 380 I = 1,N
    ACOR(I) = ACOR(I) + Y(I)
 380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
!-----------------------------------------------------------------------
! TEST FOR CONVERGENCE.  if M > 0, AN ESTIMATE OF THE CONVERGENCE
! RATE CONSTANT IS STORED IN CRATE, AND THIS IS USED IN THE TEST.
!-----------------------------------------------------------------------
 400  if (M  /=  0) CRATE = MAX(0.2E0*CRATE,DEL/DELP)
  DCON = DEL*MIN(1.0E0,1.5E0*CRATE)/(TESCO(2,NQ)*CONIT)
  if (DCON  <=  1.0E0) go to 450
  M = M + 1
  if (M  ==  3) go to 410
  if (M  >=  2 .AND. DEL  >  2.0E0*DELP) go to 410
  DELP = DEL
  call F (TN, Y, SAVF, RPAR, IPAR)
  NFE = NFE + 1
  go to 270
!-----------------------------------------------------------------------
! THE CORRECTOR ITERATION FAILED TO CONVERGE IN 3 TRIES.
! if MITER  /=  0 AND THE JACOBIAN IS OUT OF DATE, PJAC IS CALLED FOR
! THE NEXT TRY.  OTHERWISE THE YH ARRAY IS RETRACTED TO ITS VALUES
! BEFORE PREDICTION, AND H IS REDUCED, if POSSIBLE.  IF H CANNOT BE
! REDUCED OR 10 FAILURES HAVE OCCURRED, EXIT WITH KFLAG = -2.
!-----------------------------------------------------------------------
 410  if (IPUP  ==  0) go to 430
  IPUP = MITER
  go to 220
 430  TN = TOLD
  NCF = NCF + 1
  RMAX = 2.0E0
  I1 = NQNYH + 1
  DO 445 JB = 1,NQ
    I1 = I1 - NYH
    DO 440 I = I1,NQNYH
 440      YH1(I) = YH1(I) - YH1(I+NYH)
 445    CONTINUE
  if (ABS(H)  <=  HMIN*1.00001E0) go to 670
  if (NCF  ==  10) go to 670
  RH = 0.25E0
  IPUP = MITER
  IREDO = 1
  go to 170
!-----------------------------------------------------------------------
! THE CORRECTOR HAS CONVERGED.  IPUP IS SET TO -1 if MITER  /=  0,
! TO SIGNAL THAT THE JACOBIAN INVOLVED MAY NEED UPDATING LATER.
! THE LOCAL ERROR TEST IS MADE AND CONTROL PASSES TO STATEMENT 500
! if IT FAILS.
!-----------------------------------------------------------------------
 450  if (MITER  /=  0) IPUP = -1
  if (M  ==  0) DSM = DEL/TESCO(2,NQ)
  if (M  >  0) DSM = VNWRMS (N, ACOR, EWT)/TESCO(2,NQ)
  if (DSM  >  1.0E0) go to 500
!-----------------------------------------------------------------------
! AFTER A SUCCESSFUL STEP, UPDATE THE YH ARRAY.
! CONSIDER CHANGING H if IALTH = 1.  OTHERWISE DECREASE IALTH BY 1.
! if IALTH IS THEN 1 AND NQ  <  MAXORD, THEN ACOR IS SAVED FOR
! USE IN A POSSIBLE ORDER INCREASE ON THE NEXT STEP.
! if A CHANGE IN H IS CONSIDERED, AN INCREASE OR DECREASE IN ORDER
! BY ONE IS CONSIDERED ALSO.  A CHANGE IN H IS MADE ONLY if IT IS BY A
! FACTOR OF AT LEAST 1.1.  if NOT, IALTH IS SET TO 3 TO PREVENT
! TESTING FOR THAT MANY STEPS.
!-----------------------------------------------------------------------
  KFLAG = 0
  IREDO = 0
  NST = NST + 1
  HU = H
  NQU = NQ
  DO 470 J = 1,L
    DO 470 I = 1,N
 470      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
  IALTH = IALTH - 1
  if (IALTH  ==  0) go to 520
  if (IALTH  >  1) go to 690
  if (L  ==  LMAX) go to 690
  DO 490 I = 1,N
 490    YH(I,LMAX) = ACOR(I)
  go to 690
!-----------------------------------------------------------------------
! THE ERROR TEST FAILED.  KFLAG KEEPS TRACK OF MULTIPLE FAILURES.
! RESTORE TN AND THE YH ARRAY TO THEIR PREVIOUS VALUES, AND PREPARE
! TO TRY THE STEP AGAIN.  COMPUTE THE OPTIMUM STEP SIZE FOR THIS OR
! ONE LOWER ORDER.  AFTER 2 OR MORE FAILURES, H IS FORCED TO DECREASE
! BY A FACTOR OF 0.2 OR LESS.
!-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
  TN = TOLD
  I1 = NQNYH + 1
  DO 515 JB = 1,NQ
    I1 = I1 - NYH
    DO 510 I = I1,NQNYH
 510      YH1(I) = YH1(I) - YH1(I+NYH)
 515    CONTINUE
  RMAX = 2.0E0
  if (ABS(H)  <=  HMIN*1.00001E0) go to 660
  if (KFLAG  <=  -3) go to 640
  IREDO = 2
  RHUP = 0.0E0
  go to 540
!-----------------------------------------------------------------------
! REGARDLESS OF THE SUCCESS OR FAILURE OF THE STEP, FACTORS
! RHDN, RHSM, AND RHUP ARE COMPUTED, BY WHICH H COULD BE MULTIPLIED
! AT ORDER NQ - 1, ORDER NQ, OR ORDER NQ + 1, RESPECTIVELY.
! IN THE CASE OF FAILURE, RHUP = 0.0 TO AVOID AN ORDER INCREASE.
! THE LARGEST OF THESE IS DETERMINED AND THE NEW ORDER CHOSEN
! ACCORDINGLY.  if THE ORDER IS TO BE INCREASED, WE COMPUTE ONE
! ADDITIONAL SCALED DERIVATIVE.
!-----------------------------------------------------------------------
 520  RHUP = 0.0E0
  if (L  ==  LMAX) go to 540
  DO 530 I = 1,N
 530    SAVF(I) = ACOR(I) - YH(I,LMAX)
  DUP = VNWRMS (N, SAVF, EWT)/TESCO(3,NQ)
  EXUP = 1.0E0/(L+1)
  RHUP = 1.0E0/(1.4E0*DUP**EXUP + 0.0000014E0)
 540  EXSM = 1.0E0/L
  RHSM = 1.0E0/(1.2E0*DSM**EXSM + 0.0000012E0)
  RHDN = 0.0E0
  if (NQ  ==  1) go to 560
  DDN = VNWRMS (N, YH(1,L), EWT)/TESCO(1,NQ)
  EXDN = 1.0E0/NQ
  RHDN = 1.0E0/(1.3E0*DDN**EXDN + 0.0000013E0)
 560  if (RHSM  >=  RHUP) go to 570
  if (RHUP  >  RHDN) go to 590
  go to 580
 570  if (RHSM  <  RHDN) go to 580
  NEWQ = NQ
  RH = RHSM
  go to 620
 580  NEWQ = NQ - 1
  RH = RHDN
  if (KFLAG  <  0 .AND. RH  >  1.0E0) RH = 1.0E0
  go to 620
 590  NEWQ = L
  RH = RHUP
  if (RH  <  1.1E0) go to 610
  R = EL(L)/L
  DO 600 I = 1,N
 600    YH(I,NEWQ+1) = ACOR(I)*R
  go to 630
 610  IALTH = 3
  go to 690
 620  if ((KFLAG  ==  0) .AND. (RH  <  1.1E0)) go to 610
  if (KFLAG  <=  -2) RH = MIN(RH,0.2E0)
!-----------------------------------------------------------------------
! if THERE IS A CHANGE OF ORDER, RESET NQ, L, AND THE COEFFICIENTS.
! IN ANY CASE H IS RESET ACCORDING TO RH AND THE YH ARRAY IS RESCALED.
! THEN EXIT FROM 680 if THE STEP WAS OK, OR REDO THE STEP OTHERWISE.
!-----------------------------------------------------------------------
  if (NEWQ  ==  NQ) go to 170
 630  NQ = NEWQ
  L = NQ + 1
  IRET = 2
  go to 150
!-----------------------------------------------------------------------
! CONTROL REACHES THIS SECTION if 3 OR MORE FAILURES HAVE OCCURRED.
! if 10 FAILURES HAVE OCCURRED, EXIT WITH KFLAG = -1.
! IT IS ASSUMED THAT THE DERIVATIVES THAT HAVE ACCUMULATED IN THE
! YH ARRAY HAVE ERRORS OF THE WRONG ORDER.  HENCE THE FIRST
! DERIVATIVE IS RECOMPUTED, AND THE ORDER IS SET TO 1.  THEN
! H IS REDUCED BY A FACTOR OF 10, AND THE STEP IS RETRIED,
! UNTIL IT SUCCEEDS OR H REACHES HMIN.
!-----------------------------------------------------------------------
 640  if (KFLAG  ==  -10) go to 660
  RH = 0.1E0
  RH = MAX(HMIN/ABS(H),RH)
  H = H*RH
  DO 645 I = 1,N
 645    Y(I) = YH(I,1)
  call F (TN, Y, SAVF, RPAR, IPAR)
  NFE = NFE + 1
  DO 650 I = 1,N
 650    YH(I,2) = H*SAVF(I)
  IPUP = MITER
  IALTH = 5
  if (NQ  ==  1) go to 200
  NQ = 1
  L = 2
  IRET = 3
  go to 150
!-----------------------------------------------------------------------
! ALL RETURNS ARE MADE THROUGH THIS SECTION.  H IS SAVED IN HOLD
! TO ALLOW THE CALLER TO CHANGE H ON THE NEXT STEP.
!-----------------------------------------------------------------------
 660  KFLAG = -1
  go to 700
 670  KFLAG = -2
  go to 700
 680  RMAX = 10.0E0
 690  R = 1.0E0/TESCO(2,NQU)
  DO 695 I = 1,N
 695    ACOR(I) = ACOR(I)*R
 700  HOLD = H
  JSTART = 1
  return
!----------------------- END OF SUBROUTINE STOD  -----------------------
end
