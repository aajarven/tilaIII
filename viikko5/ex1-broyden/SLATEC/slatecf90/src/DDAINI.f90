subroutine DDAINI (X, Y, YPRIME, NEQ, RES, JAC, H, WT, IDID, RPAR, &
     IPAR, PHI, DELTA, E, WM, IWM, HMIN, UROUND, NONNEG, NTEMP)
!
!! DDAINI is the initialization routine for DDASSL.
!
!***LIBRARY   SLATEC (DASSL)
!***TYPE      DOUBLE PRECISION (SDAINI-S, DDAINI-D)
!***AUTHOR  Petzold, Linda R., (LLNL)
!***DESCRIPTION
!-----------------------------------------------------------------
!     DDAINI TAKES ONE STEP OF SIZE H OR SMALLER
!     WITH THE BACKWARD EULER METHOD, TO
!     FIND YPRIME.  X AND Y ARE UPDATED TO BE CONSISTENT WITH THE
!     NEW STEP.  A MODIFIED DAMPED NEWTON ITERATION IS USED TO
!     SOLVE THE CORRECTOR ITERATION.
!
!     THE INITIAL GUESS FOR YPRIME IS USED IN THE
!     PREDICTION, AND IN FORMING THE ITERATION
!     MATRIX, BUT IS NOT INVOLVED IN THE
!     ERROR TEST. THIS MAY HAVE TROUBLE
!     CONVERGING if THE INITIAL GUESS IS NO
!     GOOD, OR if G(X,Y,YPRIME) DEPENDS
!     NONLINEARLY ON YPRIME.
!
!     THE PARAMETERS REPRESENT:
!     X --         INDEPENDENT VARIABLE
!     Y --         SOLUTION VECTOR AT X
!     YPRIME --    DERIVATIVE OF SOLUTION VECTOR
!     NEQ --       NUMBER OF EQUATIONS
!     H --         STEPSIZE. IMDER MAY USE A STEPSIZE
!                  SMALLER THAN H.
!     WT --        VECTOR OF WEIGHTS FOR ERROR
!                  CRITERION
!     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS
!                  IDID= 1 -- YPRIME WAS FOUND SUCCESSFULLY
!                  IDID=-12 -- DDAINI FAILED TO FIND YPRIME
!     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS
!                  THAT ARE NOT ALTERED BY DDAINI
!     PHI --       WORK SPACE FOR DDAINI
!     DELTA,E --   WORK SPACE FOR DDAINI
!     WM,IWM --    REAL AND INTEGER ARRAYS STORING
!                  MATRIX INFORMATION
!
!-----------------------------------------------------------------
!***ROUTINES CALLED  DDAJAC, DDANRM, DDASLV
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
!   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (FNF)
!   901030  Minor corrections to declarations.  (FNF)
!***END PROLOGUE  DDAINI
!
  INTEGER  NEQ, IDID, IPAR(*), IWM(*), NONNEG, NTEMP
  DOUBLE PRECISION &
     X, Y(*), YPRIME(*), H, WT(*), RPAR(*), PHI(NEQ,*), DELTA(*), &
     E(*), WM(*), HMIN, UROUND
  EXTERNAL  RES, JAC
!
  EXTERNAL  DDAJAC, DDANRM, DDASLV
  DOUBLE PRECISION  DDANRM
!
  INTEGER  I, IER, IRES, JCALC, LNJE, LNRE, M, MAXIT, MJAC, NCF, &
     NEF, NSF
  DOUBLE PRECISION &
     CJ, DAMP, DELNRM, ERR, OLDNRM, R, RATE, S, XOLD, YNORM
  LOGICAL  CONVGD
!
  PARAMETER (LNRE=12)
  PARAMETER (LNJE=13)
!
  DATA MAXIT/10/,MJAC/5/
  DATA DAMP/0.75D0/
!
!
!---------------------------------------------------
!     BLOCK 1.
!     INITIALIZATIONS.
!---------------------------------------------------
!
!***FIRST EXECUTABLE STATEMENT  DDAINI
  IDID=1
  NEF=0
  NCF=0
  NSF=0
  XOLD=X
  YNORM=DDANRM(NEQ,Y,WT,RPAR,IPAR)
!
!     SAVE Y AND YPRIME IN PHI
  DO 100 I=1,NEQ
     PHI(I,1)=Y(I)
100      PHI(I,2)=YPRIME(I)
!
!
!----------------------------------------------------
!     BLOCK 2.
!     DO ONE BACKWARD EULER STEP.
!----------------------------------------------------
!
!     SET UP FOR START OF CORRECTOR ITERATION
200   CJ=1.0D0/H
  X=X+H
!
!     PREDICT SOLUTION AND DERIVATIVE
  DO 250 I=1,NEQ
250     Y(I)=Y(I)+H*YPRIME(I)
!
  JCALC=-1
  M=0
  CONVGD=.TRUE.
!
!
!     CORRECTOR LOOP.
300   IWM(LNRE)=IWM(LNRE)+1
  IRES=0
!
  call RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
  if (IRES < 0) go to 430
!
!
!     EVALUATE THE ITERATION MATRIX
  if (JCALC /= -1) go to 310
  IWM(LNJE)=IWM(LNJE)+1
  JCALC=0
  call DDAJAC(NEQ,X,Y,YPRIME,DELTA,CJ,H, &
     IER,WT,E,WM,IWM,RES,IRES, &
     UROUND,JAC,RPAR,IPAR,NTEMP)
!
  S=1000000.D0
  if (IRES < 0) go to 430
  if (IER /= 0) go to 430
  NSF=0
!
!
!
!     MULTIPLY RESIDUAL BY DAMPING FACTOR
310   CONTINUE
  DO 320 I=1,NEQ
320      DELTA(I)=DELTA(I)*DAMP
!
!     COMPUTE A NEW ITERATE (BACK SUBSTITUTION)
!     STORE THE CORRECTION IN DELTA
!
  call DDASLV(NEQ,DELTA,WM,IWM)
!
!     UPDATE Y AND YPRIME
  DO 330 I=1,NEQ
     Y(I)=Y(I)-DELTA(I)
330      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
!
!     TEST FOR CONVERGENCE OF THE ITERATION.
!
  DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
  if (DELNRM <= 100.D0*UROUND*YNORM) &
     go to 400
!
  if (M > 0) go to 340
     OLDNRM=DELNRM
     go to 350
!
340   RATE=(DELNRM/OLDNRM)**(1.0D0/M)
  if (RATE > 0.90D0) go to 430
  S=RATE/(1.0D0-RATE)
!
350   if (S*DELNRM  <=  0.33D0) go to 400
!
!
!     THE CORRECTOR HAS NOT YET CONVERGED. UPDATE
!     M AND AND TEST WHETHER THE MAXIMUM
!     NUMBER OF ITERATIONS HAVE BEEN TRIED.
!     EVERY MJAC ITERATIONS, GET A NEW
!     ITERATION MATRIX.
!
  M=M+1
  if (M >= MAXIT) go to 430
!
  if ((M/MJAC)*MJAC == M) JCALC=-1
  go to 300
!
!
!     THE ITERATION HAS CONVERGED.
!     CHECK NONNEGATIVITY CONSTRAINTS
400   if (NONNEG == 0) go to 450
  DO 410 I=1,NEQ
410      DELTA(I)=MIN(Y(I),0.0D0)
!
  DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
  if (DELNRM > 0.33D0) go to 430
!
  DO 420 I=1,NEQ
     Y(I)=Y(I)-DELTA(I)
420      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
  go to 450
!
!
!     EXITS FROM CORRECTOR LOOP.
430   CONVGD=.FALSE.
450   if (.NOT.CONVGD) go to 600
!
!
!
!-----------------------------------------------------
!     BLOCK 3.
!     THE CORRECTOR ITERATION CONVERGED.
!     DO ERROR TEST.
!-----------------------------------------------------
!
  DO 510 I=1,NEQ
510      E(I)=Y(I)-PHI(I,1)
  ERR=DDANRM(NEQ,E,WT,RPAR,IPAR)
!
  if (ERR <= 1.0D0) RETURN
!
!
!
!--------------------------------------------------------
!     BLOCK 4.
!     THE BACKWARD EULER STEP FAILED. RESTORE X, Y
!     AND YPRIME TO THEIR ORIGINAL VALUES.
!     REDUCE STEPSIZE AND TRY AGAIN, IF
!     POSSIBLE.
!---------------------------------------------------------
!
600   CONTINUE
  X = XOLD
  DO 610 I=1,NEQ
     Y(I)=PHI(I,1)
610      YPRIME(I)=PHI(I,2)
!
  if (CONVGD) go to 640
  if (IER == 0) go to 620
     NSF=NSF+1
     H=H*0.25D0
     if (NSF < 3.AND.ABS(H) >= HMIN) go to 690
     IDID=-12
     return
620   if (IRES > -2) go to 630
     IDID=-12
     return
630   NCF=NCF+1
  H=H*0.25D0
  if (NCF < 10.AND.ABS(H) >= HMIN) go to 690
     IDID=-12
     return
!
640   NEF=NEF+1
  R=0.90D0/(2.0D0*ERR+0.0001D0)
  R=MAX(0.1D0,MIN(0.5D0,R))
  H=H*R
  if (ABS(H) >= HMIN.AND.NEF < 10) go to 690
     IDID=-12
     return
690      go to 200
!
!-------------END OF SUBROUTINE DDAINI----------------------
end
