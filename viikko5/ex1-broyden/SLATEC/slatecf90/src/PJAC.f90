subroutine PJAC (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM, F, &
     JAC, RPAR, IPAR)
!
!! PJAC is subsidiary to DEBDF.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PJAC-S, DPJAC-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   PJAC sets up the iteration matrix (involving the Jacobian) for the
!   integration package DEBDF.
!
!***SEE ALSO  DEBDF
!***ROUTINES CALLED  SGBFA, SGEFA, VNWRMS
!***COMMON BLOCKS    DEBDF1
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!   920422  Changed DIMENSION statement.  (WRB)
!***END PROLOGUE  PJAC
!
!LLL. OPTIMIZE
  INTEGER NEQ, NYH, IWM, I, I1, I2, IER, II, IOWND, IOWNS, J, J1, &
     JJ, JSTART, KFLAG, L, LENP, MAXORD, MBA, MBAND, MEB1, MEBAND, &
     METH, MITER, ML, ML3, MU, N, NFE, NJE, NQ, NQU, NST
  EXTERNAL F, JAC
  REAL Y, YH, EWT, FTEM, SAVF, WM, &
     ROWND, ROWNS, EL0, H, HMIN, HMXI, HU, TN, UROUND, &
     CON, DI, FAC, HL0, R, R0, SRUR, YI, YJ, YJJ, VNWRMS
  DIMENSION         Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*), &
     WM(*), IWM(*), RPAR(*), IPAR(*)
  COMMON /DEBDF1/ ROWND, ROWNS(210), &
     EL0, H, HMIN, HMXI, HU, TN, UROUND, IOWND(14), IOWNS(6), &
     IER, JSTART, KFLAG, L, METH, MITER, MAXORD, N, NQ, NST, NFE, &
     NJE, NQU
!-----------------------------------------------------------------------
! PJAC IS CALLED BY STOD  TO COMPUTE AND PROCESS THE MATRIX
! P = I - H*EL(1)*J , WHERE J IS AN APPROXIMATION TO THE JACOBIAN.
! HERE J IS COMPUTED BY THE USER-SUPPLIED ROUTINE JAC IF
! MITER = 1 OR 4, OR BY FINITE DIFFERENCING if MITER = 2, 3, OR 5.
! if MITER = 3, A DIAGONAL APPROXIMATION TO J IS USED.
! J IS STORED IN WM AND REPLACED BY P.  if MITER  /=  3, P IS THEN
! SUBJECTED TO LU DECOMPOSITION IN PREPARATION FOR LATER SOLUTION
! OF LINEAR SYSTEMS WITH P AS COEFFICIENT MATRIX. THIS IS DONE
! BY SGEFA if MITER = 1 OR 2, AND BY SGBFA IF MITER = 4 OR 5.
!
! IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION
! WITH PJAC USES THE FOLLOWING..
! Y    = ARRAY CONTAINING PREDICTED VALUES ON ENTRY.
! FTEM = WORK ARRAY OF LENGTH N (ACOR IN STOD ).
! SAVF = ARRAY CONTAINING F EVALUATED AT PREDICTED Y.
! WM   = REAL WORK SPACE FOR MATRICES.  ON OUTPUT IT CONTAINS THE
!        INVERSE DIAGONAL MATRIX if MITER = 3 AND THE LU DECOMPOSITION
!        OF P if MITER IS 1, 2 , 4, OR 5.
!        STORAGE OF MATRIX ELEMENTS STARTS AT WM(3).
!        WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA..
!        WM(1) = SQRT(UROUND), USED IN NUMERICAL JACOBIAN INCREMENTS.
!        WM(2) = H*EL0, SAVED FOR LATER USE if MITER = 3.
! IWM  = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING AT
!        IWM(21), if MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS THE
!        BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) if MITER IS 4 OR 5.
! EL0  = EL(1) (INPUT).
! IER  = OUTPUT ERROR FLAG,  = 0 if NO TROUBLE,  /=  0 IF
!        P MATRIX FOUND TO BE SINGULAR.
! THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, TN, UROUND,
! MITER, N, NFE, AND NJE.
!-----------------------------------------------------------------------
!***FIRST EXECUTABLE STATEMENT  PJAC
  NJE = NJE + 1
  HL0 = H*EL0
  go to (100, 200, 300, 400, 500), MITER
! if MITER = 1, call JAC AND MULTIPLY BY SCALAR. -----------------------
 100  LENP = N*N
  DO 110 I = 1,LENP
 110    WM(I+2) = 0.0E0
  call JAC (TN, Y, WM(3), N, RPAR, IPAR)
  CON = -HL0
  DO 120 I = 1,LENP
 120    WM(I+2) = WM(I+2)*CON
  go to 240
! if MITER = 2, MAKE N CALLS TO F TO APPROXIMATE J. --------------------
 200  FAC = VNWRMS (N, SAVF, EWT)
  R0 = 1000.0E0*ABS(H)*UROUND*N*FAC
  if (R0  ==  0.0E0) R0 = 1.0E0
  SRUR = WM(1)
  J1 = 2
  DO 230 J = 1,N
    YJ = Y(J)
    R = MAX(SRUR*ABS(YJ),R0*EWT(J))
    Y(J) = Y(J) + R
    FAC = -HL0/R
    call F (TN, Y, FTEM, RPAR, IPAR)
    DO 220 I = 1,N
 220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
    Y(J) = YJ
    J1 = J1 + N
 230    CONTINUE
  NFE = NFE + N
! ADD IDENTITY MATRIX. -------------------------------------------------
 240  J = 3
  DO 250 I = 1,N
    WM(J) = WM(J) + 1.0E0
 250    J = J + (N + 1)
! DO LU DECOMPOSITION ON P. --------------------------------------------
  call SGEFA (WM(3), N, N, IWM(21), IER)
  return
! if MITER = 3, CONSTRUCT A DIAGONAL APPROXIMATION TO J AND P. ---------
 300  WM(2) = HL0
  IER = 0
  R = EL0*0.1E0
  DO 310 I = 1,N
 310    Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
  call F (TN, Y, WM(3), RPAR, IPAR)
  NFE = NFE + 1
  DO 320 I = 1,N
    R0 = H*SAVF(I) - YH(I,2)
    DI = 0.1E0*R0 - H*(WM(I+2) - SAVF(I))
    WM(I+2) = 1.0E0
    if (ABS(R0)  <  UROUND*EWT(I)) go to 320
    if (ABS(DI)  ==  0.0E0) go to 330
    WM(I+2) = 0.1E0*R0/DI
 320    CONTINUE
  return
 330  IER = -1
  return
! if MITER = 4, call JAC AND MULTIPLY BY SCALAR. -----------------------
 400  ML = IWM(1)
  MU = IWM(2)
  ML3 =  3
  MBAND = ML + MU + 1
  MEBAND = MBAND + ML
  LENP = MEBAND*N
  DO 410 I = 1,LENP
 410    WM(I+2) = 0.0E0
  call JAC (TN, Y, WM(ML3), MEBAND, RPAR, IPAR)
  CON = -HL0
  DO 420 I = 1,LENP
 420    WM(I+2) = WM(I+2)*CON
  go to 570
! if MITER = 5, MAKE MBAND CALLS TO F TO APPROXIMATE J. ----------------
 500  ML = IWM(1)
  MU = IWM(2)
  MBAND = ML + MU + 1
  MBA = MIN(MBAND,N)
  MEBAND = MBAND + ML
  MEB1 = MEBAND - 1
  SRUR = WM(1)
  FAC = VNWRMS (N, SAVF, EWT)
  R0 = 1000.0E0*ABS(H)*UROUND*N*FAC
  if (R0  ==  0.0E0) R0 = 1.0E0
  DO 560 J = 1,MBA
    DO 530 I = J,N,MBAND
      YI = Y(I)
      R = MAX(SRUR*ABS(YI),R0*EWT(I))
 530      Y(I) = Y(I) + R
    call F (TN, Y, FTEM, RPAR, IPAR)
    DO 550 JJ = J,N,MBAND
      Y(JJ) = YH(JJ,1)
      YJJ = Y(JJ)
      R = MAX(SRUR*ABS(YJJ),R0*EWT(JJ))
      FAC = -HL0/R
      I1 = MAX(JJ-MU,1)
      I2 = MIN(JJ+ML,N)
      II = JJ*MEB1 - ML + 2
      DO 540 I = I1,I2
 540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
 550      CONTINUE
 560    CONTINUE
  NFE = NFE + MBA
! ADD IDENTITY MATRIX. -------------------------------------------------
 570  II = MBAND + 2
  DO 580 I = 1,N
    WM(II) = WM(II) + 1.0E0
 580    II = II + MEBAND
! DO LU DECOMPOSITION OF P. --------------------------------------------
  call SGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)
  return
!----------------------- END OF SUBROUTINE PJAC -----------------------
end
