subroutine DPJAC (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM, DF, &
     DJAC, RPAR, IPAR)
!
!! DPJAC is subsidiary to DDEBDF.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (PJAC-S, DPJAC-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   DPJAC sets up the iteration matrix (involving the Jacobian) for the
!   integration package DDEBDF.
!
!***SEE ALSO  DDEBDF
!***ROUTINES CALLED  DGBFA, DGEFA, DVNRMS
!***COMMON BLOCKS    DDEBD1
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!   920422  Changed DIMENSION statement.  (WRB)
!***END PROLOGUE  DPJAC
!
  INTEGER I, I1, I2, IER, II, IOWND, IOWNS, IPAR, IWM, J, J1, &
        JJ, JSTART, KFLAG, L, LENP, MAXORD, MBA, MBAND, &
        MEB1, MEBAND, METH, MITER, ML, ML3, MU, N, NEQ, &
        NFE, NJE, NQ, NQU, NST, NYH
  DOUBLE PRECISION CON, DI, DVNRMS, EL0, EWT, &
        FAC, FTEM, H, HL0, HMIN, HMXI, HU, R, R0, ROWND, ROWNS, &
        RPAR, SAVF, SRUR, TN, UROUND, WM, Y, YH, YI, YJ, YJJ
  EXTERNAL DF, DJAC
  DIMENSION Y(*),YH(NYH,*),EWT(*),FTEM(*),SAVF(*),WM(*),IWM(*), &
            RPAR(*),IPAR(*)
  COMMON /DDEBD1/ ROWND,ROWNS(210),EL0,H,HMIN,HMXI,HU,TN,UROUND, &
                  IOWND(14),IOWNS(6),IER,JSTART,KFLAG,L,METH,MITER, &
                  MAXORD,N,NQ,NST,NFE,NJE,NQU
!     ------------------------------------------------------------------
!      DPJAC IS CALLED BY DSTOD  TO COMPUTE AND PROCESS THE MATRIX
!      P = I - H*EL(1)*J , WHERE J IS AN APPROXIMATION TO THE JACOBIAN.
!      HERE J IS COMPUTED BY THE USER-SUPPLIED ROUTINE DJAC IF
!      MITER = 1 OR 4, OR BY FINITE DIFFERENCING if MITER = 2, 3, OR 5.
!      if MITER = 3, A DIAGONAL APPROXIMATION TO J IS USED.
!      J IS STORED IN WM AND REPLACED BY P.  if MITER  /=  3, P IS THEN
!      SUBJECTED TO LU DECOMPOSITION IN PREPARATION FOR LATER SOLUTION
!      OF LINEAR SYSTEMS WITH P AS COEFFICIENT MATRIX. THIS IS DONE
!      BY DGEFA if MITER = 1 OR 2, AND BY DGBFA IF MITER = 4 OR 5.
!
!      IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION
!      WITH DPJAC USES THE FOLLOWING..
!      Y    = ARRAY CONTAINING PREDICTED VALUES ON ENTRY.
!      FTEM = WORK ARRAY OF LENGTH N (ACOR IN DSTOD ).
!      SAVF = ARRAY CONTAINING DF EVALUATED AT PREDICTED Y.
!      WM   = DOUBLE PRECISION WORK SPACE FOR MATRICES.  ON OUTPUT IT
!      CONTAINS THE
!             INVERSE DIAGONAL MATRIX if MITER = 3 AND THE LU
!             DECOMPOSITION OF P if MITER IS 1, 2 , 4, OR 5.
!             STORAGE OF MATRIX ELEMENTS STARTS AT WM(3).
!             WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA..
!             WM(1) = SQRT(UROUND), USED IN NUMERICAL JACOBIAN
!             INCREMENTS.  WM(2) = H*EL0, SAVED FOR LATER USE if MITER =
!             3.
!      IWM  = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING
!             AT IWM(21), if MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS
!             THE BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) if MITER
!             IS 4 OR 5.
!      EL0  = EL(1) (INPUT).
!      IER  = OUTPUT ERROR FLAG,  = 0 if NO TROUBLE,  /=  0 IF
!             P MATRIX FOUND TO BE SINGULAR.
!      THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, TN, UROUND,
!      MITER, N, NFE, AND NJE.
!-----------------------------------------------------------------------
!     BEGIN BLOCK PERMITTING ...EXITS TO 240
!        BEGIN BLOCK PERMITTING ...EXITS TO 220
!           BEGIN BLOCK PERMITTING ...EXITS TO 130
!              BEGIN BLOCK PERMITTING ...EXITS TO 70
!***FIRST EXECUTABLE STATEMENT  DPJAC
              NJE = NJE + 1
              HL0 = H*EL0
              go to (10,40,90,140,170), MITER
!                 if MITER = 1, call DJAC AND MULTIPLY BY SCALAR.
!                 -----------------------
   10             CONTINUE
              LENP = N*N
              DO 20 I = 1, LENP
                 WM(I+2) = 0.0D0
   20             CONTINUE
              call DJAC(TN,Y,WM(3),N,RPAR,IPAR)
              CON = -HL0
              DO 30 I = 1, LENP
                 WM(I+2) = WM(I+2)*CON
   30             CONTINUE
!              ...EXIT
              go to 70
!                 if MITER = 2, MAKE N CALLS TO DF TO APPROXIMATE J.
!                 --------------------
   40             CONTINUE
              FAC = DVNRMS(N,SAVF,EWT)
              R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
              if (R0  ==  0.0D0) R0 = 1.0D0
              SRUR = WM(1)
              J1 = 2
              DO 60 J = 1, N
                 YJ = Y(J)
                 R = MAX(SRUR*ABS(YJ),R0*EWT(J))
                 Y(J) = Y(J) + R
                 FAC = -HL0/R
                 call DF(TN,Y,FTEM,RPAR,IPAR)
                 DO 50 I = 1, N
                    WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
   50                CONTINUE
                 Y(J) = YJ
                 J1 = J1 + N
   60             CONTINUE
              NFE = NFE + N
   70          CONTINUE
!              ADD IDENTITY MATRIX.
!              -------------------------------------------------
           J = 3
           DO 80 I = 1, N
              WM(J) = WM(J) + 1.0D0
              J = J + (N + 1)
   80          CONTINUE
!              DO LU DECOMPOSITION ON P.
!              --------------------------------------------
           call DGEFA(WM(3),N,N,IWM(21),IER)
!     .........EXIT
           go to 240
!              if MITER = 3, CONSTRUCT A DIAGONAL APPROXIMATION TO J AND
!              P. ---------
   90          CONTINUE
           WM(2) = HL0
           IER = 0
           R = EL0*0.1D0
           DO 100 I = 1, N
              Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
  100          CONTINUE
           call DF(TN,Y,WM(3),RPAR,IPAR)
           NFE = NFE + 1
           DO 120 I = 1, N
              R0 = H*SAVF(I) - YH(I,2)
              DI = 0.1D0*R0 - H*(WM(I+2) - SAVF(I))
              WM(I+2) = 1.0D0
              if (ABS(R0)  <  UROUND*EWT(I)) go to 110
!           .........EXIT
                 if (ABS(DI)  ==  0.0D0) go to 130
                 WM(I+2) = 0.1D0*R0/DI
  110             CONTINUE
  120          CONTINUE
!     .........EXIT
           go to 240
  130       CONTINUE
        IER = -1
!     ......EXIT
        go to 240
!           if MITER = 4, call DJAC AND MULTIPLY BY SCALAR.
!           -----------------------
  140       CONTINUE
        ML = IWM(1)
        MU = IWM(2)
        ML3 = 3
        MBAND = ML + MU + 1
        MEBAND = MBAND + ML
        LENP = MEBAND*N
        DO 150 I = 1, LENP
           WM(I+2) = 0.0D0
  150       CONTINUE
        call DJAC(TN,Y,WM(ML3),MEBAND,RPAR,IPAR)
        CON = -HL0
        DO 160 I = 1, LENP
           WM(I+2) = WM(I+2)*CON
  160       CONTINUE
!        ...EXIT
        go to 220
!           if MITER = 5, MAKE MBAND CALLS TO DF TO APPROXIMATE J.
!           ----------------
  170       CONTINUE
        ML = IWM(1)
        MU = IWM(2)
        MBAND = ML + MU + 1
        MBA = MIN(MBAND,N)
        MEBAND = MBAND + ML
        MEB1 = MEBAND - 1
        SRUR = WM(1)
        FAC = DVNRMS(N,SAVF,EWT)
        R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
        if (R0  ==  0.0D0) R0 = 1.0D0
        DO 210 J = 1, MBA
           DO 180 I = J, N, MBAND
              YI = Y(I)
              R = MAX(SRUR*ABS(YI),R0*EWT(I))
              Y(I) = Y(I) + R
  180          CONTINUE
           call DF(TN,Y,FTEM,RPAR,IPAR)
           DO 200 JJ = J, N, MBAND
              Y(JJ) = YH(JJ,1)
              YJJ = Y(JJ)
              R = MAX(SRUR*ABS(YJJ),R0*EWT(JJ))
              FAC = -HL0/R
              I1 = MAX(JJ-MU,1)
              I2 = MIN(JJ+ML,N)
              II = JJ*MEB1 - ML + 2
              DO 190 I = I1, I2
                 WM(II+I) = (FTEM(I) - SAVF(I))*FAC
  190             CONTINUE
  200          CONTINUE
  210       CONTINUE
        NFE = NFE + MBA
  220    CONTINUE
!        ADD IDENTITY MATRIX.
!        -------------------------------------------------
     II = MBAND + 2
     DO 230 I = 1, N
        WM(II) = WM(II) + 1.0D0
        II = II + MEBAND
  230    CONTINUE
!        DO LU DECOMPOSITION OF P.
!        --------------------------------------------
     call DGBFA(WM(3),MEBAND,N,ML,MU,IWM(21),IER)
  240 CONTINUE
  return
!     ----------------------- END OF SUBROUTINE DPJAC
!     -----------------------
end
