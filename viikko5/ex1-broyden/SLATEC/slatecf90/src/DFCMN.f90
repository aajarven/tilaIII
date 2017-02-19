subroutine DFCMN (NDATA, XDATA, YDATA, SDDATA, NORD, NBKPT, &
     BKPTIN, NCONST, XCONST, YCONST, NDERIV, MODE, COEFF, BF, XTEMP, &
     PTEMP, BKPT, G, MDG, W, MDW, WORK, IWORK)
!
!! DFCMN is subsidiary to FC.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (FCMN-S, DFCMN-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This is a companion subprogram to DFC( ).
!     The documentation for DFC( ) has complete usage instructions.
!
!***SEE ALSO  DFC
!***ROUTINES CALLED  DAXPY, DBNDAC, DBNDSL, DCOPY, DFSPVD, DFSPVN,
!                    DLSEI, DSCAL, DSORT, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   780801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890618  Completely restructured and extensively revised (WRB & RWC)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!   900604  DP version created from SP version.  (RWC)
!***END PROLOGUE  DFCMN
  INTEGER IWORK(*), MDG, MDW, MODE, NBKPT, NCONST, NDATA, NDERIV(*), &
     NORD
  DOUBLE PRECISION BF(NORD,*), BKPT(*), BKPTIN(*), COEFF(*), &
     G(MDG,*), PTEMP(*), SDDATA(*), W(MDW,*), WORK(*), &
     XCONST(*), XDATA(*), XTEMP(*), YCONST(*), YDATA(*)
!
  EXTERNAL DAXPY, DBNDAC, DBNDSL, DCOPY, DFSPVD, DFSPVN, DLSEI, &
     DSCAL, DSORT, XERMSG
!
  DOUBLE PRECISION DUMMY, PRGOPT(10), RNORM, RNORME, RNORML, XMAX, &
     XMIN, XVAL, YVAL
  INTEGER I, IDATA, IDERIV, ILEFT, INTRVL, INTW1, IP, IR, IROW, &
     ITYPE, IW1, IW2, L, LW, MT, N, NB, NEQCON, NINCON, NORDM1, &
     NORDP1, NP1
  LOGICAL BAND, NEW, VAR
  CHARACTER*8 XERN1
!
!***FIRST EXECUTABLE STATEMENT  DFCMN
!
!     Analyze input.
!
  if (NORD < 1 .OR. NORD > 20) THEN
     call XERMSG ('SLATEC', 'DFCMN', &
        'IN DFC, THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.', &
        2, 1)
     MODE = -1
     return
!
  ELSEIF (NBKPT < 2*NORD) THEN
     call XERMSG ('SLATEC', 'DFCMN', &
        'IN DFC, THE NUMBER OF KNOTS MUST BE AT LEAST TWICE ' // &
        'THE B-SPLINE ORDER.', 2, 1)
     MODE = -1
     return
  end if
!
  if (NDATA < 0) THEN
     call XERMSG ('SLATEC', 'DFCMN', &
        'IN DFC, THE NUMBER OF DATA POINTS MUST BE NONNEGATIVE.', &
        2, 1)
     MODE = -1
     return
  end if
!
!     Amount of storage allocated for W(*), IW(*).
!
  IW1 = IWORK(1)
  IW2 = IWORK(2)
  NB = (NBKPT-NORD+3)*(NORD+1) + 2*MAX(NDATA,NBKPT) + NBKPT + &
       NORD**2
!
!     See if sufficient storage has been allocated.
!
  if (IW1 < NB) THEN
     WRITE (XERN1, '(I8)') NB
     call XERMSG ('SLATEC', 'DFCMN', &
        'IN DFC, INSUFFICIENT STORAGE FOR W(*).  CHECK NB = ' // &
        XERN1, 2, 1)
     MODE = -1
     return
  end if
!
  if (MODE == 1) THEN
     BAND = .TRUE.
     VAR = .FALSE.
     NEW = .TRUE.
  ELSEIF (MODE == 2) THEN
     BAND = .FALSE.
     VAR = .TRUE.
     NEW = .TRUE.
  ELSEIF (MODE == 3) THEN
     BAND = .TRUE.
     VAR = .FALSE.
     NEW = .FALSE.
  ELSEIF (MODE == 4) THEN
     BAND = .FALSE.
     VAR = .TRUE.
     NEW = .FALSE.
  ELSE
     call XERMSG ('SLATEC', 'DFCMN', &
        'IN DFC, INPUT VALUE OF MODE MUST BE 1-4.', 2, 1)
     MODE = -1
     return
  end if
  MODE = 0
!
!     Sort the breakpoints.
!
  call DCOPY (NBKPT, BKPTIN, 1, BKPT, 1)
  call DSORT (BKPT, DUMMY, NBKPT, 1)
!
!     Initialize variables.
!
  NEQCON = 0
  NINCON = 0
  DO 100 I = 1,NCONST
     L = NDERIV(I)
     ITYPE = MOD(L,4)
     if (ITYPE < 2) THEN
        NINCON = NINCON + 1
     ELSE
        NEQCON = NEQCON + 1
     ENDIF
  100 CONTINUE
!
!     Compute the number of variables.
!
  N = NBKPT - NORD
  NP1 = N + 1
  LW = NB + (NP1+NCONST)*NP1 + 2*(NEQCON+NP1) + (NINCON+NP1) + &
       (NINCON+2)*(NP1+6)
  INTW1 = NINCON + 2*NP1
!
!     Save interval containing knots.
!
  XMIN = BKPT(NORD)
  XMAX = BKPT(NP1)
!
!     Find the smallest referenced independent variable value in any
!     constraint.
!
  DO 110 I = 1,NCONST
     XMIN = MIN(XMIN,XCONST(I))
     XMAX = MAX(XMAX,XCONST(I))
  110 CONTINUE
  NORDM1 = NORD - 1
  NORDP1 = NORD + 1
!
!     Define the option vector PRGOPT(1-10) for use in DLSEI( ).
!
  PRGOPT(1) = 4
!
!     Set the covariance matrix computation flag.
!
  PRGOPT(2) = 1
  if (VAR) THEN
     PRGOPT(3) = 1
  ELSE
     PRGOPT(3) = 0
  end if
!
!     Increase the rank determination tolerances for both equality
!     constraint equations and least squares equations.
!
  PRGOPT(4) = 7
  PRGOPT(5) = 4
  PRGOPT(6) = 1.D-4
!
  PRGOPT(7) = 10
  PRGOPT(8) = 5
  PRGOPT(9) = 1.D-4
!
  PRGOPT(10) = 1
!
!     Turn off work array length checking in DLSEI( ).
!
  IWORK(1) = 0
  IWORK(2) = 0
!
!     Initialize variables and analyze input.
!
  if (NEW) THEN
!
!        To process least squares equations sort data and an array of
!        pointers.
!
     call DCOPY (NDATA, XDATA, 1, XTEMP, 1)
     DO 120 I = 1,NDATA
        PTEMP(I) = I
  120    CONTINUE
!
     if (NDATA > 0) THEN
        call DSORT (XTEMP, PTEMP, NDATA, 2)
        XMIN = MIN(XMIN,XTEMP(1))
        XMAX = MAX(XMAX,XTEMP(NDATA))
     ENDIF
!
!        Fix breakpoint array if needed.
!
     DO 130 I = 1,NORD
        BKPT(I) = MIN(BKPT(I),XMIN)
  130    CONTINUE
!
     DO 140 I = NP1,NBKPT
        BKPT(I) = MAX(BKPT(I),XMAX)
  140    CONTINUE
!
!        Initialize parameters of banded matrix processor, DBNDAC( ).
!
     MT = 0
     IP = 1
     IR = 1
     ILEFT = NORD
     DO 160 IDATA = 1,NDATA
!
!           Sorted indices are in PTEMP(*).
!
        L = PTEMP(IDATA)
        XVAL = XDATA(L)
!
!           When interval changes, process equations in the last block.
!
        if (XVAL >= BKPT(ILEFT+1)) THEN
           call DBNDAC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
           MT = 0
!
!              Move pointer up to have BKPT(ILEFT) <= XVAL,
!                 ILEFT < NP1.
!
  150          if (XVAL >= BKPT(ILEFT+1) .AND. ILEFT < N) THEN
              ILEFT = ILEFT + 1
              go to 150
           ENDIF
        ENDIF
!
!           Obtain B-spline function value.
!
        call DFSPVN (BKPT, NORD, 1, XVAL, ILEFT, BF)
!
!           Move row into place.
!
        IROW = IR + MT
        MT = MT + 1
        call DCOPY (NORD, BF, 1, G(IROW,1), MDG)
        G(IROW,NORDP1) = YDATA(L)
!
!           Scale data if uncertainty is nonzero.
!
        if (SDDATA(L) /= 0.D0) call DSCAL (NORDP1, 1.D0/SDDATA(L), &
                                    G(IROW,1), MDG)
!
!           When staging work area is exhausted, process rows.
!
        if (IROW == MDG-1) THEN
           call DBNDAC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
           MT = 0
        ENDIF
  160    CONTINUE
!
!        Process last block of equations.
!
     call DBNDAC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
!
!        Last call to adjust block positioning.
!
     call dinit ( NORDP1, 0.D0, G(IR,1), MDG)
     call DBNDAC (G, MDG, NORD, IP, IR, 1, NP1)
  end if
!
  BAND = BAND .AND. NCONST == 0
  DO 170 I = 1,N
     BAND = BAND .AND. G(I,1) /= 0.D0
  170 CONTINUE
!
!     Process banded least squares equations.
!
  if (BAND) THEN
     call DBNDSL (1, G, MDG, NORD, IP, IR, COEFF, N, RNORM)
     return
  end if
!
!     Check further for sufficient storage in working arrays.
!
  if (IW1 < LW) THEN
     WRITE (XERN1, '(I8)') LW
     call XERMSG ('SLATEC', 'DFCMN', &
        'IN DFC, INSUFFICIENT STORAGE FOR W(*).  CHECK LW = ' // &
        XERN1, 2, 1)
     MODE = -1
     return
  end if
!
  if (IW2 < INTW1) THEN
     WRITE (XERN1, '(I8)') INTW1
     call XERMSG ('SLATEC', 'DFCMN', &
        'IN DFC, INSUFFICIENT STORAGE FOR IW(*).  CHECK IW1 = ' // &
        XERN1, 2, 1)
     MODE = -1
     return
  end if
!
!     Write equality constraints.
!     Analyze constraint indicators for an equality constraint.
!
  NEQCON = 0
  DO 220 IDATA = 1,NCONST
     L = NDERIV(IDATA)
     ITYPE = MOD(L,4)
     if (ITYPE > 1) THEN
        IDERIV = L/4
        NEQCON = NEQCON + 1
        ILEFT = NORD
        XVAL = XCONST(IDATA)
!
  180       if (XVAL < BKPT(ILEFT+1) .OR. ILEFT >= N) go to 190
        ILEFT = ILEFT + 1
        go to 180
!
  190       call DFSPVD (BKPT, NORD, XVAL, ILEFT, BF, IDERIV+1)
        call dinit (NP1, 0.D0, W(NEQCON,1), MDW)
        call DCOPY (NORD, BF(1,IDERIV+1), 1, W(NEQCON,ILEFT-NORDM1), &
                    MDW)
!
        if (ITYPE == 2) THEN
           W(NEQCON,NP1) = YCONST(IDATA)
        ELSE
           ILEFT = NORD
           YVAL = YCONST(IDATA)
!
  200          if (YVAL < BKPT(ILEFT+1) .OR. ILEFT >= N) go to 210
           ILEFT = ILEFT + 1
           go to 200
!
  210          call DFSPVD (BKPT, NORD, YVAL, ILEFT, BF, IDERIV+1)
           call DAXPY (NORD, -1.D0, BF(1, IDERIV+1), 1, &
                       W(NEQCON, ILEFT-NORDM1), MDW)
        ENDIF
     ENDIF
  220 CONTINUE
!
!     Transfer least squares data.
!
  DO 230 I = 1,NP1
     IROW = I + NEQCON
     call dinit ( N, 0.D0, W(IROW,1), MDW)
     call DCOPY (MIN(NP1-I, NORD), G(I,1), MDG, W(IROW,I), MDW)
     W(IROW,NP1) = G(I,NORDP1)
  230 CONTINUE
!
!     Write inequality constraints.
!     Analyze constraint indicators for inequality constraints.
!
  NINCON = 0
  DO 260 IDATA = 1,NCONST
     L = NDERIV(IDATA)
     ITYPE = MOD(L,4)
     if (ITYPE < 2) THEN
        IDERIV = L/4
        NINCON = NINCON + 1
        ILEFT = NORD
        XVAL = XCONST(IDATA)
!
  240       if (XVAL < BKPT(ILEFT+1) .OR. ILEFT >= N) go to 250
        ILEFT = ILEFT + 1
        go to 240
!
  250       call DFSPVD (BKPT, NORD, XVAL, ILEFT, BF, IDERIV+1)
        IROW = NEQCON + NP1 + NINCON
        call dinit ( N, 0.D0, W(IROW,1), MDW)
        INTRVL = ILEFT - NORDM1
        call DCOPY (NORD, BF(1, IDERIV+1), 1, W(IROW, INTRVL), MDW)
!
        if (ITYPE == 1) THEN
           W(IROW,NP1) = YCONST(IDATA)
        ELSE
           W(IROW,NP1) = -YCONST(IDATA)
           call DSCAL (NORD, -1.D0, W(IROW, INTRVL), MDW)
        ENDIF
     ENDIF
  260 CONTINUE
!
!     Solve constrained least squares equations.
!
  call DLSEI(W, MDW, NEQCON, NP1, NINCON, N, PRGOPT, COEFF, RNORME, &
            RNORML, MODE, WORK, IWORK)
  return
end
