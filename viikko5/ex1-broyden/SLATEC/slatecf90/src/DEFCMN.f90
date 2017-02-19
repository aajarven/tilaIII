subroutine DEFCMN (NDATA, XDATA, YDATA, SDDATA, NORD, NBKPT, &
     BKPTIN, MDEIN, MDEOUT, COEFF, BF, XTEMP, PTEMP, BKPT, G, MDG, &
     W, MDW, LW)
!
!! DEFCMN is subsidiary to DEFC.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (EFCMN-S, DEFCMN-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!     This is a companion subprogram to DEFC( ).
!     This subprogram does weighted least squares fitting of data by
!     B-spline curves.
!     The documentation for DEFC( ) has complete usage instructions.
!
!***SEE ALSO  DEFC
!***ROUTINES CALLED  DBNDAC, DBNDSL, DCOPY, DFSPVN, DSCAL, DSORT, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890618  Completely restructured and extensively revised (WRB & RWC)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!   900604  DP version created from SP version.  (RWC)
!***END PROLOGUE  DEFCMN
  INTEGER LW, MDEIN, MDEOUT, MDG, MDW, NBKPT, NDATA, NORD
  DOUBLE PRECISION BF(NORD,*), BKPT(*), BKPTIN(*), COEFF(*), &
     G(MDG,*), PTEMP(*), SDDATA(*), W(MDW,*), XDATA(*), XTEMP(*), &
     YDATA(*)
!
  EXTERNAL DBNDAC, DBNDSL, DCOPY, DFSPVN, DSCAL, DSORT, XERMSG
!
  DOUBLE PRECISION DUMMY, RNORM, XMAX, XMIN, XVAL
  INTEGER I, IDATA, ILEFT, INTSEQ, IP, IR, IROW, L, MT, N, NB, &
     NORDM1, NORDP1, NP1
  CHARACTER*8 XERN1, XERN2
!
!***FIRST EXECUTABLE STATEMENT  DEFCMN
!
!     Initialize variables and analyze input.
!
  N = NBKPT - NORD
  NP1 = N + 1
!
!     Initially set all output coefficients to zero.
!
  call dinit (N, 0.D0, COEFF, 1)
  MDEOUT = -1
  if (NORD < 1 .OR. NORD > 20) THEN
     call XERMSG ('SLATEC', 'DEFCMN', &
        'IN DEFC, THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.', &
        3, 1)
     return
  end if
!
  if (NBKPT < 2*NORD) THEN
     call XERMSG ('SLATEC', 'DEFCMN', &
        'IN DEFC, THE NUMBER OF KNOTS MUST BE AT LEAST TWICE ' // &
        'THE B-SPLINE ORDER.', 4, 1)
     return
  end if
!
  if (NDATA < 0) THEN
     call XERMSG ('SLATEC', 'DEFCMN', &
        'IN DEFC, THE NUMBER OF DATA POINTS MUST BE NONNEGATIVE.', &
        5, 1)
     return
  end if
!
  NB = (NBKPT-NORD+3)*(NORD+1) + (NBKPT+1)*(NORD+1) + &
       2*MAX(NBKPT,NDATA) + NBKPT + NORD**2
  if (LW  <  NB) THEN
     WRITE (XERN1, '(I8)') NB
     WRITE (XERN2, '(I8)') LW
     call XERMSG ('SLATEC', 'DEFCMN', &
        'IN DEFC, INSUFFICIENT STORAGE FOR W(*).  CHECK FORMULA ' // &
        'THAT READS LW >=  ... .  NEED = ' // XERN1 // &
        ' GIVEN = ' // XERN2, 6, 1)
     MDEOUT = -1
     return
  end if
!
  if (MDEIN /= 1 .AND. MDEIN /= 2) THEN
     call XERMSG ('SLATEC', 'DEFCMN', &
        'IN DEFC, INPUT VALUE OF MDEIN MUST BE 1-2.', 7, 1)
     return
  end if
!
!     Sort the breakpoints.
!
  call DCOPY (NBKPT, BKPTIN, 1, BKPT, 1)
  call DSORT (BKPT, DUMMY, NBKPT, 1)
!
!     Save interval containing knots.
!
  XMIN = BKPT(NORD)
  XMAX = BKPT(NP1)
  NORDM1 = NORD - 1
  NORDP1 = NORD + 1
!
!     Process least squares equations.
!
!     Sort data and an array of pointers.
!
  call DCOPY (NDATA, XDATA, 1, XTEMP, 1)
  DO 100 I = 1,NDATA
     PTEMP(I) = I
  100 CONTINUE
!
  if (NDATA > 0) THEN
     call DSORT (XTEMP, PTEMP, NDATA, 2)
     XMIN = MIN(XMIN,XTEMP(1))
     XMAX = MAX(XMAX,XTEMP(NDATA))
  end if
!
!     Fix breakpoint array if needed. This should only involve very
!     minor differences with the input array of breakpoints.
!
  DO 110 I = 1,NORD
     BKPT(I) = MIN(BKPT(I),XMIN)
  110 CONTINUE
!
  DO 120 I = NP1,NBKPT
     BKPT(I) = MAX(BKPT(I),XMAX)
  120 CONTINUE
!
!     Initialize parameters of banded matrix processor, DBNDAC( ).
!
  MT = 0
  IP = 1
  IR = 1
  ILEFT = NORD
  INTSEQ = 1
  DO 150 IDATA = 1,NDATA
!
!        Sorted indices are in PTEMP(*).
!
     L = PTEMP(IDATA)
     XVAL = XDATA(L)
!
!        When interval changes, process equations in the last block.
!
     if (XVAL >= BKPT(ILEFT+1)) THEN
        call DBNDAC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
        MT = 0
!
!           Move pointer up to have BKPT(ILEFT) <= XVAL, ILEFT <= N.
!
        DO 130 ILEFT = ILEFT,N
           if (XVAL < BKPT(ILEFT+1)) go to 140
           if (MDEIN == 2) THEN
!
!                 Data is being sequentially accumulated.
!                 Transfer previously accumulated rows from W(*,*) to
!                 G(*,*) and process them.
!
              call DCOPY (NORDP1, W(INTSEQ,1), MDW, G(IR,1), MDG)
              call DBNDAC (G, MDG, NORD, IP, IR, 1, INTSEQ)
              INTSEQ = INTSEQ + 1
           ENDIF
  130       CONTINUE
     ENDIF
!
!        Obtain B-spline function value.
!
  140    call DFSPVN (BKPT, NORD, 1, XVAL, ILEFT, BF)
!
!        Move row into place.
!
     IROW = IR + MT
     MT = MT + 1
     call DCOPY (NORD, BF, 1, G(IROW,1), MDG)
     G(IROW,NORDP1) = YDATA(L)
!
!        Scale data if uncertainty is nonzero.
!
     if (SDDATA(L) /= 0.D0) call DSCAL (NORDP1, 1.D0/SDDATA(L), &
                                 G(IROW,1), MDG)
!
!        When staging work area is exhausted, process rows.
!
     if (IROW == MDG-1) THEN
        call DBNDAC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
        MT = 0
     ENDIF
  150 CONTINUE
!
!     Process last block of equations.
!
  call DBNDAC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
!
!     Finish processing any previously accumulated rows from W(*,*)
!     to G(*,*).
!
  if (MDEIN == 2) THEN
     DO 160 I = INTSEQ,NP1
        call DCOPY (NORDP1, W(I,1), MDW, G(IR,1), MDG)
        call DBNDAC (G, MDG, NORD, IP, IR, 1, MIN(N,I))
  160    CONTINUE
  end if
!
!     Last call to adjust block positioning.
!
  call dinit ( NORDP1, 0.D0, G(IR,1), MDG)
  call DBNDAC (G, MDG, NORD, IP, IR, 1, NP1)
!
!     Transfer accumulated rows from G(*,*) to W(*,*) for
!     possible later sequential accumulation.
!
  DO 170 I = 1,NP1
     call DCOPY (NORDP1, G(I,1), MDG, W(I,1), MDW)
  170 CONTINUE
!
!     Solve for coefficients when possible.
!
  DO 180 I = 1,N
     if (G(I,1) == 0.D0) THEN
        MDEOUT = 2
        return
     ENDIF
  180 CONTINUE
!
!     All the diagonal terms in the accumulated triangular
!     matrix are nonzero.  The solution can be computed but
!     it may be unsuitable for further use due to poor
!     conditioning or the lack of constraints.  No checking
!     for either of these is done here.
!
  call DBNDSL (1, G, MDG, NORD, IP, IR, COEFF, N, RNORM)
  MDEOUT = 1
  return
end
