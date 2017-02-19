subroutine PCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
!
!! PCHCS adjusts derivative values for PCHIC.
!
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      SINGLE PRECISION (PCHCS-S, DPCHCS-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!         PCHCS:  PCHIC Monotonicity Switch Derivative Setter.
!
!     Called by  PCHIC  to adjust the values of D in the vicinity of a
!     switch in direction of monotonicity, to produce a more "visually
!     pleasing" curve than that given by  PCHIM .
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IERR
!        REAL  SWITCH, H(N), SLOPE(N), D(INCFD,N)
!
!        call  PCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
!
!   Parameters:
!
!     SWITCH -- (input) indicates the amount of control desired over
!           local excursions from data.
!
!     N -- (input) number of data points.  (assumes N > 2 .)
!
!     H -- (input) real array of interval lengths.
!     SLOPE -- (input) real array of data slopes.
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!     D -- (input) real array of derivative values at the data points,
!           as determined by PCHCI.
!          (output) derivatives in the vicinity of switches in direction
!           of monotonicity may be adjusted to produce a more "visually
!           pleasing" curve.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in D.
!           This argument is provided primarily for 2-D applications.
!
!     IERR -- (output) error flag.  should be zero.
!           If negative, trouble in PCHSW.  (should never happen.)
!
!    -------
!    WARNING:  This routine does no validity-checking of arguments.
!    -------
!
!  Fortran intrinsics used:  ABS, MAX, MIN.
!
!***SEE ALSO  PCHIC
!***ROUTINES CALLED  PCHST, PCHSW
!***REVISION HISTORY  (YYMMDD)
!   820218  DATE WRITTEN
!   820617  Redesigned to (1) fix  problem with lack of continuity
!           approaching a flat-topped peak (2) be cleaner and
!           easier to verify.
!           Eliminated subroutines PCHSA and PCHSX in the process.
!   820622  1. Limited fact to not exceed one, so computed D is a
!             convex combination of PCHCI value and PCHSD value.
!           2. Changed fudge from 1 to 4 (based on experiments).
!   820623  Moved PCHSD to an inline function (eliminating MSWTYP).
!   820805  Converted to SLATEC library version.
!   870813  Minor cosmetic changes.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR section in prologue.  (WRB)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  PCHCS
!
!  Programming notes:
!     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.
!**End
!
!  DECLARE ARGUMENTS.
!
  INTEGER  N, INCFD, IERR
  REAL  SWITCH, H(*), SLOPE(*), D(INCFD,*)
!
!  DECLARE LOCAL VARIABLES.
!
  INTEGER  I, INDX, K, NLESS1
  REAL  DEL(3), DEXT, DFLOC, DFMX, FACT, FUDGE, ONE, SLMAX, &
        WTAVE(2), ZERO
  SAVE ZERO, ONE, FUDGE
  REAL  PCHST
!
!  DEFINE INLINE FUNCTION FOR WEIGHTED AVERAGE OF SLOPES.
!
  REAL  PCHSD, S1, S2, H1, H2
  PCHSD(S1,S2,H1,H2) = (H2/(H1+H2))*S1 + (H1/(H1+H2))*S2
!
!  INITIALIZE.
!
  DATA  ZERO /0./,  ONE /1./
  DATA  FUDGE /4./
!***FIRST EXECUTABLE STATEMENT  PCHCS
  IERR = 0
  NLESS1 = N - 1
!
!  LOOP OVER SEGMENTS.
!
  DO 900  I = 2, NLESS1
     if ( PCHST(SLOPE(I-1),SLOPE(I)) )  100, 300, 900
!             --------------------------
!
  100    CONTINUE
!
!....... SLOPE SWITCHES MONOTONICITY AT I-TH POINT .....................
!
!           DO NOT CHANGE D if 'UP-DOWN-UP'.
        if (I  >  2)  THEN
           if ( PCHST(SLOPE(I-2),SLOPE(I))  >  ZERO)  go to 900
!                   --------------------------
        ENDIF
        if (I  <  NLESS1)  THEN
           if ( PCHST(SLOPE(I+1),SLOPE(I-1))  >  ZERO)  go to 900
!                   ----------------------------
        ENDIF
!
!   ....... COMPUTE PROVISIONAL VALUE FOR D(1,I).
!
        DEXT = PCHSD (SLOPE(I-1), SLOPE(I), H(I-1), H(I))
!
!   ....... DETERMINE WHICH INTERVAL CONTAINS THE EXTREMUM.
!
        if ( PCHST(DEXT, SLOPE(I-1)) )  200, 900, 250
!                -----------------------
!
  200       CONTINUE
!              DEXT AND SLOPE(I-1) HAVE OPPOSITE SIGNS --
!                        EXTREMUM IS IN (X(I-1),X(I)).
           K = I-1
!              SET UP TO COMPUTE NEW VALUES FOR D(1,I-1) AND D(1,I).
           WTAVE(2) = DEXT
           if (K  >  1) &
              WTAVE(1) = PCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
           go to 400
!
  250       CONTINUE
!              DEXT AND SLOPE(I) HAVE OPPOSITE SIGNS --
!                        EXTREMUM IS IN (X(I),X(I+1)).
           K = I
!              SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
           WTAVE(1) = DEXT
           if (K  <  NLESS1) &
              WTAVE(2) = PCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))
           go to 400
!
  300    CONTINUE
!
!....... AT LEAST ONE OF SLOPE(I-1) AND SLOPE(I) IS ZERO --
!                     CHECK FOR FLAT-TOPPED PEAK .......................
!
        if (I  ==  NLESS1)  go to 900
        if ( PCHST(SLOPE(I-1), SLOPE(I+1))  >=  ZERO)  go to 900
!                -----------------------------
!
!           WE HAVE FLAT-TOPPED PEAK ON (X(I),X(I+1)).
        K = I
!           SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
        WTAVE(1) = PCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
        WTAVE(2) = PCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))
!
  400    CONTINUE
!
!....... AT THIS POINT WE HAVE DETERMINED THAT THERE WILL BE AN EXTREMUM
!        ON (X(K),X(K+1)), WHERE K=I OR I-1, AND HAVE SET ARRAY WTAVE--
!           WTAVE(1) IS A WEIGHTED AVERAGE OF SLOPE(K-1) AND SLOPE(K),
!                    if K > 1
!           WTAVE(2) IS A WEIGHTED AVERAGE OF SLOPE(K) AND SLOPE(K+1),
!                    if K < N-1
!
     SLMAX = ABS(SLOPE(K))
     if (K  >  1)    SLMAX = MAX( SLMAX, ABS(SLOPE(K-1)) )
     if (K < NLESS1) SLMAX = MAX( SLMAX, ABS(SLOPE(K+1)) )
!
     if (K  >  1)  DEL(1) = SLOPE(K-1) / SLMAX
     DEL(2) = SLOPE(K) / SLMAX
     if (K < NLESS1)  DEL(3) = SLOPE(K+1) / SLMAX
!
     if ((K > 1) .AND. (K < NLESS1))  THEN
!           NORMAL CASE -- EXTREMUM IS NOT IN A BOUNDARY INTERVAL.
        FACT = FUDGE* ABS(DEL(3)*(DEL(1)-DEL(2))*(WTAVE(2)/SLMAX))
        D(1,K) = D(1,K) + MIN(FACT,ONE)*(WTAVE(1) - D(1,K))
        FACT = FUDGE* ABS(DEL(1)*(DEL(3)-DEL(2))*(WTAVE(1)/SLMAX))
        D(1,K+1) = D(1,K+1) + MIN(FACT,ONE)*(WTAVE(2) - D(1,K+1))
     ELSE
!           SPECIAL CASE K=1 (WHICH CAN OCCUR ONLY if I=2) OR
!                        K=NLESS1 (WHICH CAN OCCUR ONLY if I=NLESS1).
        FACT = FUDGE* ABS(DEL(2))
        D(1,I) = MIN(FACT,ONE) * WTAVE(I-K+1)
!              NOTE THAT I-K+1 = 1 if K=I  (=NLESS1),
!                        I-K+1 = 2 if K=I-1(=1).
     ENDIF
!
!
!....... ADJUST if NECESSARY TO LIMIT EXCURSIONS FROM DATA.
!
     if (SWITCH  <=  ZERO)  go to 900
!
     DFLOC = H(K)*ABS(SLOPE(K))
     if (K  >  1)    DFLOC = MAX( DFLOC, H(K-1)*ABS(SLOPE(K-1)) )
     if (K < NLESS1) DFLOC = MAX( DFLOC, H(K+1)*ABS(SLOPE(K+1)) )
     DFMX = SWITCH*DFLOC
     INDX = I-K+1
!        INDX = 1 if K=I, 2 IF K=I-1.
!        ---------------------------------------------------------------
     call PCHSW (DFMX, INDX, D(1,K), D(1,K+1), H(K), SLOPE(K), IERR)
!        ---------------------------------------------------------------
     if (IERR  /=  0)  return
!
!....... END OF SEGMENT LOOP.
!
  900 CONTINUE
!
  return
end
