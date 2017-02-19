subroutine PCHCI (N, H, SLOPE, D, INCFD)
!
!! PCHCI sets interior derivatives for PCHIC.
!
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      SINGLE PRECISION (PCHCI-S, DPCHCI-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!          PCHCI:  PCHIC Initial Derivative Setter.
!
!    Called by PCHIC to set derivatives needed to determine a monotone
!    piecewise cubic Hermite interpolant to the data.
!
!    Default boundary conditions are provided which are compatible
!    with monotonicity.  If the data are only piecewise monotonic, the
!    interpolant will have an extremum at each point where monotonicity
!    switches direction.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D array.
!
!    The resulting piecewise cubic Hermite function should be identical
!    (within roundoff error) to that produced by PCHIM.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N
!        REAL  H(N), SLOPE(N), D(INCFD,N)
!
!        call  PCHCI (N, H, SLOPE, D, INCFD)
!
!   Parameters:
!
!     N -- (input) number of data points.
!           If N=2, simply does linear interpolation.
!
!     H -- (input) real array of interval lengths.
!     SLOPE -- (input) real array of data slopes.
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!     D -- (output) real array of derivative values at the data points.
!           If the data are monotonic, these values will determine a
!           a monotone cubic Hermite function.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in D.
!           This argument is provided primarily for 2-D applications.
!
!    -------
!    WARNING:  This routine does no validity-checking of arguments.
!    -------
!
!  Fortran intrinsics used:  ABS, MAX, MIN.
!
!***SEE ALSO  PCHIC
!***ROUTINES CALLED  PCHST
!***REVISION HISTORY  (YYMMDD)
!   820218  DATE WRITTEN
!   820601  Modified end conditions to be continuous functions of
!           data when monotonicity switches in next interval.
!   820602  1. Modified formulas so end conditions are less prone
!             to over/underflow problems.
!           2. Minor modification to HSUM calculation.
!   820805  Converted to SLATEC library version.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR section in prologue.  (WRB)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  PCHCI
!
!  Programming notes:
!     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.
!**End
!
!  DECLARE ARGUMENTS.
!
  INTEGER  N, INCFD
  REAL  H(*), SLOPE(*), D(INCFD,*)
!
!  DECLARE LOCAL VARIABLES.
!
  INTEGER  I, NLESS1
  REAL  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, HSUM, HSUMT3, THREE, &
        W1, W2, ZERO
  SAVE ZERO, THREE
  REAL  PCHST
!
!  INITIALIZE.
!
  DATA  ZERO /0./,  THREE /3./
!***FIRST EXECUTABLE STATEMENT  PCHCI
  NLESS1 = N - 1
  DEL1 = SLOPE(1)
!
!  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!
  if (NLESS1  >  1)  go to 10
  D(1,1) = DEL1
  D(1,N) = DEL1
  go to 5000
!
!  NORMAL CASE  (N  >=  3).
!
   10 CONTINUE
  DEL2 = SLOPE(2)
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!     SHAPE-PRESERVING.
!
  HSUM = H(1) + H(2)
  W1 = (H(1) + HSUM)/HSUM
  W2 = -H(1)/HSUM
  D(1,1) = W1*DEL1 + W2*DEL2
  if ( PCHST(D(1,1),DEL1)  <=  ZERO)  THEN
     D(1,1) = ZERO
  ELSE if ( PCHST(DEL1,DEL2)  <  ZERO)  THEN
!        NEED DO THIS CHECK ONLY if MONOTONICITY SWITCHES.
     DMAX = THREE*DEL1
     if (ABS(D(1,1))  >  ABS(DMAX))  D(1,1) = DMAX
  end if
!
!  LOOP THROUGH INTERIOR POINTS.
!
  DO 50  I = 2, NLESS1
     if (I  ==  2)  go to 40
!
     HSUM = H(I-1) + H(I)
     DEL1 = DEL2
     DEL2 = SLOPE(I)
   40    CONTINUE
!
!        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!
     D(1,I) = ZERO
     if ( PCHST(DEL1,DEL2)  <=  ZERO)  go to 50
!
!        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
     HSUMT3 = HSUM+HSUM+HSUM
     W1 = (HSUM + H(I-1))/HSUMT3
     W2 = (HSUM + H(I)  )/HSUMT3
     DMAX = MAX( ABS(DEL1), ABS(DEL2) )
     DMIN = MIN( ABS(DEL1), ABS(DEL2) )
     DRAT1 = DEL1/DMAX
     DRAT2 = DEL2/DMAX
     D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
!
   50 CONTINUE
!
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!     SHAPE-PRESERVING.
!
  W1 = -H(N-1)/HSUM
  W2 = (H(N-1) + HSUM)/HSUM
  D(1,N) = W1*DEL1 + W2*DEL2
  if ( PCHST(D(1,N),DEL2)  <=  ZERO)  THEN
     D(1,N) = ZERO
  ELSE if ( PCHST(DEL1,DEL2)  <  ZERO)  THEN
!        NEED DO THIS CHECK ONLY if MONOTONICITY SWITCHES.
     DMAX = THREE*DEL2
     if (ABS(D(1,N))  >  ABS(DMAX))  D(1,N) = DMAX
  end if
!
!  NORMAL RETURN.
!
 5000 CONTINUE
  return
end
