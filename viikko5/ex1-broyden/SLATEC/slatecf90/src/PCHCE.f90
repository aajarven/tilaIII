subroutine PCHCE (IC, VC, N, X, H, SLOPE, D, INCFD, IERR)
!
!! PCHCE sets boundary conditions for PCHIC.
!
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      SINGLE PRECISION (PCHCE-S, DPCHCE-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!          PCHCE:  PCHIC End Derivative Setter.
!
!    Called by PCHIC to set end derivatives as requested by the user.
!    It must be called after interior derivative values have been set.
!                      -----
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D array.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  IC(2), N, IERR
!        REAL  VC(2), X(N), H(N), SLOPE(N), D(INCFD,N)
!
!        call  PCHCE (IC, VC, N, X, H, SLOPE, D, INCFD, IERR)
!
!   Parameters:
!
!     IC -- (input) integer array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!           ( see prologue to PCHIC for details. )
!
!     VC -- (input) real array of length 2 specifying desired boundary
!           values.  VC(1) need be set only if IC(1) = 2 or 3 .
!                    VC(2) need be set only if IC(2) = 2 or 3 .
!
!     N -- (input) number of data points.  (assumes N >= 2)
!
!     X -- (input) real array of independent variable values.  (the
!           elements of X are assumed to be strictly increasing.)
!
!     H -- (input) real array of interval lengths.
!     SLOPE -- (input) real array of data slopes.
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!     D -- (input) real array of derivative values at the data points.
!           The value corresponding to X(I) must be stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!          (output) the value of D at X(1) and/or X(N) is changed, if
!           necessary, to produce the requested boundary conditions.
!           no other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in D.
!           This argument is provided primarily for 2-D applications.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if IBEG < 0 and D(1) had to be adjusted for
!                        monotonicity.
!              IERR = 2  if IEND < 0 and D(1+(N-1)*INCFD) had to be
!                        adjusted for monotonicity.
!              IERR = 3  if both of the above are true.
!
!    -------
!    WARNING:  This routine does no validity-checking of arguments.
!    -------
!
!  Fortran intrinsics used:  ABS.
!
!***SEE ALSO  PCHIC
!***ROUTINES CALLED  PCHDF, PCHST, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820218  DATE WRITTEN
!   820805  Converted to SLATEC library version.
!   870707  Minor corrections made to prologue..
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR section in prologue.  (WRB)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  PCHCE
!
!  Programming notes:
!     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.
!     2. One could reduce the number of arguments and amount of local
!        storage, at the expense of reduced code clarity, by passing in
!        the array WK (rather than splitting it into H and SLOPE) and
!        increasing its length enough to incorporate STEMP and XTEMP.
!     3. The two monotonicity checks only use the sufficient conditions.
!        Thus, it is possible (but unlikely) for a boundary condition to
!        be changed, even though the original interpolant was monotonic.
!        (At least the result is a continuous function of the data.)
!**End
!
!  DECLARE ARGUMENTS.
!
  INTEGER  IC(2), N, INCFD, IERR
  REAL  VC(2), X(*), H(*), SLOPE(*), D(INCFD,*)
!
!  DECLARE LOCAL VARIABLES.
!
  INTEGER  IBEG, IEND, IERF, INDEX, J, K
  REAL  HALF, STEMP(3), THREE, TWO, XTEMP(4), ZERO
  SAVE ZERO, HALF, TWO, THREE
  REAL  PCHDF, PCHST
!
!  INITIALIZE.
!
  DATA  ZERO /0./,  HALF /0.5/,  TWO /2./,  THREE /3./
!
!***FIRST EXECUTABLE STATEMENT  PCHCE
  IBEG = IC(1)
  IEND = IC(2)
  IERR = 0
!
!  SET TO DEFAULT BOUNDARY CONDITIONS if N IS TOO SMALL.
!
  if ( ABS(IBEG) > N )  IBEG = 0
  if ( ABS(IEND) > N )  IEND = 0
!
!  TREAT BEGINNING BOUNDARY CONDITION.
!
  if (IBEG  ==  0)  go to 2000
  K = ABS(IBEG)
  if (K  ==  1)  THEN
!        BOUNDARY VALUE PROVIDED.
     D(1,1) = VC(1)
  ELSE if (K  ==  2)  THEN
!        BOUNDARY SECOND DERIVATIVE PROVIDED.
     D(1,1) = HALF*( (THREE*SLOPE(1) - D(1,2)) - HALF*VC(1)*H(1) )
  ELSE if (K  <  5)  THEN
!        USE K-POINT DERIVATIVE FORMULA.
!        PICK UP FIRST K POINTS, IN REVERSE ORDER.
     DO 10  J = 1, K
        INDEX = K-J+1
!           INDEX RUNS FROM K DOWN TO 1.
        XTEMP(J) = X(INDEX)
        if (J  <  K)  STEMP(J) = SLOPE(INDEX-1)
   10    CONTINUE
!                 -----------------------------
     D(1,1) = PCHDF (K, XTEMP, STEMP, IERF)
!                 -----------------------------
     if (IERF  /=  0)  go to 5001
  ELSE
!        USE 'NOT A KNOT' CONDITION.
     D(1,1) = ( THREE*(H(1)*SLOPE(2) + H(2)*SLOPE(1)) &
               - TWO*(H(1)+H(2))*D(1,2) - H(1)*D(1,3) ) / H(2)
  end if
!
  if (IBEG  >  0)  go to 2000
!
!  CHECK D(1,1) FOR COMPATIBILITY WITH MONOTONICITY.
!
  if (SLOPE(1)  ==  ZERO)  THEN
     if (D(1,1)  /=  ZERO)  THEN
        D(1,1) = ZERO
        IERR = IERR + 1
     ENDIF
  ELSE if ( PCHST(D(1,1),SLOPE(1))  <  ZERO)  THEN
     D(1,1) = ZERO
     IERR = IERR + 1
  ELSE if ( ABS(D(1,1))  >  THREE*ABS(SLOPE(1)) )  THEN
     D(1,1) = THREE*SLOPE(1)
     IERR = IERR + 1
  end if
!
!  TREAT END BOUNDARY CONDITION.
!
 2000 CONTINUE
  if (IEND  ==  0)  go to 5000
  K = ABS(IEND)
  if (K  ==  1)  THEN
!        BOUNDARY VALUE PROVIDED.
     D(1,N) = VC(2)
  ELSE if (K  ==  2)  THEN
!        BOUNDARY SECOND DERIVATIVE PROVIDED.
     D(1,N) = HALF*( (THREE*SLOPE(N-1) - D(1,N-1)) + &
                                             HALF*VC(2)*H(N-1) )
  ELSE if (K  <  5)  THEN
!        USE K-POINT DERIVATIVE FORMULA.
!        PICK UP LAST K POINTS.
     DO 2010  J = 1, K
        INDEX = N-K+J
!           INDEX RUNS FROM N+1-K UP TO N.
        XTEMP(J) = X(INDEX)
        if (J  <  K)  STEMP(J) = SLOPE(INDEX)
 2010    CONTINUE
!                 -----------------------------
     D(1,N) = PCHDF (K, XTEMP, STEMP, IERF)
!                 -----------------------------
     if (IERF  /=  0)  go to 5001
  ELSE
!        USE 'NOT A KNOT' CONDITION.
     D(1,N) = ( THREE*(H(N-1)*SLOPE(N-2) + H(N-2)*SLOPE(N-1)) &
               - TWO*(H(N-1)+H(N-2))*D(1,N-1) - H(N-1)*D(1,N-2) ) &
                                                           / H(N-2)
  end if
!
  if (IEND  >  0)  go to 5000
!
!  CHECK D(1,N) FOR COMPATIBILITY WITH MONOTONICITY.
!
  if (SLOPE(N-1)  ==  ZERO)  THEN
     if (D(1,N)  /=  ZERO)  THEN
        D(1,N) = ZERO
        IERR = IERR + 2
     ENDIF
  ELSE if ( PCHST(D(1,N),SLOPE(N-1))  <  ZERO)  THEN
     D(1,N) = ZERO
     IERR = IERR + 2
  ELSE if ( ABS(D(1,N))  >  THREE*ABS(SLOPE(N-1)) )  THEN
     D(1,N) = THREE*SLOPE(N-1)
     IERR = IERR + 2
  end if
!
!  NORMAL RETURN.
!
 5000 CONTINUE
  return
!
!  ERROR RETURN.
!
 5001 CONTINUE
!     ERROR RETURN FROM PCHDF.
!   *** THIS CASE SHOULD NEVER OCCUR ***
  IERR = -1
  call XERMSG ('SLATEC', 'PCHCE', 'ERROR RETURN FROM PCHDF', IERR, &
     1)
  return
end
