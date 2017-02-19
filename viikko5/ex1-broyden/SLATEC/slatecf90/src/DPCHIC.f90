subroutine DPCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, &
     IERR)
!
!! DPCHIC sets derivatives needed to determine a piecewise monotone ...
!            piecewise cubic Hermite interpolant to given data.
!            User control is available over boundary conditions and/or
!            treatment of points where monotonicity switches direction.
!
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E1A
!***TYPE      DOUBLE PRECISION (PCHIC-S, DPCHIC-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
!             PCHIP, PIECEWISE CUBIC INTERPOLATION,
!             SHAPE-PRESERVING INTERPOLATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!         DPCHIC:  Piecewise Cubic Hermite Interpolation Coefficients.
!
!     Sets derivatives needed to determine a piecewise monotone piece-
!     wise cubic interpolant to the data given in X and F satisfying the
!     boundary conditions specified by IC and VC.
!
!     The treatment of points where monotonicity switches direction is
!     controlled by argument SWITCH.
!
!     To facilitate two-dimensional applications, includes an increment
!     between successive values of the F- and D arrays.
!
!     The resulting piecewise cubic Hermite function may be evaluated
!     by DPCHFE or DPCHFD.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  IC(2), N, NWK, IERR
!        DOUBLE PRECISION  VC(2), SWITCH, X(N), F(INCFD,N), D(INCFD,N),
!                          WK(NWK)
!
!        call DPCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, IERR)
!
!   Parameters:
!
!     IC -- (input) integer array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!
!           IBEG = 0  for the default boundary condition (the same as
!                     used by DPCHIM).
!           If IBEG /= 0, then its sign indicates whether the boundary
!                     derivative is to be adjusted, if necessary, to be
!                     compatible with monotonicity:
!              IBEG > 0  if no adjustment is to be performed.
!              IBEG < 0  if the derivative is to be adjusted for
!                     monotonicity.
!
!           Allowable values for the magnitude of IBEG are:
!           IBEG = 1  if first derivative at X(1) is given in VC(1).
!           IBEG = 2  if second derivative at X(1) is given in VC(1).
!           IBEG = 3  to use the 3-point difference formula for D(1).
!                     (Reverts to the default b.c. if N < 3 .)
!           IBEG = 4  to use the 4-point difference formula for D(1).
!                     (Reverts to the default b.c. if N < 4 .)
!           IBEG = 5  to set D(1) so that the second derivative is con-
!              tinuous at X(2). (Reverts to the default b.c. if N < 4.)
!              This option is somewhat analogous to the "not a knot"
!              boundary condition provided by DPCHSP.
!
!          NOTES (IBEG):
!           1. An error return is taken if ABS(IBEG) > 5 .
!           2. Only in case  IBEG <= 0  is it guaranteed that the
!              interpolant will be monotonic in the first interval.
!              If the returned value of D(1) lies between zero and
!              3*SLOPE(1), the interpolant will be monotonic.  This
!              is **NOT** checked if IBEG > 0 .
!           3. If IBEG < 0 and D(1) had to be changed to achieve mono-
!              tonicity, a warning error is returned.
!
!           IEND may take on the same values as IBEG, but applied to
!           derivative at X(N).  In case IEND = 1 or 2, the value is
!           given in VC(2).
!
!          NOTES (IEND):
!           1. An error return is taken if ABS(IEND) > 5 .
!           2. Only in case  IEND <= 0  is it guaranteed that the
!              interpolant will be monotonic in the last interval.
!              If the returned value of D(1+(N-1)*INCFD) lies between
!              zero and 3*SLOPE(N-1), the interpolant will be monotonic.
!              This is **NOT** checked if IEND > 0 .
!           3. If IEND < 0 and D(1+(N-1)*INCFD) had to be changed to
!              achieve monotonicity, a warning error is returned.
!
!     VC -- (input) real*8 array of length 2 specifying desired boundary
!           values, as indicated above.
!           VC(1) need be set only if IC(1) = 1 or 2 .
!           VC(2) need be set only if IC(2) = 1 or 2 .
!
!     SWITCH -- (input) indicates desired treatment of points where
!           direction of monotonicity switches:
!           Set SWITCH to zero if interpolant is required to be mono-
!           tonic in each interval, regardless of monotonicity of data.
!             NOTES:
!              1. This will cause D to be set to zero at all switch
!                 points, thus forcing extrema there.
!              2. The result of using this option with the default boun-
!                 dary conditions will be identical to using DPCHIM, but
!                 will generally cost more compute time.
!                 This option is provided only to facilitate comparison
!                 of different switch and/or boundary conditions.
!           Set SWITCH nonzero to use a formula based on the 3-point
!              difference formula in the vicinity of switch points.
!           If SWITCH is positive, the interpolant on each interval
!              containing an extremum is controlled to not deviate from
!              the data by more than SWITCH*DFLOC, where DFLOC is the
!              maximum of the change of F on this interval and its two
!              immediate neighbors.
!           If SWITCH is negative, no such control is to be imposed.
!
!     N -- (input) number of data points.  (Error return if N < 2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1)  <  X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of dependent variable values to be
!           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!           X(I).
!
!     D -- (output) real*8 array of derivative values at the data
!           points.  These values will determine a monotone cubic
!           Hermite function on each subinterval on which the data
!           are monotonic, except possibly adjacent to switches in
!           monotonicity. The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD < 1 .)
!
!     WK -- (scratch) real*8 array of working storage.  The user may
!           wish to know that the returned values are:
!              WK(I)     = H(I)     = X(I+1) - X(I) ;
!              WK(N-1+I) = SLOPE(I) = (F(1,I+1) - F(1,I)) / H(I)
!           for  I = 1(1)N-1.
!
!     NWK -- (input) length of work array.
!           (Error return if  NWK < 2*(N-1) .)
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
!           "Recoverable" errors:
!              IERR = -1  if N < 2 .
!              IERR = -2  if INCFD < 1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if ABS(IBEG) > 5 .
!              IERR = -5  if ABS(IEND) > 5 .
!              IERR = -6  if both of the above are true.
!              IERR = -7  if NWK < 2*(N-1) .
!             (The D array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  1. F. N. Fritsch, Piecewise Cubic Hermite Interpolation
!                 Package, Report UCRL-87285, Lawrence Livermore Natio-
!                 nal Laboratory, July 1982.  [Poster presented at the
!                 SIAM 30th Anniversary Meeting, 19-23 July 1982.]
!               2. F. N. Fritsch and J. Butland, A method for construc-
!                 ting local monotone piecewise cubic interpolants, SIAM
!                 Journal on Scientific and Statistical Computing 5, 2
!                 (June 1984), pp. 300-304.
!               3. F. N. Fritsch and R. E. Carlson, Monotone piecewise
!                 cubic interpolation, SIAM Journal on Numerical Ana-
!                 lysis 17, 2 (April 1980), pp. 238-246.
!***ROUTINES CALLED  DPCHCE, DPCHCI, DPCHCS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820218  DATE WRITTEN
!   820804  Converted to SLATEC library version.
!   870707  Corrected XERROR calls for d.p. name(s).
!   870813  Updated Reference 2.
!   890206  Corrected XERROR calls.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890703  Corrected category record.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920429  Revised format and order of references.  (WRB,FNF)
!***END PROLOGUE  DPCHIC
!  Programming notes:
!
!     To produce a single precision version, simply:
!        a. Change DPCHIC to PCHIC wherever it occurs,
!        b. Change DPCHCE to PCHCE wherever it occurs,
!        c. Change DPCHCI to PCHCI wherever it occurs,
!        d. Change DPCHCS to PCHCS wherever it occurs,
!        e. Change the double precision declarations to real, and
!        f. Change the constant  ZERO  to single precision.
!
!  DECLARE ARGUMENTS.
!
  INTEGER  IC(2), N, INCFD, NWK, IERR
  DOUBLE PRECISION  VC(2), SWITCH, X(*), F(INCFD,*), D(INCFD,*), &
   WK(NWK)
!
!  DECLARE LOCAL VARIABLES.
!
  INTEGER  I, IBEG, IEND, NLESS1
  DOUBLE PRECISION  ZERO
  SAVE ZERO
  DATA  ZERO /0.D0/
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  DPCHIC
  if ( N < 2 )  go to 5001
  if ( INCFD < 1 )  go to 5002
  DO 1  I = 2, N
     if ( X(I) <= X(I-1) )  go to 5003
    1 CONTINUE
!
  IBEG = IC(1)
  IEND = IC(2)
  IERR = 0
  if (ABS(IBEG)  >  5)  IERR = IERR - 1
  if (ABS(IEND)  >  5)  IERR = IERR - 2
  if (IERR  <  0)  go to 5004
!
!  FUNCTION DEFINITION IS OK -- GO ON.
!
  NLESS1 = N - 1
  if ( NWK  <  2*NLESS1 )  go to 5007
!
!  SET UP H AND SLOPE ARRAYS.
!
  DO 20  I = 1, NLESS1
     WK(I) = X(I+1) - X(I)
     WK(NLESS1+I) = (F(1,I+1) - F(1,I)) / WK(I)
   20 CONTINUE
!
!  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!
  if (NLESS1  >  1)  go to 1000
  D(1,1) = WK(2)
  D(1,N) = WK(2)
  go to 3000
!
!  NORMAL CASE  (N  >=  3) .
!
 1000 CONTINUE
!
!  SET INTERIOR DERIVATIVES AND DEFAULT END CONDITIONS.
!
!     --------------------------------------
  call DPCHCI (N, WK(1), WK(N), D, INCFD)
!     --------------------------------------
!
!  SET DERIVATIVES AT POINTS WHERE MONOTONICITY SWITCHES DIRECTION.
!
  if (SWITCH  ==  ZERO)  go to 3000
!     ----------------------------------------------------
  call DPCHCS (SWITCH, N, WK(1), WK(N), D, INCFD, IERR)
!     ----------------------------------------------------
  if (IERR  /=  0)  go to 5008
!
!  SET END CONDITIONS.
!
 3000 CONTINUE
  if ( (IBEG == 0) .AND. (IEND == 0) )  go to 5000
!     -------------------------------------------------------
  call DPCHCE (IC, VC, N, X, WK(1), WK(N), D, INCFD, IERR)
!     -------------------------------------------------------
  if (IERR  <  0)  go to 5009
!
!  NORMAL RETURN.
!
 5000 CONTINUE
  return
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N < 2 RETURN.
  IERR = -1
  call XERMSG ('SLATEC', 'DPCHIC', &
     'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  return
!
 5002 CONTINUE
!     INCFD < 1 RETURN.
  IERR = -2
  call XERMSG ('SLATEC', 'DPCHIC', 'INCREMENT LESS THAN ONE', IERR, &
     1)
  return
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  call XERMSG ('SLATEC', 'DPCHIC', &
     'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  return
!
 5004 CONTINUE
!     IC OUT OF RANGE RETURN.
  IERR = IERR - 3
  call XERMSG ('SLATEC', 'DPCHIC', 'IC OUT OF RANGE', IERR, 1)
  return
!
 5007 CONTINUE
!     NWK  <  2*(N-1)  return.
  IERR = -7
  call XERMSG ('SLATEC', 'DPCHIC', 'WORK ARRAY TOO SMALL', IERR, 1)
  return
!
 5008 CONTINUE
!     ERROR RETURN FROM DPCHCS.
  IERR = -8
  call XERMSG ('SLATEC', 'DPCHIC', 'ERROR RETURN FROM DPCHCS', &
     IERR, 1)
  return
!
 5009 CONTINUE
!     ERROR RETURN FROM DPCHCE.
!   *** THIS CASE SHOULD NEVER OCCUR ***
  IERR = -9
  call XERMSG ('SLATEC', 'DPCHIC', 'ERROR RETURN FROM DPCHCE', &
     IERR, 1)
  return
!------------- LAST LINE OF DPCHIC FOLLOWS -----------------------------
end
