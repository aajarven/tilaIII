subroutine PCHCM (N, X, F, D, INCFD, SKIP, ISMON, IERR)
!
!! PCHCM checks a cubic Hermite function for monotonicity.
!
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3
!***TYPE      SINGLE PRECISION (PCHCM-S, DPCHCM-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
!             PCHIP, PIECEWISE CUBIC INTERPOLATION, UTILITY ROUTINE
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Computing & Mathematics Research Division
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
! *Usage:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, ISMON(N), IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N)
!        LOGICAL  SKIP
!
!        call  PCHCM (N, X, F, D, INCFD, SKIP, ISMON, IERR)
!
! *Arguments:
!
!     N:IN  is the number of data points.  (Error return if N < 2 .)
!
!     X:IN  is a real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1)  <  X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F:IN  is a real array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D:IN  is a real array of derivative values.  D(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     INCFD:IN  is the increment between successive values in F and D.
!           (Error return if  INCFD < 1 .)
!
!     SKIP:INOUT  is a logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed.
!           SKIP will be set to .TRUE. on normal return.
!
!     ISMON:OUT  is an integer array indicating on which intervals the
!           PCH function defined by  N, X, F, D  is monotonic.
!           For data interval [X(I),X(I+1)],
!             ISMON(I) = -3  if function is probably decreasing;
!             ISMON(I) = -1  if function is strictly decreasing;
!             ISMON(I) =  0  if function is constant;
!             ISMON(I) =  1  if function is strictly increasing;
!             ISMON(I) =  2  if function is non-monotonic;
!             ISMON(I) =  3  if function is probably increasing.
!                If ABS(ISMON)=3, this means that the D values are near
!                the boundary of the monotonicity region.  A small
!                increase produces non-monotonicity; decrease, strict
!                monotonicity.
!           The above applies to I=1(1)N-1.  ISMON(N) indicates whether
!              the entire function is monotonic on [X(1),X(N)].
!
!     IERR:OUT  is an error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N < 2 .
!              IERR = -2  if INCFD < 1 .
!              IERR = -3  if the X-array is not strictly increasing.
!          (The ISMON-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
! *Description:
!
!          PCHCM:  Piecewise Cubic Hermite -- Check Monotonicity.
!
!     Checks the piecewise cubic Hermite function defined by  N,X,F,D
!     for monotonicity.
!
!     To provide compatibility with PCHIM and PCHIC, includes an
!     increment between successive values of the F- and D arrays.
!
! *Cautions:
!     This provides the same capability as old PCHMC, except that a
!     new output value, -3, was added February 1989.  (Formerly, -3
!     and +3 were lumped together in the single value 3.)  Codes that
!     flag nonmonotonicity by "IF (ISMON == 2)" need not be changed.
!     Codes that check via "IF (ISMON >= 3)" should change the test to
!     "IF (IABS(ISMON) >= 3)".  Codes that declare monotonicity via
!     "IF (ISMON <= 1)" should change to "IF (IABS(ISMON) <= 1)".
!
!***REFERENCES  F. N. Fritsch and R. E. Carlson, Monotone piecewise
!                 cubic interpolation, SIAM Journal on Numerical Ana-
!                 lysis 17, 2 (April 1980), pp. 238-246.
!***ROUTINES CALLED  CHFCM, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820518  DATE WRITTEN
!   820804  Converted to SLATEC library version.
!   831201  Reversed order of subscripts of F and D, so that the
!           routine will work properly when INCFD > 1 .  (Bug!)
!   870707  Minor cosmetic changes to prologue.
!   890208  Added possible ISMON value of -3 and modified code so
!           that 1,3,-1 produces ISMON(N)=2, rather than 3.
!   890306  Added caution about changed output.
!   890407  Changed name from PCHMC to PCHCM, as requested at the
!           March 1989 SLATEC CML meeting, and made a few other
!           minor modifications necessitated by this change.
!   890407  Converted to new SLATEC format.
!   890407  Modified DESCRIPTION to LDOC format.
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920429  Revised format and order of references.  (WRB,FNF)
!***END PROLOGUE  PCHCM
!
!  Fortran intrinsics used:  ISIGN.
!  Other routines used:  CHFCM, XERMSG.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     An alternate organization would have separate loops for computing
!     ISMON(i), i=1,...,NSEG, and for the computation of ISMON(N).  The
!     first loop can be readily parallelized, since the NSEG calls to
!     CHFCM are independent.  The second loop can be cut short if
!     ISMON(N) is ever equal to 2, for it cannot be changed further.
!
!     To produce a double precision version, simply:
!        a. Change PCHCM to DPCHCM wherever it occurs,
!        b. Change CHFCM to DCHFCM wherever it occurs, and
!        c. Change the real declarations to double precision.
!
!  DECLARE ARGUMENTS.
!
  INTEGER  N, INCFD, ISMON(N), IERR
  REAL  X(N), F(INCFD,N), D(INCFD,N)
  LOGICAL  SKIP
!
!  DECLARE LOCAL VARIABLES.
!
  INTEGER  I, NSEG
  REAL  DELTA
  INTEGER  CHFCM
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  PCHCM
  if (SKIP)  go to 5
!
  if ( N < 2 )  go to 5001
  if ( INCFD < 1 )  go to 5002
  DO 1  I = 2, N
     if ( X(I) <= X(I-1) )  go to 5003
    1 CONTINUE
  SKIP = .TRUE.
!
!  FUNCTION DEFINITION IS OK -- GO ON.
!
    5 CONTINUE
  NSEG = N - 1
  DO 90  I = 1, NSEG
     DELTA = (F(1,I+1)-F(1,I))/(X(I+1)-X(I))
!                   -------------------------------
     ISMON(I) = CHFCM (D(1,I), D(1,I+1), DELTA)
!                   -------------------------------
     if (I  ==  1)  THEN
        ISMON(N) = ISMON(1)
     ELSE
!           Need to figure out cumulative monotonicity from following
!           "multiplication table":
!
!                    +        I S M O N (I)
!                     +  -3  -1   0   1   3   2
!                      +------------------------+
!               I   -3 I -3  -3  -3   2   2   2 I
!               S   -1 I -3  -1  -1   2   2   2 I
!               M    0 I -3  -1   0   1   3   2 I
!               O    1 I  2   2   1   1   3   2 I
!               N    3 I  2   2   3   3   3   2 I
!              (N)   2 I  2   2   2   2   2   2 I
!                      +------------------------+
!           Note that the 2 row and column are out of order so as not
!           to obscure the symmetry in the rest of the table.
!
!           No change needed if equal or constant on this interval or
!           already declared nonmonotonic.
        if ( (ISMON(I) /= ISMON(N)) .AND. (ISMON(I) /= 0) &
                                    .AND. (ISMON(N) /= 2) )  THEN
           if ( (ISMON(I) == 2) .OR. (ISMON(N) == 0) )  THEN
              ISMON(N) =  ISMON(I)
           ELSE if (ISMON(I)*ISMON(N)  <  0)  THEN
!                 This interval has opposite sense from curve so far.
              ISMON(N) = 2
           ELSE
!                 At this point, both are nonzero with same sign, and
!                 we have already eliminated case both +-1.
              ISMON(N) = ISIGN (3, ISMON(N))
           ENDIF
        ENDIF
     ENDIF
   90 CONTINUE
!
!  NORMAL RETURN.
!
  IERR = 0
  return
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N < 2 RETURN.
  IERR = -1
  call XERMSG ('SLATEC', 'PCHCM', &
     'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  return
!
 5002 CONTINUE
!     INCFD < 1 RETURN.
  IERR = -2
  call XERMSG ('SLATEC', 'PCHCM', 'INCREMENT LESS THAN ONE', IERR, &
     1)
  return
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  call XERMSG ('SLATEC', 'PCHCM', 'X-ARRAY NOT STRICTLY INCREASING' &
     , IERR, 1)
  return
end
