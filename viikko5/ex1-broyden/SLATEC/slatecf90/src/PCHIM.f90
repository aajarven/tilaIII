subroutine PCHIM (N, X, F, D, INCFD, IERR)
!
!! PCHIM sets derivatives needed to determine a monotone piecewise ...
!            cubic Hermite interpolant to given data.  Boundary values
!            are provided which are compatible with monotonicity.  The
!            interpolant will have an extremum at each point where mono-
!            tonicity switches direction.  (See PCHIC if user control is
!            desired over boundary or switch conditions.)
!
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E1A
!***TYPE      SINGLE PRECISION (PCHIM-S, DPCHIM-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
!             PCHIP, PIECEWISE CUBIC INTERPOLATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          PCHIM:  Piecewise Cubic Hermite Interpolation to
!                  Monotone data.
!
!     Sets derivatives needed to determine a monotone piecewise cubic
!     Hermite interpolant to the data given in X and F.
!
!     Default boundary conditions are provided which are compatible
!     with monotonicity.  (See PCHIC if user control of boundary con-
!     ditions is desired.)
!
!     If the data are only piecewise monotonic, the interpolant will
!     have an extremum at each point where monotonicity switches direc-
!     tion.  (See PCHIC if user control is desired in such cases.)
!
!     To facilitate two-dimensional applications, includes an increment
!     between successive values of the F- and D arrays.
!
!     The resulting piecewise cubic Hermite function may be evaluated
!     by PCHFE or PCHFD.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N)
!
!        call  PCHIM (N, X, F, D, INCFD, IERR)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N < 2 .)
!           If N=2, simply does linear interpolation.
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1)  <  X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of dependent variable values to be inter-
!           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
!           PCHIM is designed for monotonic data, but it will work for
!           any F-array.  It will force extrema at points where mono-
!           tonicity switches direction.  If some other treatment of
!           switch points is desired, PCHIC should be used instead.
!                                     -----
!     D -- (output) real array of derivative values at the data points.
!           If the data are monotonic, these values will determine a
!           a monotone cubic Hermite function.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD < 1 .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR > 0  means that IERR switches in the direction
!                 of monotonicity were detected.
!           "Recoverable" errors:
!              IERR = -1  if N < 2 .
!              IERR = -2  if INCFD < 1 .
!              IERR = -3  if the X-array is not strictly increasing.
!             (The D array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  1. F. N. Fritsch and J. Butland, A method for construc-
!                 ting local monotone piecewise cubic interpolants, SIAM
!                 Journal on Scientific and Statistical Computing 5, 2
!                 (June 1984), pp. 300-304.
!               2. F. N. Fritsch and R. E. Carlson, Monotone piecewise
!                 cubic interpolation, SIAM Journal on Numerical Ana-
!                 lysis 17, 2 (April 1980), pp. 238-246.
!***ROUTINES CALLED  PCHST, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811103  DATE WRITTEN
!   820201  1. Introduced  PCHST  to reduce possible over/under-
!             flow problems.
!           2. Rearranged derivative formula for same reason.
!   820602  1. Modified end conditions to be continuous functions
!             of data when monotonicity switches in next interval.
!           2. Modified formulas so end conditions are less prone
!             of over/underflow problems.
!   820803  Minor cosmetic changes for release 1.
!   870813  Updated Reference 1.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890703  Corrected category record.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920429  Revised format and order of references.  (WRB,FNF)
!***END PROLOGUE  PCHIM
!  Programming notes:
!
!     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.
!     2. To produce a double precision version, simply:
!        a. Change PCHIM to DPCHIM wherever it occurs,
!        b. Change PCHST to DPCHST wherever it occurs,
!        c. Change all references to the Fortran intrinsics to their
!           double precision equivalents,
!        d. Change the real declarations to double precision, and
!        e. Change the constants ZERO and THREE to double precision.
!
!  DECLARE ARGUMENTS.
!
  INTEGER  N, INCFD, IERR
  REAL  X(*), F(INCFD,*), D(INCFD,*)
!
!  DECLARE LOCAL VARIABLES.
!
  INTEGER  I, NLESS1
  REAL  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, DSAVE, &
        H1, H2, HSUM, HSUMT3, THREE, W1, W2, ZERO
  SAVE ZERO, THREE
  REAL  PCHST
  DATA  ZERO /0./,  THREE /3./
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  PCHIM
  if ( N < 2 )  go to 5001
  if ( INCFD < 1 )  go to 5002
  DO 1  I = 2, N
     if ( X(I) <= X(I-1) )  go to 5003
    1 CONTINUE
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
  IERR = 0
  NLESS1 = N - 1
  H1 = X(2) - X(1)
  DEL1 = (F(1,2) - F(1,1))/H1
  DSAVE = DEL1
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
  H2 = X(3) - X(2)
  DEL2 = (F(1,3) - F(1,2))/H2
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!     SHAPE-PRESERVING.
!
  HSUM = H1 + H2
  W1 = (H1 + HSUM)/HSUM
  W2 = -H1/HSUM
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
     H1 = H2
     H2 = X(I+1) - X(I)
     HSUM = H1 + H2
     DEL1 = DEL2
     DEL2 = (F(1,I+1) - F(1,I))/H2
   40    CONTINUE
!
!        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!
     D(1,I) = ZERO
     if ( PCHST(DEL1,DEL2) )  42, 41, 45
!
!        COUNT NUMBER OF CHANGES IN DIRECTION OF MONOTONICITY.
!
   41    CONTINUE
     if (DEL2  ==  ZERO)  go to 50
     if ( PCHST(DSAVE,DEL2)  <  ZERO)  IERR = IERR + 1
     DSAVE = DEL2
     go to 50
!
   42    CONTINUE
     IERR = IERR + 1
     DSAVE = DEL2
     go to 50
!
!        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
   45    CONTINUE
     HSUMT3 = HSUM+HSUM+HSUM
     W1 = (HSUM + H1)/HSUMT3
     W2 = (HSUM + H2)/HSUMT3
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
  W1 = -H2/HSUM
  W2 = (H2 + HSUM)/HSUM
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
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N < 2 RETURN.
  IERR = -1
  call XERMSG ('SLATEC', 'PCHIM', &
     'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  return
!
 5002 CONTINUE
!     INCFD < 1 RETURN.
  IERR = -2
  call XERMSG ('SLATEC', 'PCHIM', 'INCREMENT LESS THAN ONE', IERR, &
     1)
  return
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  call XERMSG ('SLATEC', 'PCHIM', 'X-ARRAY NOT STRICTLY INCREASING' &
     , IERR, 1)
  return
!------------- LAST LINE OF PCHIM FOLLOWS ------------------------------
end
