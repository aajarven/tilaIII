subroutine PCHFE (N, X, F, D, INCFD, SKIP, NE, XE, FE, IERR)
!
!! PCHFE evaluates a piecewise cubic Hermite function at an array of points. ...
!  May be used by itself for Hermite interpolation,
!            or as an evaluator for PCHIM or PCHIC.
!
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3
!***TYPE      SINGLE PRECISION (PCHFE-S, DPCHFE-D)
!***KEYWORDS  CUBIC HERMITE EVALUATION, HERMITE INTERPOLATION, PCHIP,
!             PIECEWISE CUBIC EVALUATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          PCHFE:  Piecewise Cubic Hermite Function Evaluator
!
!     Evaluates the cubic Hermite function defined by  N, X, F, D  at
!     the points  XE(J), J=1(1)NE.
!
!     To provide compatibility with PCHIM and PCHIC, includes an
!     increment between successive values of the F- and D arrays.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, NE, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE)
!        LOGICAL  SKIP
!
!        call  PCHFE (N, X, F, D, INCFD, SKIP, NE, XE, FE, IERR)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N < 2 .)
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1)  <  X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD < 1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in PCHIM or PCHIC).
!           SKIP will be set to .TRUE. on normal return.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE < 1 .)
!
!     XE -- (input) real array of points at which the function is to be
!           evaluated.
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XE are increasing relative to X;
!              that is,   XE(J)  >=  X(I)
!              implies    XE(K)  >=  X(I),  all K >= J .
!           2. If any of the XE are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     FE -- (output) real array of values of the cubic Hermite function
!           defined by  N, X, F, D  at the points  XE.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR > 0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N < 2 .
!              IERR = -2  if INCFD < 1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NE < 1 .
!             (The FE-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CHFEV, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811020  DATE WRITTEN
!   820803  Minor cosmetic changes for release 1.
!   870707  Minor cosmetic changes to prologue.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  PCHFE
!  Programming notes:
!
!     1. To produce a double precision version, simply:
!        a. Change PCHFE to DPCHFE, and CHFEV to DCHFEV, wherever they
!           occur,
!        b. Change the real declaration to double precision,
!
!     2. Most of the coding between the call to CHFEV and the end of
!        the IR-loop could be eliminated if it were permissible to
!        assume that XE is ordered relative to X.
!
!     3. CHFEV does not assume that X1 is less than X2.  thus, it would
!        be possible to write a version of PCHFE that assumes a strict-
!        ly decreasing X-array by simply running the IR-loop backwards
!        (and reversing the order of appropriate tests).
!
!     4. The present code has a minor bug, which I have decided is not
!        worth the effort that would be required to fix it.
!        If XE contains points in [X(N-1),X(N)], followed by points  <
!        X(N-1), followed by points  > X(N), the extrapolation points
!        will be counted (at least) twice in the total returned in IERR.
!
!  DECLARE ARGUMENTS.
!
  INTEGER  N, INCFD, NE, IERR
  REAL  X(*), F(INCFD,*), D(INCFD,*), XE(*), FE(*)
  LOGICAL  SKIP
!
!  DECLARE LOCAL VARIABLES.
!
  INTEGER  I, IERC, IR, J, JFIRST, NEXT(2), NJ
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  PCHFE
  if (SKIP)  go to 5
!
  if ( N < 2 )  go to 5001
  if ( INCFD < 1 )  go to 5002
  DO 1  I = 2, N
     if ( X(I) <= X(I-1) )  go to 5003
    1 CONTINUE
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
    5 CONTINUE
  if ( NE < 1 )  go to 5004
  IERR = 0
  SKIP = .TRUE.
!
!  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
!                              ( INTERVAL IS X(IL) <= X < X(IR) . )
  JFIRST = 1
  IR = 2
   10 CONTINUE
!
!     SKIP OUT OF LOOP if HAVE PROCESSED ALL EVALUATION POINTS.
!
     if (JFIRST  >  NE)  go to 5000
!
!     LOCATE ALL POINTS IN INTERVAL.
!
     DO 20  J = JFIRST, NE
        if (XE(J)  >=  X(IR))  go to 30
   20    CONTINUE
     J = NE + 1
     go to 40
!
!     HAVE LOCATED FIRST POINT BEYOND INTERVAL.
!
   30    CONTINUE
     if (IR  ==  N)  J = NE + 1
!
   40    CONTINUE
     NJ = J - JFIRST
!
!     SKIP EVALUATION if NO POINTS IN INTERVAL.
!
     if (NJ  ==  0)  go to 50
!
!     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
!
!       ----------------------------------------------------------------
    call CHFEV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR), &
                NJ, XE(JFIRST), FE(JFIRST), NEXT, IERC)
!       ----------------------------------------------------------------
     if (IERC  <  0)  go to 5005
!
     if (NEXT(2)  ==  0)  go to 42
!        if (NEXT(2)  >  0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
!           RIGHT OF X(IR).
!
        if (IR  <  N)  go to 41
!           if (IR  ==  N)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
           IERR = IERR + NEXT(2)
           go to 42
   41       CONTINUE
!           ELSE
!              WE SHOULD NEVER HAVE GOTTEN HERE.
           go to 5005
!           ENDIF
!        ENDIF
   42    CONTINUE
!
     if (NEXT(1)  ==  0)  go to 49
!        if (NEXT(1)  >  0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
!           LEFT OF X(IR-1).
!
        if (IR  >  2)  go to 43
!           if (IR  ==  2)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
           IERR = IERR + NEXT(1)
           go to 49
   43       CONTINUE
!           ELSE
!              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
!              EVALUATION INTERVAL.
!
!              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
           DO 44  I = JFIRST, J-1
              if (XE(I)  <  X(IR-1))  go to 45
   44          CONTINUE
!              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
!                     IN CHFEV.
           go to 5005
!
   45          CONTINUE
!              RESET J.  (THIS WILL BE THE NEW JFIRST.)
           J = I
!
!              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
           DO 46  I = 1, IR-1
              if (XE(J)  <  X(I)) go to 47
   46          CONTINUE
!              NB-- CAN NEVER DROP THROUGH HERE, SINCE XE(J) < X(IR-1).
!
   47          CONTINUE
!              AT THIS POINT, EITHER  XE(J)  <  X(1)
!                 OR      X(I-1)  <=  XE(J)  <  X(I) .
!              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
!              CYCLING.
           IR = MAX(1, I-1)
!           ENDIF
!        ENDIF
   49    CONTINUE
!
     JFIRST = J
!
!     END OF IR-LOOP.
!
   50 CONTINUE
  IR = IR + 1
  if (IR  <=  N)  go to 10
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
  call XERMSG ('SLATEC', 'PCHFE', &
     'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  return
!
 5002 CONTINUE
!     INCFD < 1 RETURN.
  IERR = -2
  call XERMSG ('SLATEC', 'PCHFE', 'INCREMENT LESS THAN ONE', IERR, &
     1)
  return
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  call XERMSG ('SLATEC', 'PCHFE', 'X-ARRAY NOT STRICTLY INCREASING' &
     , IERR, 1)
  return
!
 5004 CONTINUE
!     NE < 1 RETURN.
  IERR = -4
  call XERMSG ('SLATEC', 'PCHFE', &
     'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
  return
!
 5005 CONTINUE
!     ERROR RETURN FROM CHFEV.
!   *** THIS CASE SHOULD NEVER OCCUR ***
  IERR = -5
  call XERMSG ('SLATEC', 'PCHFE', &
     'ERROR RETURN FROM CHFEV -- FATAL', IERR, 2)
  return
!------------- LAST LINE OF PCHFE FOLLOWS ------------------------------
end
