FUNCTION PCHIA (N, X, F, D, INCFD, SKIP, A, B, IERR)
!
!! PCHIA evaluates the definite integral of a piecewise cubic Hermite ...
!  function over an arbitrary interval.
!
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3, H2A1B2
!***TYPE      SINGLE PRECISION (PCHIA-S, DPCHIA-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, NUMERICAL INTEGRATION, PCHIP,
!             QUADRATURE
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          PCHIA:  Piecewise Cubic Hermite Integrator, Arbitrary limits
!
!     Evaluates the definite integral of the cubic Hermite function
!     defined by  N, X, F, D  over the interval [A, B].
!
!     To provide compatibility with PCHIM and PCHIC, includes an
!     increment between successive values of the F- and D arrays.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N), A, B
!        REAL  VALUE, PCHIA
!        LOGICAL  SKIP
!
!        VALUE = PCHIA (N, X, F, D, INCFD, SKIP, A, B, IERR)
!
!   Parameters:
!
!     VALUE -- (output) value of the requested integral.
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
!           SKIP will be set to .TRUE. on return with IERR >= 0 .
!
!     A,B -- (input) the limits of integration.
!           NOTE:  There is no requirement that [A,B] be contained in
!                  [X(1),X(N)].  However, the resulting integral value
!                  will be highly suspect, if not.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if  A  is outside the interval [X(1),X(N)].
!              IERR = 2  if  B  is outside the interval [X(1),X(N)].
!              IERR = 3  if both of the above are true.  (Note that this
!                        means that either [A,B] contains data interval
!                        or the intervals do not intersect at all.)
!           "Recoverable" errors:
!              IERR = -1  if N < 2 .
!              IERR = -2  if INCFD < 1 .
!              IERR = -3  if the X-array is not strictly increasing.
!                (VALUE will be zero in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!              IERR = -4  in case of an error return from PCHID (which
!                         should never occur).
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CHFIE, PCHID, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820730  DATE WRITTEN
!   820804  Converted to SLATEC library version.
!   870707  Corrected double precision conversion instructions.
!   870813  Minor cosmetic changes.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890703  Corrected category record.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   930503  Corrected to set VALUE=0 when IERR.lt.0.  (FNF)
!   930504  Changed CHFIV to CHFIE.  (FNF)
!***END PROLOGUE  PCHIA
!
!  Programming notes:
!  1. The error flag from PCHID is tested, because a logic flaw
!     could conceivably result in IERD=-4, which should be reported.
!**End
!
!  DECLARE ARGUMENTS.
!
  INTEGER  N, INCFD, IERR
  REAL PCHIA
  REAL  X(*), F(INCFD,*), D(INCFD,*), A, B
  LOGICAL  SKIP
!
!  DECLARE LOCAL VARIABLES.
!
  INTEGER  I, IA, IB, IERD, IL, IR
  REAL  VALUE, XA, XB, ZERO
  SAVE ZERO
  REAL  CHFIE, PCHID
!
!  INITIALIZE.
!
  DATA  ZERO /0./
!***FIRST EXECUTABLE STATEMENT  PCHIA
  VALUE = ZERO
!
!  VALIDITY-CHECK ARGUMENTS.
!
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
  SKIP = .TRUE.
  IERR = 0
  if ( (A < X(1)) .OR. (A > X(N)) )  IERR = IERR + 1
  if ( (B < X(1)) .OR. (B > X(N)) )  IERR = IERR + 2
!
!  COMPUTE INTEGRAL VALUE.
!
  if (A  /=  B)  THEN
     XA = MIN (A, B)
     XB = MAX (A, B)
     if (XB  <=  X(2))  THEN
!           INTERVAL IS TO LEFT OF X(2), SO USE FIRST CUBIC.
!                   --------------------------------------
        VALUE = CHFIE (X(1),X(2), F(1,1),F(1,2), &
                                  D(1,1),D(1,2), A, B)
!                   --------------------------------------
     ELSE if (XA  >=  X(N-1))  THEN
!           INTERVAL IS TO RIGHT OF X(N-1), SO USE LAST CUBIC.
!                   -----------------------------------------
        VALUE = CHFIE(X(N-1),X(N), F(1,N-1),F(1,N), &
                                   D(1,N-1),D(1,N), A, B)
!                   -----------------------------------------
     ELSE
!           'NORMAL' CASE -- XA < XB, XA < X(N-1), XB > X(2).
!      ......LOCATE IA AND IB SUCH THAT
!               X(IA-1) < XA <= X(IA) <= X(IB) <= XB <= X(IB+1)
        IA = 1
        DO 10  I = 1, N-1
           if (XA  >  X(I))  IA = I + 1
   10       CONTINUE
!             IA = 1 IMPLIES XA < X(1) .  OTHERWISE,
!             IA IS LARGEST INDEX SUCH THAT X(IA-1) < XA,.
!
        IB = N
        DO 20  I = N, IA, -1
           if (XB  <  X(I))  IB = I - 1
   20       CONTINUE
!             IB = N IMPLIES XB > X(N) .  OTHERWISE,
!             IB IS SMALLEST INDEX SUCH THAT XB < X(IB+1) .
!
!     ......COMPUTE THE INTEGRAL.
        if (IB  <  IA)  THEN
!              THIS MEANS IB = IA-1 AND
!                 (A,B) IS A SUBSET OF (X(IB),X(IA)).
!                      ------------------------------------------
           VALUE = CHFIE (X(IB),X(IA), F(1,IB),F(1,IA), &
                                       D(1,IB),D(1,IA), A, B)
!                      ------------------------------------------
        ELSE
!
!              FIRST COMPUTE INTEGRAL OVER (X(IA),X(IB)).
!                (Case (IB  ==  IA) is taken care of by initialization
!                 of VALUE to ZERO.)
           if (IB  >  IA)  THEN
!                         ---------------------------------------------
              VALUE = PCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERD)
!                         ---------------------------------------------
              if (IERD  <  0)  go to 5004
           ENDIF
!
!              THEN ADD ON INTEGRAL OVER (XA,X(IA)).
           if (XA  <  X(IA))  THEN
              IL = MAX(1, IA-1)
              IR = IL + 1
!                                 -------------------------------------
              VALUE = VALUE + CHFIE (X(IL),X(IR), F(1,IL),F(1,IR), &
                                        D(1,IL),D(1,IR), XA, X(IA))
!                                 -------------------------------------
           ENDIF
!
!              THEN ADD ON INTEGRAL OVER (X(IB),XB).
           if (XB  >  X(IB))  THEN
              IR = MIN (IB+1, N)
              IL = IR - 1
!                                 -------------------------------------
              VALUE = VALUE + CHFIE (X(IL),X(IR), F(1,IL),F(1,IR), &
                                        D(1,IL),D(1,IR), X(IB), XB)
!                                 -------------------------------------
           ENDIF
!
!              FINALLY, ADJUST SIGN if NECESSARY.
           if (A  >  B)  VALUE = -VALUE
        ENDIF
     ENDIF
  end if
!
!  NORMAL RETURN.
!
 5000 CONTINUE
  PCHIA = VALUE
  return
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N < 2 RETURN.
  IERR = -1
  call XERMSG ('SLATEC', 'PCHIA', &
     'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  go to 5000
!
 5002 CONTINUE
!     INCFD < 1 RETURN.
  IERR = -2
  call XERMSG ('SLATEC', 'PCHIA', 'INCREMENT LESS THAN ONE', IERR, &
     1)
  go to 5000
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  call XERMSG ('SLATEC', 'PCHIA', &
     'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  go to 5000
!
 5004 CONTINUE
!     TROUBLE IN PCHID.  (SHOULD NEVER OCCUR.)
  IERR = -4
  call XERMSG ('SLATEC', 'PCHIA', 'TROUBLE IN PCHID', IERR, 1)
  go to 5000
!------------- LAST LINE OF PCHIA FOLLOWS ------------------------------
end
