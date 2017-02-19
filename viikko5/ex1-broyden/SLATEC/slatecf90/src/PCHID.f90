FUNCTION PCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERR)
!
!! PCHID evaluates the definite integral of a piecewise cubic Hermite ...
!  function over an interval whose endpoints are data points.
!
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3, H2A1B2
!***TYPE      SINGLE PRECISION (PCHID-S, DPCHID-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, NUMERICAL INTEGRATION, PCHIP,
!             QUADRATURE
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          PCHID:  Piecewise Cubic Hermite Integrator, Data limits
!
!     Evaluates the definite integral of the cubic Hermite function
!     defined by  N, X, F, D  over the interval [X(IA), X(IB)].
!
!     To provide compatibility with PCHIM and PCHIC, includes an
!     increment between successive values of the F- and D arrays.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IA, IB, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N)
!        LOGICAL  SKIP
!
!        VALUE = PCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERR)
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
!           SKIP will be set to .TRUE. on return with IERR = 0 or -4.
!
!     IA,IB -- (input) indices in X-array for the limits of integration.
!           both must be in the range [1,N].  (Error return if not.)
!           No restrictions on their relative values.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N < 2 .
!              IERR = -2  if INCFD < 1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if IA or IB is out of range.
!                (VALUE will be zero in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820723  DATE WRITTEN
!   820804  Converted to SLATEC library version.
!   870813  Minor cosmetic changes.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890703  Corrected category record.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   930504  Corrected to set VALUE=0 when IERR.ne.0.  (FNF)
!***END PROLOGUE  PCHID
!
!  Programming notes:
!  1. This routine uses a special formula that is valid only for
!     integrals whose limits coincide with data values.  This is
!     mathematically equivalent to, but much more efficient than,
!     calls to CHFIE.
!**End
!
!  DECLARE ARGUMENTS.
!
  REAL PCHID
  INTEGER  N, INCFD, IA, IB, IERR
  REAL  X(*), F(INCFD,*), D(INCFD,*)
  LOGICAL  SKIP
!
!  DECLARE LOCAL VARIABLES.
!
  INTEGER  I, IUP, LOW
  REAL  H, HALF, SIX, SUM, VALUE, ZERO
  SAVE ZERO, HALF, SIX
!
!  INITIALIZE.
!
  DATA  ZERO /0./,  HALF /0.5/,  SIX /6./
!***FIRST EXECUTABLE STATEMENT  PCHID
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
  if ((IA < 1) .OR. (IA > N))  go to 5004
  if ((IB < 1) .OR. (IB > N))  go to 5004
  IERR = 0
!
!  COMPUTE INTEGRAL VALUE.
!
  if (IA  /=  IB)  THEN
     LOW = MIN(IA, IB)
     IUP = MAX(IA, IB) - 1
     SUM = ZERO
     DO 10  I = LOW, IUP
        H = X(I+1) - X(I)
        SUM = SUM + H*( (F(1,I) + F(1,I+1)) + &
                        (D(1,I) - D(1,I+1))*(H/SIX) )
   10    CONTINUE
     VALUE = HALF * SUM
     if (IA  >  IB)  VALUE = -VALUE
  end if
!
!  NORMAL RETURN.
!
 5000 CONTINUE
  PCHID = VALUE
  return
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N < 2 RETURN.
  IERR = -1
  call XERMSG ('SLATEC', 'PCHID', &
     'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  go to 5000
!
 5002 CONTINUE
!     INCFD < 1 RETURN.
  IERR = -2
  call XERMSG ('SLATEC', 'PCHID', 'INCREMENT LESS THAN ONE', IERR, &
     1)
  go to 5000
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  call XERMSG ('SLATEC', 'PCHID', &
     'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  go to 5000
!
 5004 CONTINUE
!     IA OR IB OUT OF RANGE RETURN.
  IERR = -4
  call XERMSG ('SLATEC', 'PCHID', 'IA OR IB OUT OF RANGE', IERR, 1)
  go to 5000
!------------- LAST LINE OF PCHID FOLLOWS ------------------------------
end
