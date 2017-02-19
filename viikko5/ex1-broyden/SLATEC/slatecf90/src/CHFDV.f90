subroutine CHFDV (X1, X2, F1, F2, D1, D2, NE, XE, FE, DE, NEXT, &
     IERR)
!
!! CHFDV evaluates a cubic polynomial given in Hermite form and its ...
!            first derivative at an array of points.  While designed for ...
!            use by PCHFD, it may be useful directly as an evaluator ...
!            for a piecewise cubic Hermite function in applications, ...
!            such as graphing, where the interval is known in advance. ...
!            If only function values are required, use CHFEV instead.
!
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3, H1
!***TYPE      SINGLE PRECISION (CHFDV-S, DCHFDV-D)
!***KEYWORDS  CUBIC HERMITE DIFFERENTIATION, CUBIC HERMITE EVALUATION,
!             CUBIC POLYNOMIAL EVALUATION, PCHIP
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!        CHFDV:  Cubic Hermite Function and Derivative Evaluator
!
!     Evaluates the cubic polynomial determined by function values
!     F1,F2 and derivatives D1,D2 on interval (X1,X2), together with
!     its first derivative, at the points  XE(J), J=1(1)NE.
!
!     If only function values are required, use CHFEV, instead.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        INTEGER  NE, NEXT(2), IERR
!        REAL  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE), DE(NE)
!
!        call  CHFDV (X1,X2, F1,F2, D1,D2, NE, XE, FE, DE, NEXT, IERR)
!
!   Parameters:
!
!     X1,X2 -- (input) endpoints of interval of definition of cubic.
!           (Error return if  X1 == X2 .)
!
!     F1,F2 -- (input) values of function at X1 and X2, respectively.
!
!     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE < 1 .)
!
!     XE -- (input) real array of points at which the functions are to
!           be evaluated.  If any of the XE are outside the interval
!           [X1,X2], a warning error is returned in NEXT.
!
!     FE -- (output) real array of values of the cubic function defined
!           by  X1,X2, F1,F2, D1,D2  at the points  XE.
!
!     DE -- (output) real array of values of the first derivative of
!           the same function at the points  XE.
!
!     NEXT -- (output) integer array indicating number of extrapolation
!           points:
!            NEXT(1) = number of evaluation points to left of interval.
!            NEXT(2) = number of evaluation points to right of interval.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if NE < 1 .
!              IERR = -2  if X1 == X2 .
!                (Output arrays have not been changed in either case.)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811019  DATE WRITTEN
!   820803  Minor cosmetic changes for release 1.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  CHFDV
!  Programming notes:
!
!     To produce a double precision version, simply:
!        a. Change CHFDV to DCHFDV wherever it occurs,
!        b. Change the real declaration to double precision, and
!        c. Change the constant ZERO to double precision.
!
!  DECLARE ARGUMENTS.
!
  INTEGER  NE, NEXT(2), IERR
  REAL  X1, X2, F1, F2, D1, D2, XE(*), FE(*), DE(*)
!
!  DECLARE LOCAL VARIABLES.
!
  INTEGER  I
  REAL  C2, C2T2, C3, C3T3, DEL1, DEL2, DELTA, H, X, XMI, XMA, ZERO
  SAVE ZERO
  DATA  ZERO /0./
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  CHFDV
  if (NE  <  1)  go to 5001
  H = X2 - X1
  if (H  ==  ZERO)  go to 5002
!
!  INITIALIZE.
!
  IERR = 0
  NEXT(1) = 0
  NEXT(2) = 0
  XMI = MIN(ZERO, H)
  XMA = MAX(ZERO, H)
!
!  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
!
  DELTA = (F2 - F1)/H
  DEL1 = (D1 - DELTA)/H
  DEL2 = (D2 - DELTA)/H
!                                           (DELTA IS NO LONGER NEEDED.)
  C2 = -(DEL1+DEL1 + DEL2)
  C2T2 = C2 + C2
  C3 = (DEL1 + DEL2)/H
!                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
  C3T3 = C3+C3+C3
!
!  EVALUATION LOOP.
!
  DO 500  I = 1, NE
     X = XE(I) - X1
     FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
     DE(I) = D1 + X*(C2T2 + X*C3T3)
!          COUNT EXTRAPOLATION POINTS.
     if ( X < XMI )  NEXT(1) = NEXT(1) + 1
     if ( X > XMA )  NEXT(2) = NEXT(2) + 1
!        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
  500 CONTINUE
!
!  NORMAL RETURN.
!
  return
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     NE < 1 RETURN.
  IERR = -1
  call XERMSG ('SLATEC', 'CHFDV', &
     'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
  return
!
 5002 CONTINUE
!     X1 == X2 RETURN.
  IERR = -2
  call XERMSG ('SLATEC', 'CHFDV', 'INTERVAL ENDPOINTS EQUAL', IERR, &
     1)
  return
!------------- LAST LINE OF CHFDV FOLLOWS ------------------------------
end
