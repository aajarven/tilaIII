subroutine INTRV (XT, LXT, X, ILO, ILEFT, MFLAG)
!
!! INTRV computes the largest integer ILEFT in 1  <=  ILEFT  <=  LXT ...
!            such that XT(ILEFT)  <=  X where XT(*) is a subdivision ...
!            of the X interval.
!
!***LIBRARY   SLATEC
!***CATEGORY  E3, K6
!***TYPE      SINGLE PRECISION (INTRV-S, DINTRV-D)
!***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Abstract
!         INTRV is the INTERV routine of the reference.
!
!         INTRV computes the largest integer ILEFT in 1  <=  ILEFT  <=
!         LXT such that XT(ILEFT)  <=  X where XT(*) is a subdivision of
!         the X interval.  Precisely,
!
!                      X  <  XT(1)                1         -1
!         if  XT(I)  <=  X  <  XT(I+1)  then  ILEFT=I  , MFLAG=0
!           XT(LXT)  <=  X                         LXT        1,
!
!         That is, when multiplicities are present in the break point
!         to the left of X, the largest index is taken for ILEFT.
!
!     Description of Arguments
!         Input
!          XT      - XT is a knot or break point vector of length LXT
!          LXT     - length of the XT vector
!          X       - argument
!          ILO     - an initialization parameter which must be set
!                    to 1 the first time the spline array XT is
!                    processed by INTRV.
!
!         Output
!          ILO     - ILO contains information for efficient process-
!                    ing after the initial call, and ILO must not be
!                    changed by the user.  Distinct splines require
!                    distinct ILO parameters.
!          ILEFT   - largest integer satisfying XT(ILEFT)  <=  X
!          MFLAG   - signals when X lies out of bounds
!
!     Error Conditions
!         None
!
!***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  INTRV
!
  INTEGER IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE
  REAL X, XT
  DIMENSION XT(*)
!***FIRST EXECUTABLE STATEMENT  INTRV
  IHI = ILO + 1
  if (IHI < LXT) go to 10
  if (X >= XT(LXT)) go to 110
  if (LXT <= 1) go to 90
  ILO = LXT - 1
  IHI = LXT
!
   10 if (X >= XT(IHI)) go to 40
  if (X >= XT(ILO)) go to 100
!
! *** NOW X  <  XT(IHI) . FIND LOWER BOUND
  ISTEP = 1
   20 IHI = ILO
  ILO = IHI - ISTEP
  if (ILO <= 1) go to 30
  if (X >= XT(ILO)) go to 70
  ISTEP = ISTEP*2
  go to 20
   30 ILO = 1
  if (X < XT(1)) go to 90
  go to 70
! *** NOW X  >=  XT(ILO) . FIND UPPER BOUND
   40 ISTEP = 1
   50 ILO = IHI
  IHI = ILO + ISTEP
  if (IHI >= LXT) go to 60
  if (X < XT(IHI)) go to 70
  ISTEP = ISTEP*2
  go to 50
   60 if (X >= XT(LXT)) go to 110
  IHI = LXT
!
! *** NOW XT(ILO)  <=  X  <  XT(IHI) . NARROW THE INTERVAL
   70 MIDDLE = (ILO+IHI)/2
  if (MIDDLE == ILO) go to 100
!     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
  if (X < XT(MIDDLE)) go to 80
  ILO = MIDDLE
  go to 70
   80 IHI = MIDDLE
  go to 70
! *** SET OUTPUT AND RETURN
   90 MFLAG = -1
  ILEFT = 1
  return
  100 MFLAG = 0
  ILEFT = ILO
  return
  110 MFLAG = 1
  ILEFT = LXT
  return
end
