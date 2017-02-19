FUNCTION C9LN2R (Z)
!
!! C9LN2R evaluates LOG(1+Z) from second order relative accuracy so ...
!  that  LOG(1+Z) = Z - Z**2/2 + Z**3*C9LN2R(Z).
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4B
!***TYPE      COMPLEX (R9LN2R-S, D9LN2R-D, C9LN2R-C)
!***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM, SECOND ORDER
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate  LOG(1+Z)  from 2-nd order with relative error accuracy so
! that     LOG(1+Z) = Z - Z**2/2 + Z**3*C9LN2R(Z).
!
! Now  LOG(1+Z) = 0.5*LOG(1+2*X+ABS(Z)**2) + I*CARG(1+Z),
! where X = REAL(Z)  and  Y = AIMAG(Z).
! We find
!     Z**3 * C9LN2R(Z) = -X*ABS(Z)**2 - 0.25*ABS(Z)**4
!        + (2*X+ABS(Z)**2)**3 * R9LN2R(2*X+ABS(Z)**2)
!        + I * (CARG(1+Z) + (X-1)*Y)
! The imaginary part must be evaluated carefully as
!     (ATAN(Y/(1+X)) - Y/(1+X)) + Y/(1+X) - (1-X)*Y
!       = (Y/(1+X))**3 * R9ATN1(Y/(1+X)) + X**2*Y/(1+X)
!
! Now we divide through by Z**3 carefully.  Write
!     1/Z**3 = (X-I*Y)/ABS(Z)**3 * (1/ABS(Z)**3)
! then   C9LN2R(Z) = ((X-I*Y)/ABS(Z))**3 * (-X/ABS(Z) - ABS(Z)/4
!        + 0.5*((2*X+ABS(Z)**2)/ABS(Z))**3 * R9LN2R(2*X+ABS(Z)**2)
!        + I*Y/(ABS(Z)*(1+X)) * ((X/ABS(Z))**2 +
!          + (Y/(ABS(Z)*(1+X)))**2 * R9ATN1(Y/(1+X)) ) )
!
! If we let  XZ = X/ABS(Z)  and  YZ = Y/ABS(Z)  we may write
!     C9LN2R(Z) = (XZ-I*YZ)**3 * (-XZ - ABS(Z)/4
!        + 0.5*(2*XZ+ABS(Z))**3 * R9LN2R(2*X+ABS(Z)**2)
!        + I*YZ/(1+X) * (XZ**2 + (YZ/(1+X))**2*R9ATN1(Y/(1+X)) ))
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R9ATN1, R9LN2R
!***REVISION HISTORY  (YYMMDD)
!   780401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  C9LN2R
  COMPLEX C9LN2R
  COMPLEX Z
!***FIRST EXECUTABLE STATEMENT  C9LN2R
  X = REAL (Z)
  Y = AIMAG (Z)
!
  CABSZ = ABS(Z)
  if (CABSZ > 0.8125) go to 20
!
  C9LN2R = CMPLX (1.0/3.0, 0.0)
  if (CABSZ == 0.0) RETURN
!
  XZ = X/CABSZ
  YZ = Y/CABSZ
!
  ARG = 2.0*XZ + CABSZ
  RPART = 0.5*ARG**3*R9LN2R(CABSZ*ARG) - XZ - 0.25*CABSZ
  Y1X = YZ/(1.0+X)
  AIPART = Y1X * (XZ**2 + Y1X**2*R9ATN1(CABSZ*Y1X) )
!
  C9LN2R = CMPLX(XZ,-YZ)**3 * CMPLX(RPART,AIPART)
  return
!
 20   C9LN2R = (LOG(1.0+Z) - Z*(1.0-0.5*Z)) / Z**3
  return
!
end
