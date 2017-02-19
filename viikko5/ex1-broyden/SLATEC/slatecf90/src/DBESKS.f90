subroutine DBESKS (XNU, X, NIN, BK)
!
!! DBESKS computes a sequence of modified Bessel functions of the ...
!            third kind of fractional order.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B3
!***TYPE      DOUBLE PRECISION (BESKS-S, DBESKS-D)
!***KEYWORDS  FNLIB, FRACTIONAL ORDER, MODIFIED BESSEL FUNCTION,
!             SEQUENCE OF BESSEL FUNCTIONS, SPECIAL FUNCTIONS,
!             THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBESKS computes a sequence of modified Bessel functions of the third
! kind of order XNU + I at X, where X  >  0, XNU lies in (-1,1),
! and I = 0, 1, ... , NIN - 1, if NIN is positive and I = 0, 1, ... ,
! NIN + 1, if NIN is negative.  On return, the vector BK(.) contains
! the results at X for order starting at XNU.  XNU, X, and BK are
! double precision.  NIN is an integer.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DBSKES, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBESKS
  DOUBLE PRECISION XNU, X, BK(*), EXPXI, XMAX, D1MACH
  SAVE XMAX
  DATA XMAX / 0.D0 /
!***FIRST EXECUTABLE STATEMENT  DBESKS
  if (XMAX == 0.D0) XMAX = -LOG (D1MACH(1))
!
  if (X  >  XMAX) call XERMSG ('SLATEC', 'DBESKS', &
     'X SO BIG BESSEL K UNDERFLOWS', 1, 2)
!
  call DBSKES (XNU, X, NIN, BK)
!
  EXPXI = EXP (-X)
  N = ABS (NIN)
  DO 20 I=1,N
    BK(I) = EXPXI * BK(I)
 20   CONTINUE
!
  return
end
