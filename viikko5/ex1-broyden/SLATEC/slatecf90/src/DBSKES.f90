subroutine DBSKES (XNU, X, NIN, BKE)
!
!! DBSKES computes a sequence of exponentially scaled modified Bessel ...
!  functions of the third kind of fractional order.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B3
!***TYPE      DOUBLE PRECISION (BESKES-S, DBSKES-D)
!***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, FRACTIONAL ORDER,
!             MODIFIED BESSEL FUNCTION, SEQUENCE OF BESSEL FUNCTIONS,
!             SPECIAL FUNCTIONS, THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBSKES(XNU,X,NIN,BKE) computes a double precision sequence
! of exponentially scaled modified Bessel functions
! of the third kind of order XNU + I at X, where X  >  0,
! XNU lies in (-1,1), and I = 0, 1, ... , NIN - 1, if NIN is positive
! and I = 0, -1, ... , NIN + 1, if NIN is negative.  On return, the
! vector BKE(.) contains the results at X for order starting at XNU.
! XNU, X, and BKE are double precision.  NIN is integer.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, D9KNUS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBSKES
  DOUBLE PRECISION XNU, X, BKE(*), BKNU1, V, VINCR, VEND, ALNBIG, &
    D1MACH, DIRECT
  SAVE ALNBIG
  DATA ALNBIG / 0.D0 /
!***FIRST EXECUTABLE STATEMENT  DBSKES
  if (ALNBIG == 0.D0) ALNBIG = LOG (D1MACH(2))
!
  V = ABS(XNU)
  N = ABS(NIN)
!
  if (V  >=  1.D0) call XERMSG ('SLATEC', 'DBSKES', &
     'ABS(XNU) MUST BE LT 1', 2, 2)
  if (X  <=  0.D0) call XERMSG ('SLATEC', 'DBSKES', 'X IS LE 0', 3, &
     2)
  if (N  ==  0) call XERMSG ('SLATEC', 'DBSKES', &
     'N THE NUMBER IN THE SEQUENCE IS 0', 4, 2)
!
  call D9KNUS (V, X, BKE(1), BKNU1, ISWTCH)
  if (N == 1) RETURN
!
  VINCR = SIGN (1.0, REAL(NIN))
  DIRECT = VINCR
  if (XNU /= 0.D0) DIRECT = VINCR*SIGN(1.D0, XNU)
  if (ISWTCH  ==  1 .AND. DIRECT  >  0.) call XERMSG ('SLATEC', &
     'DBSKES', 'X SO SMALL BESSEL K-SUB-XNU+1 OVERFLOWS', 5, 2)
  BKE(2) = BKNU1
!
  if (DIRECT < 0.) call D9KNUS (ABS(XNU+VINCR), X, BKE(2), BKNU1, &
    ISWTCH)
  if (N == 2) RETURN
!
  VEND = ABS (XNU+NIN) - 1.0D0
  if ((VEND-.5D0)*LOG(VEND)+0.27D0-VEND*(LOG(X)-.694D0)  >  &
     ALNBIG) call XERMSG ('SLATEC', 'DBSKES', &
        'X SO SMALL OR ABS(NU) SO BIG THAT BESSEL K-SUB-NU ' // &
        'OVERFLOWS', 5, 2)
!
  V = XNU
  DO 10 I=3,N
    V = V + VINCR
    BKE(I) = 2.0D0*V*BKE(I-1)/X + BKE(I-2)
 10   CONTINUE
!
  return
end
