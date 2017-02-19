subroutine MPUNFL (X)
!
!! MPUNFL is called to handle multiple precision underflow.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DQDOTA and DQDOTI
!***LIBRARY   SLATEC
!***TYPE      ALL (MPUNFL-A)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
! Called on multiple-precision underflow, i.e.  when the
! exponent of 'mp' number X would be less than -M.
!
!***SEE ALSO  DQDOTA, DQDOTI
!***ROUTINES CALLED  MPCHK
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  MPUNFL
  INTEGER X(*)
!***FIRST EXECUTABLE STATEMENT  MPUNFL
  call MPCHK (1, 4)
! THE UNDERFLOWING NUMBER IS SET TO ZERO
! AN ALTERNATIVE WOULD BE TO call MPMINR (X) AND RETURN,
! POSSIBLY UPDATING A COUNTER AND TERMINATING EXECUTION
! AFTER A PRESET NUMBER OF UNDERFLOWS.  ACTION COULD EASILY
! BE DETERMINED BY A FLAG IN LABELLED COMMON.
  X(1) = 0
  return
end
