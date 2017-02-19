subroutine MPERR
!
!! MPERR is subsidiary to DQDOTA and DQDOTI.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (MPERR-A)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!  This routine is called when a fatal error condition is
!  encountered, and after a message has been written on
!  logical unit LUN.
!
!***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    MPCOM
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   ??????  Modified for use with BLAS.  Blank COMMON changed to named
!           COMMON.  R given dimension 12.
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
!***END PROLOGUE  MPERR
  COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
  INTEGER B, T, R
!***FIRST EXECUTABLE STATEMENT  MPERR
  call XERMSG('SLATEC', 'MPERR', &
     ' *** EXECUTION TERMINATED BY call TO MPERR' // &
     ' IN MP VERSION 770217 ***', 1, 2)
!
! AT PRESENT JUST STOP, BUT COULD DUMP B, T, ETC. HERE.
! ACTION COULD EASILY BE CONTROLLED BY A FLAG IN LABELLED COMMON.
! ANSI VERSION USES STOP, UNIVAC 1108 VERSION USES
! RETURN 0 IN ORDER TO GIVE A TRACE-BACK.
! FOR DEBUGGING PURPOSES IT MAY BE USEFUL SIMPLY TO
! RETURN HERE.  MOST MP ROUTINES RETURN WITH RESULT
! ZERO AFTER CALLING MPERR.
  STOP
end
