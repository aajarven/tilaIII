function PIMACH (DUM)
!
!! PIMACH supplies the value of PI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PIMACH-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subprogram supplies the value of the constant PI correct to
!     machine precision where
!
!     PI=3.1415926535897932384626433832795028841971693993751058209749446
!
!***SEE ALSO  HSTCSP, HSTSSP, HWSCSP
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PIMACH
!
!***FIRST EXECUTABLE STATEMENT  PIMACH
  PIMACH = 3.14159265358979
  return
end
