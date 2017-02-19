subroutine MACON
!
!! MACON is subsidiary to BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (MACON-S, DMACON-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!    Sets up machine constants using R1MACH
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  R1MACH
!***COMMON BLOCKS    ML5MCO
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  MACON
  COMMON /ML5MCO/ URO,SRU,EPS,SQOVFL,TWOU,FOURU,LPAR
!***FIRST EXECUTABLE STATEMENT  MACON
  URO=R1MACH(4)
  SRU=SQRT(URO)
  DD=-LOG10(URO)
  LPAR=0.5*DD
  KE=0.5+0.75*DD
  EPS=10.**(-2*KE)
  SQOVFL=SQRT(R1MACH(2))
  TWOU=2.0*URO
  FOURU=4.0*URO
  return
end
