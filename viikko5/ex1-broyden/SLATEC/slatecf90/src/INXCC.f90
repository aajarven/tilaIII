subroutine INXCC (I, IR, IDXC, NC)
!
!! INXCC is subsidiary to CBLKTR.
!
!***LIBRARY   SLATEC
!***TYPE      INTEGER (INXCC-I)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  CBLKTR
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    CCBLK
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  INXCC
  COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        , &
                  NM         ,NCMPLX     ,IK
!***FIRST EXECUTABLE STATEMENT  INXCC
  NC = 2**IR
  IDXC = I
  if (IDXC+NC-1-NM) 102,102,101
  101 NC = 0
  102 RETURN
end
