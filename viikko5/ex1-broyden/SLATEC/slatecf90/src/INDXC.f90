subroutine INDXC (I, IR, IDXC, NC)
!
!! INDXC is subsidiary to BLKTRI.
!
!***LIBRARY   SLATEC
!***TYPE      INTEGER (INDXC-I)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  BLKTRI
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    CBLKT
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  INDXC
  COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        , &
                  NM         ,NCMPLX     ,IK
!***FIRST EXECUTABLE STATEMENT  INDXC
  NC = 2**IR
  IDXC = I
  if (IDXC+NC-1-NM) 102,102,101
  101 NC = 0
  102 RETURN
end
