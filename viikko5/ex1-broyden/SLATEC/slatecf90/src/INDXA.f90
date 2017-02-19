subroutine INDXA (I, IR, IDXA, NA)
!
!! INDXA is subsidiary to BLKTRI.
!
!***LIBRARY   SLATEC
!***TYPE      INTEGER (INDXA-I)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  BLKTRI
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    CBLKT
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  INDXA
  COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        , &
                  NM         ,NCMPLX     ,IK
!***FIRST EXECUTABLE STATEMENT  INDXA
  NA = 2**IR
  IDXA = I-NA+1
  if (I-NM) 102,102,101
  101 NA = 0
  102 RETURN
end
