subroutine INXCA (I, IR, IDXA, NA)
!
!! INXCA is subsidiary to CBLKTR.
!
!***LIBRARY   SLATEC
!***TYPE      INTEGER (INXCA-I)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  CBLKTR
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    CCBLK
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  INXCA
  COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        , &
                  NM         ,NCMPLX     ,IK
!***FIRST EXECUTABLE STATEMENT  INXCA
  NA = 2**IR
  IDXA = I-NA+1
  if (I-NM) 102,102,101
  101 NA = 0
  102 RETURN
end
