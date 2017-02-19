subroutine INXCB (I, IR, IDX, IDP)
!
!! INXCB is subsidiary to CBLKTR.
!
!***LIBRARY   SLATEC
!***TYPE      INTEGER (INXCB-I)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  CBLKTR
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    CCBLK
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  INXCB
!
  COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        , &
                  NM         ,NCMPLX     ,IK
!***FIRST EXECUTABLE STATEMENT  INXCB
  IDP = 0
  if (IR) 107,101,103
  101 if (I-NM) 102,102,107
  102 IDX = I
  IDP = 1
  return
  103 IZH = 2**IR
  ID = I-IZH-IZH
  IDX = ID+ID+(IR-1)*IK+IR+(IK-I)/IZH+4
  IPL = IZH-1
  IDP = IZH+IZH-1
  if (I-IPL-NM) 105,105,104
  104 IDP = 0
  return
  105 if (I+IPL-NM) 107,107,106
  106 IDP = NM+IPL-I+1
  107 RETURN
end
