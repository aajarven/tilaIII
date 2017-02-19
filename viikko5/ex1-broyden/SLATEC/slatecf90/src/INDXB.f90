subroutine INDXB (I, IR, IDX, IDP)
!
!! INDXB is subsidiary to BLKTRI.
!
!***LIBRARY   SLATEC
!***TYPE      INTEGER (INDXB-I)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  BLKTRI
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    CBLKT
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   920422  Added statement so IDX would always be defined.  (WRB)
!***END PROLOGUE  INDXB
!
  COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        , &
                  NM         ,NCMPLX     ,IK
!***FIRST EXECUTABLE STATEMENT  INDXB
  IDX = I
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
