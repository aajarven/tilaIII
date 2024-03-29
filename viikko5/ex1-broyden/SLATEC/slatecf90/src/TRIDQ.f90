subroutine TRIDQ (MR, A, B, C, Y, D)
!
!! TRIDQ is subsidiary to POIS3D.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to POIS3D
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (TRIDQ-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  POIS3D
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900308  Renamed routine from TRID to TRIDQ.  (WRB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  TRIDQ
  DIMENSION       A(*)       ,B(*)       ,C(*)       ,Y(*)       , &
                  D(*)
!***FIRST EXECUTABLE STATEMENT  TRIDQ
  M = MR
  MM1 = M-1
  Z = 1./B(1)
  D(1) = C(1)*Z
  Y(1) = Y(1)*Z
  DO 101 I=2,MM1
     Z = 1./(B(I)-A(I)*D(I-1))
     D(I) = C(I)*Z
     Y(I) = (Y(I)-A(I)*Y(I-1))*Z
  101 CONTINUE
  Z = B(M)-A(M)*D(MM1)
  if (Z  /=  0.) go to 102
  Y(M) = 0.
  go to 103
  102 Y(M) = (Y(M)-A(M)*Y(MM1))/Z
  103 CONTINUE
  DO 104 IP=1,MM1
     I = M-IP
     Y(I) = Y(I)-D(I)*Y(I+1)
  104 CONTINUE
  return
end
