subroutine C1MERG (TCOS, I1, M1, I2, M2, I3)
!
!! C1MERG merges two strings of complex numbers.  Each string is ...
!            ascending by the real part.
!
!***LIBRARY   SLATEC
!***TYPE      COMPLEX (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   This subroutine merges two ascending strings of numbers in the
!   array TCOS.  The first string is of length M1 and starts at
!   TCOS(I1+1).  The second string is of length M2 and starts at
!   TCOS(I2+1).  The merged string goes into TCOS(I3+1).  The ordering
!   is on the real part.
!
!***SEE ALSO  CMGNBN
!***ROUTINES CALLED  CCOPY
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   910408  Modified to use IF-THEN-ELSE.  Make it look like MERGE
!           which was modified earlier due to compiler problems on
!           the IBM RS6000.  (RWC)
!   920130  Code name changed from CMPMRG to C1MERG.  (WRB)
!***END PROLOGUE  C1MERG
  INTEGER I1, I2, I3, M1, M2
  COMPLEX TCOS(*)
!
  INTEGER J1, J2, J3
!
!***FIRST EXECUTABLE STATEMENT  C1MERG
  if (M1 == 0 .AND. M2 == 0) RETURN
!
  if (M1 == 0 .AND. M2 /= 0) THEN
     call CCOPY (M2, TCOS(I2+1), 1, TCOS(I3+1), 1)
     return
  end if
!
  if (M1 /= 0 .AND. M2 == 0) THEN
     call CCOPY (M1, TCOS(I1+1), 1, TCOS(I3+1), 1)
     return
  end if
!
  J1 = 1
  J2 = 1
  J3 = 1
!
   10 if (REAL(TCOS(J1+I1))  <=  REAL(TCOS(I2+J2))) THEN
     TCOS(I3+J3) = TCOS(I1+J1)
     J1 = J1+1
     if (J1  >  M1) THEN
        call CCOPY (M2-J2+1, TCOS(I2+J2), 1, TCOS(I3+J3+1), 1)
        return
     ENDIF
  ELSE
     TCOS(I3+J3) = TCOS(I2+J2)
     J2 = J2+1
     if (J2  >  M2) THEN
        call CCOPY (M1-J1+1, TCOS(I1+J1), 1, TCOS(I3+J3+1), 1)
        return
     ENDIF
  end if
  J3 = J3+1
  go to 10
end
