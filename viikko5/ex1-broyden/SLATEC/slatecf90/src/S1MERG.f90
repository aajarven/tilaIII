subroutine S1MERG (TCOS, I1, M1, I2, M2, I3)
!
!! S1MERG merges two strings of ascending real numbers.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   This subroutine merges two ascending strings of numbers in the
!   array TCOS.  The first string is of length M1 and starts at
!   TCOS(I1+1).  The second string is of length M2 and starts at
!   TCOS(I2+1).  The merged string goes into TCOS(I3+1).
!
!***SEE ALSO  GENBUN
!***ROUTINES CALLED  SCOPY
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   901120  Modified to use IF-THEN-ELSE.  Previous spaghetti code did
!           not compile correctly with optimization on the IBM RS6000.
!           (RWC)
!   920130  Code name changed from MERGE to S1MERG.  (WRB)
!***END PROLOGUE  S1MERG
  INTEGER I1, I2, I3, M1, M2
  REAL TCOS(*)
!
  INTEGER J1, J2, J3
!
!***FIRST EXECUTABLE STATEMENT  S1MERG
  if (M1 == 0 .AND. M2 == 0) RETURN
!
  if (M1 == 0 .AND. M2 /= 0) THEN
     call SCOPY (M2, TCOS(I2+1), 1, TCOS(I3+1), 1)
     return
  end if
!
  if (M1 /= 0 .AND. M2 == 0) THEN
     call SCOPY (M1, TCOS(I1+1), 1, TCOS(I3+1), 1)
     return
  end if
!
  J1 = 1
  J2 = 1
  J3 = 1
!
   10 if (TCOS(I1+J1)  <=  TCOS(I2+J2)) THEN
     TCOS(I3+J3) = TCOS(I1+J1)
     J1 = J1+1
     if (J1  >  M1) THEN
        call SCOPY (M2-J2+1, TCOS(I2+J2), 1, TCOS(I3+J3+1), 1)
        return
     ENDIF
  ELSE
     TCOS(I3+J3) = TCOS(I2+J2)
     J2 = J2+1
     if (J2  >  M2) THEN
        call SCOPY (M1-J1+1, TCOS(I1+J1), 1, TCOS(I3+J3+1), 1)
        return
     ENDIF
  end if
  J3 = J3+1
  go to 10
end
