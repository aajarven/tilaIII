subroutine D1MERG (TCOS, I1, M1, I2, M2, I3)
!
!! D1MERG merges two strings of ascending double precision numbers.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (S1MERG-S, D1MERG-D, CMERGE-C, I1MERG-I)
!***AUTHOR  Boland, W. Robert, (LANL)
!           Clemens, Reginald, (PLK)
!***DESCRIPTION
!
!   This subroutine merges two ascending strings of numbers in the
!   array TCOS.  The first string is of length M1 and starts at
!   TCOS(I1+1).  The second string is of length M2 and starts at
!   TCOS(I2+1).  The merged string goes into TCOS(I3+1).
!
!   This routine is currently unused, but was added to complete
!   the set of routines S1MERG and C1MERG (both of which are used).
!
!***ROUTINES CALLED  DCOPY
!***REVISION HISTORY  (YYMMDD)
!   910819  DATE WRITTEN
!***END PROLOGUE  D1MERG
  INTEGER I1, I2, I3, M1, M2
  DOUBLE PRECISION TCOS(*)
!
  INTEGER J1, J2, J3
!
!***FIRST EXECUTABLE STATEMENT  D1MERG
  if (M1 == 0 .AND. M2 == 0) RETURN
!
  if (M1 == 0 .AND. M2 /= 0) THEN
     call DCOPY (M2, TCOS(I2+1), 1, TCOS(I3+1), 1)
     return
  end if
!
  if (M1 /= 0 .AND. M2 == 0) THEN
     call DCOPY (M1, TCOS(I1+1), 1, TCOS(I3+1), 1)
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
        call DCOPY (M2-J2+1, TCOS(I2+J2), 1, TCOS(I3+J3+1), 1)
        return
     ENDIF
  ELSE
     TCOS(I3+J3) = TCOS(I2+J2)
     J2 = J2+1
     if (J2  >  M2) THEN
        call DCOPY (M1-J1+1, TCOS(I1+J1), 1, TCOS(I3+J3+1), 1)
        return
     ENDIF
  end if
  J3 = J3+1
  go to 10
end
