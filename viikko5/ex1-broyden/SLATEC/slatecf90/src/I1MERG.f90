subroutine I1MERG (ICOS, I1, M1, I2, M2, I3)
!
!! I1MERG merges two strings of ascending integers.
!
!***LIBRARY   SLATEC
!***TYPE      INTEGER (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I)
!***AUTHOR  Boland, W. Robert, (LANL)
!           Clemens, Reginald, (PLK)
!***DESCRIPTION
!
!   This subroutine merges two ascending strings of integers in the
!   array ICOS.  The first string is of length M1 and starts at
!   ICOS(I1+1).  The second string is of length M2 and starts at
!   ICOS(I2+1).  The merged string goes into ICOS(I3+1).
!
!***ROUTINES CALLED  ICOPY
!***REVISION HISTORY  (YYMMDD)
!   920202  DATE WRITTEN
!***END PROLOGUE  I1MERG
  INTEGER I1, I2, I3, M1, M2
  REAL ICOS(*)
!
  INTEGER J1, J2, J3
!
!***FIRST EXECUTABLE STATEMENT  I1MERG
  if (M1 == 0 .AND. M2 == 0) RETURN
!
  if (M1 == 0 .AND. M2 /= 0) THEN
     call ICOPY (M2, ICOS(I2+1), 1, ICOS(I3+1), 1)
     return
  end if
!
  if (M1 /= 0 .AND. M2 == 0) THEN
     call ICOPY (M1, ICOS(I1+1), 1, ICOS(I3+1), 1)
     return
  end if
!
  J1 = 1
  J2 = 1
  J3 = 1
!
   10 if (ICOS(I1+J1)  <=  ICOS(I2+J2)) THEN
     ICOS(I3+J3) = ICOS(I1+J1)
     J1 = J1+1
     if (J1  >  M1) THEN
        call ICOPY (M2-J2+1, ICOS(I2+J2), 1, ICOS(I3+J3+1), 1)
        return
     ENDIF
  ELSE
     ICOS(I3+J3) = ICOS(I2+J2)
     J2 = J2+1
     if (J2  >  M2) THEN
        call ICOPY (M1-J1+1, ICOS(I1+J1), 1, ICOS(I3+J3+1), 1)
        return
     ENDIF
  end if
  J3 = J3+1
  go to 10
end
