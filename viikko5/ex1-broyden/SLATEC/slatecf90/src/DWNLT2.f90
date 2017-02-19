  LOGICAL FUNCTION DWNLT2 (ME, MEND, IR, FACTOR, TAU, SCALE, WIC)
!
!! DWNLT2 is subsidiary to WNLIT.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (WNLT2-S, DWNLT2-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     To test independence of incoming column.
!
!     Test the column IC to determine if it is linearly independent
!     of the columns already in the basis.  In the initial tri. step,
!     we usually want the heavy weight ALAMDA to be included in the
!     test for independence.  In this case, the value of FACTOR will
!     have been set to 1.E0 before this procedure is invoked.
!     In the potentially rank deficient problem, the value of FACTOR
!     will have been set to ALSQ=ALAMDA**2 to remove the effect of the
!     heavy weight from the test for independence.
!
!     Write new column as partitioned vector
!           (A1)  number of components in solution so far = NIV
!           (A2)  M-NIV components
!     And compute  SN = inverse weighted length of A1
!                  RN = inverse weighted length of A2
!     Call the column independent when RN  >  TAU*SN
!
!***SEE ALSO  DWNLIT
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
!   900604  DP version created from SP version.  (RWC)
!***END PROLOGUE  DWNLT2
  DOUBLE PRECISION FACTOR, SCALE(*), TAU, WIC(*)
  INTEGER IR, ME, MEND
!
  DOUBLE PRECISION RN, SN, T
  INTEGER J
!
!***FIRST EXECUTABLE STATEMENT  DWNLT2
  SN = 0.E0
  RN = 0.E0
  DO 10 J=1,MEND
     T = SCALE(J)
     if (J <= ME) T = T/FACTOR
     T = T*WIC(J)**2
!
     if (J < IR) THEN
        SN = SN + T
     ELSE
        RN = RN + T
     ENDIF
   10 CONTINUE
  DWNLT2 = RN  >  SN*TAU**2
  return
end
