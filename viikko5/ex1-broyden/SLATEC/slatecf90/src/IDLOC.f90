  INTEGER FUNCTION IDLOC (LOC, SX, IX)
!
!! IDLOC is subsidiary to DSPLP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (IPLOC-S, IDLOC-D)
!***KEYWORDS  RELATIVE ADDRESS DETERMINATION FUNCTION, SLATEC
!***AUTHOR  Boland, W. Robert, (LANL)
!           Nicol, Tom, (University of British Columbia)
!***DESCRIPTION
!
!   Given a "virtual" location,  IDLOC returns the relative working
!   address of the vector component stored in SX, IX.  Any necessary
!   page swaps are performed automatically for the user in this
!   function subprogram.
!
!   LOC       is the "virtual" address of the data to be retrieved.
!   SX ,IX    represent the matrix where the data is stored.
!
!***SEE ALSO  DSPLP
!***ROUTINES CALLED  DPRWPG, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   890606  DATE WRITTEN
!   890606  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   910731  Added code to set IDLOC to 0 if LOC is non-positive.  (WRB)
!***END PROLOGUE  IDLOC
  DOUBLE PRECISION SX(*)
  INTEGER IX(*)
!***FIRST EXECUTABLE STATEMENT  IDLOC
  if (LOC <= 0) THEN
     call XERMSG ('SLATEC', 'IDLOC', &
       'A value of LOC, the first argument,  <=  0 was encountered', &
       55, 1)
     IDLOC = 0
     return
  end if
!
!     Two cases exist:  (1 <= LOC <= K) .OR. (LOC > K).
!
  K = IX(3) + 4
  LMX = IX(1)
  LMXM1 = LMX - 1
  if (LOC <= K) THEN
     IDLOC = LOC
     return
  end if
!
!     Compute length of the page, starting address of the page, page
!     number and relative working address.
!
  LPG = LMX-K
  ITEMP = LOC - K - 1
  IPAGE = ITEMP/LPG + 1
  IDLOC = MOD(ITEMP,LPG) + K + 1
  NP = ABS(IX(LMXM1))
!
!     Determine if a page fault has occurred.  If so, write page NP
!     and read page IPAGE.  Write the page only if it has been
!     modified.
!
  if (IPAGE /= NP) THEN
     if (SX(LMX) == 1.0) THEN
        SX(LMX) = 0.0
        KEY = 2
        call DPRWPG (KEY, NP, LPG, SX, IX)
     ENDIF
     KEY = 1
     call DPRWPG (KEY, IPAGE, LPG, SX, IX)
  end if
  return
end
