  INTEGER FUNCTION IPLOC (LOC, SX, IX)
!
!! IPLOC is subsidiary to SPLP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (IPLOC-S, IDLOC-D)
!***KEYWORDS  RELATIVE ADDRESS DETERMINATION FUNCTION, SLATEC
!***AUTHOR  Hanson, R. J., (SNLA)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   Given a "virtual" location,  IPLOC returns the relative working
!   address of the vector component stored in SX, IX.  Any necessary
!   page swaps are performed automatically for the user in this
!   function subprogram.
!
!   LOC       is the "virtual" address of the data to be retrieved.
!   SX ,IX    represent the matrix where the data is stored.
!
!***SEE ALSO  SPLP
!***ROUTINES CALLED  PRWPGE, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   810306  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890606  Restructured to match double precision version.  (WRB)
!   890606  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   910731  Added code to set IPLOC to 0 if LOC is non-positive.  (WRB)
!***END PROLOGUE  IPLOC
  REAL SX(*)
  INTEGER IX(*)
!***FIRST EXECUTABLE STATEMENT  IPLOC
  if (LOC <= 0) THEN
     call XERMSG ('SLATEC', 'IPLOC', &
       'A value of LOC, the first argument,  <=  0 was encountered', &
       55, 1)
     IPLOC = 0
     return
  end if
!
!     Two cases exist:  (1 <= LOC <= K) .OR. (LOC > K).
!
  K = IX(3) + 4
  LMX = IX(1)
  LMXM1 = LMX - 1
  if (LOC <= K) THEN
     IPLOC = LOC
     return
  end if
!
!     Compute length of the page, starting address of the page, page
!     number and relative working address.
!
  LPG = LMX-K
  ITEMP = LOC - K - 1
  IPAGE = ITEMP/LPG + 1
  IPLOC = MOD(ITEMP,LPG) + K + 1
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
        call PRWPGE (KEY, NP, LPG, SX, IX)
     ENDIF
     KEY = 1
     call PRWPGE (KEY, IPAGE, LPG, SX, IX)
  end if
  return
end
