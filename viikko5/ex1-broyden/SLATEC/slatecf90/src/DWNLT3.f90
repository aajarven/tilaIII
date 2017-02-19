subroutine DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
!
!! DWNLT3 is subsidiary to WNLIT.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (WNLT3-S, DWNLT3-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     Perform column interchange.
!     Exchange elements of permuted index vector and perform column
!     interchanges.
!
!***SEE ALSO  DWNLIT
!***ROUTINES CALLED  DSWAP
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
!   900604  DP version created from SP version.  (RWC)
!***END PROLOGUE  DWNLT3
  INTEGER I, IMAX, IPIVOT(*), M, MDW
  DOUBLE PRECISION H(*), W(MDW,*)
!
  EXTERNAL DSWAP
!
  DOUBLE PRECISION T
  INTEGER ITEMP
!
!***FIRST EXECUTABLE STATEMENT  DWNLT3
  if (IMAX /= I) THEN
     ITEMP        = IPIVOT(I)
     IPIVOT(I)    = IPIVOT(IMAX)
     IPIVOT(IMAX) = ITEMP
!
     call DSWAP(M, W(1,IMAX), 1, W(1,I), 1)
!
     T       = H(IMAX)
     H(IMAX) = H(I)
     H(I)    = T
  end if
  return
end
