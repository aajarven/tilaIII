subroutine WNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
!
!! WNLT3 is subsidiary to WNLIT.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (WNLT3-S, DWNLT3-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     Perform column interchange.
!     Exchange elements of permuted index vector and perform column
!     interchanges.
!
!***SEE ALSO  WNLIT
!***ROUTINES CALLED  SSWAP
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890620  Code extracted from WNLT and made a subroutine.  (RWC))
!***END PROLOGUE  WNLT3
  INTEGER I, IMAX, IPIVOT(*), M, MDW
  REAL             H(*), W(MDW,*)
!
  EXTERNAL SSWAP
!
  REAL             T
  INTEGER ITEMP
!
!***FIRST EXECUTABLE STATEMENT  WNLT3
  if (IMAX /= I) THEN
     ITEMP        = IPIVOT(I)
     IPIVOT(I)    = IPIVOT(IMAX)
     IPIVOT(IMAX) = ITEMP
!
     call SSWAP(M, W(1,IMAX), 1, W(1,I), 1)
!
     T       = H(IMAX)
     H(IMAX) = H(I)
     H(I)    = T
  end if
  return
end
