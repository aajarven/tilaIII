subroutine BSPLVN (T, JHIGH, INDEX, X, ILEFT, VNIKX)
!
!! BSPLVN is subsidiary to FC.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (BSPLVN-S, DFSPVN-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
! Calculates the value of all possibly nonzero B-splines at *X* of
!  order MAX(JHIGH,(J+1)(INDEX-1)) on *T*.
!
!***SEE ALSO  FC
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   780801  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  BSPLVN
  DIMENSION T(*),VNIKX(*)
  DIMENSION DELTAM(20),DELTAP(20)
  SAVE J, DELTAM, DELTAP
  DATA J/1/,(DELTAM(I),I=1,20),(DELTAP(I),I=1,20)/40*0./
!***FIRST EXECUTABLE STATEMENT  BSPLVN
                                   go to (10,20),INDEX
   10 J = 1
  VNIKX(1) = 1.
  if (J  >=  JHIGH)                go to 99
!
   20    IPJ = ILEFT+J
     DELTAP(J) = T(IPJ) - X
     IMJP1 = ILEFT-J+1
     DELTAM(J) = X - T(IMJP1)
     VMPREV = 0.
     JP1 = J+1
     DO 26 L=1,J
        JP1ML = JP1-L
        VM = VNIKX(L)/(DELTAP(L) + DELTAM(JP1ML))
        VNIKX(L) = VM*DELTAP(L) + VMPREV
   26       VMPREV = VM*DELTAM(JP1ML)
     VNIKX(JP1) = VMPREV
     J = JP1
     if (J  <  JHIGH)             go to 20
!
   99                                  return
end
