subroutine DWNLT1 (I, LEND, MEND, IR, MDW, RECALC, IMAX, HBAR, H, &
     SCALE, W)
!
!! DWNLT1 is subsidiary to WNLIT.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (WNLT1-S, DWNLT1-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     To update the column Sum Of Squares and find the pivot column.
!     The column Sum of Squares Vector will be updated at each step.
!     When numerically necessary, these values will be recomputed.
!
!***SEE ALSO  DWNLIT
!***ROUTINES CALLED  IDAMAX
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
!   900604  DP version created from SP version.  (RWC)
!***END PROLOGUE  DWNLT1
  INTEGER I, IMAX, IR, LEND, MDW, MEND
  DOUBLE PRECISION H(*), HBAR, SCALE(*), W(MDW,*)
  LOGICAL RECALC
!
  EXTERNAL IDAMAX
  INTEGER IDAMAX
!
  INTEGER J, K
!
!***FIRST EXECUTABLE STATEMENT  DWNLT1
  if (IR /= 1 .AND. (.NOT.RECALC)) THEN
!
!        Update column SS=sum of squares.
!
     DO 10 J=I,LEND
        H(J) = H(J) - SCALE(IR-1)*W(IR-1,J)**2
   10    CONTINUE
!
!        Test for numerical accuracy.
!
     IMAX = IDAMAX(LEND-I+1, H(I), 1) + I - 1
     RECALC = (HBAR+1.E-3*H(IMAX))  ==  HBAR
  end if
!
!     If required, recalculate column SS, using rows IR through MEND.
!
  if (RECALC) THEN
     DO 30 J=I,LEND
        H(J) = 0.D0
        DO 20 K=IR,MEND
           H(J) = H(J) + SCALE(K)*W(K,J)**2
   20       CONTINUE
   30    CONTINUE
!
!        Find column with largest SS.
!
     IMAX = IDAMAX(LEND-I+1, H(I), 1) + I - 1
     HBAR = H(IMAX)
  end if
  return
end
