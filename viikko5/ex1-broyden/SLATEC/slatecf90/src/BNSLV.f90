subroutine BNSLV (W, NROWW, NROW, NBANDL, NBANDU, B)
!
!! BNSLV is subsidiary to BINT4 and BINTK.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (BNSLV-S, DBNSLV-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!  BNSLV is the BANSLV routine from
!        * A Practical Guide to Splines *  by C. de Boor
!
!  Companion routine to  BNFAC . It returns the solution  X  of the
!  linear system  A*X = B  in place of  B , given the LU-factorization
!  for  A  in the work array  W from BNFAC.
!
! *****  I N P U T  ******
!  W, NROWW,NROW,NBANDL,NBANDU.....Describe the LU-factorization of a
!        banded matrix  A  of order  NROW  as constructed in  BNFAC .
!        For details, see  BNFAC .
!  B.....Right side of the system to be solved .
!
! *****  O U T P U T  ******
!  B.....Contains the solution  X , of order  NROW .
!
! *****  M E T H O D  ******
!     (With  A = L*U, as stored in  W,) the unit lower triangular system
!  L(U*X) = B  is solved for  Y = U*X, and  Y  stored in  B . Then the
!  upper triangular system  U*X = Y  is solved for  X  . The calcul-
!  ations are so arranged that the innermost loops stay within columns.
!
!***SEE ALSO  BINT4, BINTK
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  BNSLV
!
  INTEGER NBANDL, NBANDU, NROW, NROWW, I, J, JMAX, MIDDLE, NROWM1
  REAL W(NROWW,*), B(*)
!***FIRST EXECUTABLE STATEMENT  BNSLV
  MIDDLE = NBANDU + 1
  if (NROW == 1) go to 80
  NROWM1 = NROW - 1
  if (NBANDL == 0) go to 30
!                                 FORWARD PASS
!            FOR I=1,2,...,NROW-1, SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
!            OF  L )  FROM RIGHT SIDE  (BELOW I-TH ROW) .
  DO 20 I=1,NROWM1
    JMAX = MIN(NBANDL,NROW-I)
    DO 10 J=1,JMAX
      B(I+J) = B(I+J) - B(I)*W(MIDDLE+J,I)
   10   CONTINUE
   20 CONTINUE
!                                 BACKWARD PASS
!            FOR I=NROW,NROW-1,...,1, DIVIDE RIGHT SIDE(I) BY I-TH DIAG-
!            ONAL ENTRY OF  U, THEN SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
!            OF  U)  FROM RIGHT SIDE  (ABOVE I-TH ROW).
   30 if (NBANDU > 0) go to 50
!                                A  IS LOWER TRIANGULAR .
  DO 40 I=1,NROW
    B(I) = B(I)/W(1,I)
   40 CONTINUE
  return
   50 I = NROW
   60 B(I) = B(I)/W(MIDDLE,I)
  JMAX = MIN(NBANDU,I-1)
  DO 70 J=1,JMAX
    B(I-J) = B(I-J) - B(I)*W(MIDDLE-J,I)
   70 CONTINUE
  I = I - 1
  if (I > 1) go to 60
   80 B(1) = B(1)/W(MIDDLE,1)
  return
end
