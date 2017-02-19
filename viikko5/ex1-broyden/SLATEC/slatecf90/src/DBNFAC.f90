subroutine DBNFAC (W, NROWW, NROW, NBANDL, NBANDU, IFLAG)
!
!! DBNFAC is subsidiary to DBINT4 and DBINTK.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (BNFAC-S, DBNFAC-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!  DBNFAC is the BANFAC routine from
!        * A Practical Guide to Splines *  by C. de Boor
!
!  DBNFAC is a double precision routine
!
!  Returns in  W  the LU-factorization (without pivoting) of the banded
!  matrix  A  of order  NROW  with  (NBANDL + 1 + NBANDU) bands or diag-
!  onals in the work array  W .
!
! *****  I N P U T  ****** W is double precision
!  W.....Work array of size  (NROWW,NROW)  containing the interesting
!        part of a banded matrix  A , with the diagonals or bands of  A
!        stored in the rows of  W , while columns of  A  correspond to
!        columns of  W . This is the storage mode used in  LINPACK  and
!        results in efficient innermost loops.
!           Explicitly,  A  has  NBANDL  bands below the diagonal
!                            +     1     (main) diagonal
!                            +   NBANDU  bands above the diagonal
!        and thus, with    MIDDLE = NBANDU + 1,
!          A(I+J,J)  is in  W(I+MIDDLE,J)  for I=-NBANDU,...,NBANDL
!                                              J=1,...,NROW .
!        For example, the interesting entries of A (1,2)-banded matrix
!        of order  9  would appear in the first  1+1+2 = 4  rows of  W
!        as follows.
!                          13243546576879
!                       1223344556677889
!                    112233445566778899
!                    2132435465768798
!
!        All other entries of  W  not identified in this way with an en-
!        try of  A  are never referenced .
!  NROWW.....Row dimension of the work array  W .
!        must be   >=   NBANDL + 1 + NBANDU  .
!  NBANDL.....Number of bands of  A  below the main diagonal
!  NBANDU.....Number of bands of  A  above the main diagonal .
!
! *****  O U T P U T  ****** W is double precision
!  IFLAG.....Integer indicating success( = 1) or failure ( = 2) .
!     If  IFLAG = 1, then
!  W.....contains the LU-factorization of  A  into a unit lower triangu-
!        lar matrix  L  and an upper triangular matrix  U (both banded)
!        and stored in customary fashion over the corresponding entries
!        of  A . This makes it possible to solve any particular linear
!        system  A*X = B  for  X  by a
!              call DBNSLV ( W, NROWW, NROW, NBANDL, NBANDU, B )
!        with the solution X  contained in  B  on return .
!     If  IFLAG = 2, then
!        one of  NROW-1, NBANDL,NBANDU failed to be nonnegative, or else
!        one of the potential pivots was found to be zero indicating
!        that  A  does not have an LU-factorization. This implies that
!        A  is singular in case it is totally positive .
!
! *****  M E T H O D  ******
!     Gauss elimination  W I T H O U T  pivoting is used. The routine is
!  intended for use with matrices  A  which do not require row inter-
!  changes during factorization, especially for the  T O T A L L Y
!  P O S I T I V E  matrices which occur in spline calculations.
!     The routine should NOT be used for an arbitrary banded matrix.
!
!***SEE ALSO  DBINT4, DBINTK
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DBNFAC
!
  INTEGER IFLAG, NBANDL, NBANDU, NROW, NROWW, I, IPK, J, JMAX, K, &
   KMAX, MIDDLE, MIDMK, NROWM1
  DOUBLE PRECISION W(NROWW,*), FACTOR, PIVOT
!
!***FIRST EXECUTABLE STATEMENT  DBNFAC
  IFLAG = 1
  MIDDLE = NBANDU + 1
!                         W(MIDDLE,.) CONTAINS THE MAIN DIAGONAL OF  A .
  NROWM1 = NROW - 1
  if (NROWM1) 120, 110, 10
   10 if (NBANDL > 0) go to 30
!                A IS UPPER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO .
  DO 20 I=1,NROWM1
    if (W(MIDDLE,I) == 0.0D0) go to 120
   20 CONTINUE
  go to 110
   30 if (NBANDU > 0) go to 60
!              A IS LOWER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO AND
!                 DIVIDE EACH COLUMN BY ITS DIAGONAL .
  DO 50 I=1,NROWM1
    PIVOT = W(MIDDLE,I)
    if (PIVOT == 0.0D0) go to 120
    JMAX = MIN(NBANDL,NROW-I)
    DO 40 J=1,JMAX
      W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
   40   CONTINUE
   50 CONTINUE
  return
!
!        A  IS NOT JUST A TRIANGULAR MATRIX. CONSTRUCT LU FACTORIZATION
   60 DO 100 I=1,NROWM1
!                                  W(MIDDLE,I)  IS PIVOT FOR I-TH STEP .
    PIVOT = W(MIDDLE,I)
    if (PIVOT == 0.0D0) go to 120
!                 JMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN COLUMN  I
!                     BELOW THE DIAGONAL .
    JMAX = MIN(NBANDL,NROW-I)
!              DIVIDE EACH ENTRY IN COLUMN  I  BELOW DIAGONAL BY PIVOT .
    DO 70 J=1,JMAX
      W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
   70   CONTINUE
!                 KMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN ROW  I  TO
!                     THE RIGHT OF THE DIAGONAL .
    KMAX = MIN(NBANDU,NROW-I)
!                  SUBTRACT  A(I,I+K)*(I-TH COLUMN) FROM (I+K)-TH COLUMN
!                  (BELOW ROW  I ) .
    DO 90 K=1,KMAX
      IPK = I + K
      MIDMK = MIDDLE - K
      FACTOR = W(MIDMK,IPK)
      DO 80 J=1,JMAX
        W(MIDMK+J,IPK) = W(MIDMK+J,IPK) - W(MIDDLE+J,I)*FACTOR
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
!                                       CHECK THE LAST DIAGONAL ENTRY .
  110 if (W(MIDDLE,NROW) /= 0.0D0) RETURN
  120 IFLAG = 2
  return
end
