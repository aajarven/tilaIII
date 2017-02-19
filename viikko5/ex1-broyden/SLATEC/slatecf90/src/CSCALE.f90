subroutine CSCALE (A, NRDA, NROW, NCOL, COLS, COLSAV, ROWS, &
     ROWSAV, ANORM, SCALES, ISCALE, IC)
!
!! CSCALE is subsidiary to BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (CSCALE-S, DCSCAL-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     This routine scales the matrix A by columns when needed
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  SDOT
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  CSCALE
  DIMENSION A(NRDA,*),COLS(*),COLSAV(*),SCALES(*), &
            ROWS(*),ROWSAV(*)
!
  SAVE TEN4, TEN20
  DATA TEN4,TEN20/1.E+4,1.E+20/
!
!***FIRST EXECUTABLE STATEMENT  CSCALE
  if (ISCALE  /=  (-1)) go to 25
!
  if (IC  ==  0) go to 10
  DO 5 K=1,NCOL
    5    COLS(K)=SDOT(NROW,A(1,K),1,A(1,K),1)
!
   10 ASCALE=ANORM/NCOL
  DO 20 K=1,NCOL
     CS=COLS(K)
     if ((CS  >  TEN4*ASCALE) .OR. (TEN4*CS  <  ASCALE)) go to 50
     if ((CS  <  1./TEN20) .OR. (CS  >  TEN20)) go to 50
   20 CONTINUE
!
   25 DO 30 K=1,NCOL
   30    SCALES(K)=1.
  return
!
   50 ALOG2=LOG(2.)
  ANORM=0.
  DO 100 K=1,NCOL
     CS=COLS(K)
     if (CS  /=  0.) go to 60
     SCALES(K)=1.
     go to 100
   60    P=LOG(CS)/ALOG2
     IP=-0.5*P
     S=2.**IP
     SCALES(K)=S
     if (IC  ==  1) go to 70
     COLS(K)=S*S*COLS(K)
     ANORM=ANORM+COLS(K)
     COLSAV(K)=COLS(K)
   70    DO 80 J=1,NROW
   80       A(J,K)=S*A(J,K)
  100 CONTINUE
!
  if (IC  ==  0) RETURN
!
  DO 200 K=1,NROW
     ROWS(K)=SDOT(NCOL,A(K,1),NRDA,A(K,1),NRDA)
     ROWSAV(K)=ROWS(K)
  200    ANORM=ANORM+ROWS(K)
  return
end
