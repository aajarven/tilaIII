subroutine DCSCAL (A, NRDA, NROW, NCOL, COLS, COLSAV, ROWS, &
     ROWSAV, ANORM, SCALES, ISCALE, IC)
!
!! DCSCAL is subsidiary to DBVSUP and DSUDS.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (CSCALE-S, DCSCAL-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     This routine scales the matrix A by columns when needed.
!
!***SEE ALSO  DBVSUP, DSUDS
!***ROUTINES CALLED  DDOT
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DCSCAL
  DOUBLE PRECISION DDOT
  INTEGER IC, IP, ISCALE, J, K, NCOL, NRDA, NROW
  DOUBLE PRECISION A(NRDA,*), ALOG2, ANORM, ASCALE, COLS(*), &
       COLSAV(*), CS, P, ROWS(*), ROWSAV(*), S, &
       SCALES(*), TEN20, TEN4
!
  SAVE TEN4, TEN20
  DATA TEN4,TEN20 /1.0D4,1.0D20/
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 130
!        BEGIN BLOCK PERMITTING ...EXITS TO 60
!***FIRST EXECUTABLE STATEMENT  DCSCAL
        if (ISCALE  /=  (-1)) go to 40
!
           if (IC  ==  0) go to 20
              DO 10 K = 1, NCOL
                 COLS(K) = DDOT(NROW,A(1,K),1,A(1,K),1)
   10             CONTINUE
   20          CONTINUE
!
           ASCALE = ANORM/NCOL
           DO 30 K = 1, NCOL
              CS = COLS(K)
!        .........EXIT
              if ((CS  >  TEN4*ASCALE) .OR. (TEN4*CS  <  ASCALE)) &
                 go to 60
!        .........EXIT
              if ((CS  <  1.0D0/TEN20) .OR. (CS  >  TEN20)) &
                 go to 60
   30          CONTINUE
   40       CONTINUE
!
        DO 50 K = 1, NCOL
           SCALES(K) = 1.0D0
   50       CONTINUE
!     ......EXIT
        go to 130
   60    CONTINUE
!
     ALOG2 = LOG(2.0D0)
     ANORM = 0.0D0
     DO 110 K = 1, NCOL
        CS = COLS(K)
        if (CS  /=  0.0D0) go to 70
           SCALES(K) = 1.0D0
        go to 100
   70       CONTINUE
           P = LOG(CS)/ALOG2
           IP = -0.5D0*P
           S = 2.0D0**IP
           SCALES(K) = S
           if (IC  ==  1) go to 80
              COLS(K) = S*S*COLS(K)
              ANORM = ANORM + COLS(K)
              COLSAV(K) = COLS(K)
   80          CONTINUE
           DO 90 J = 1, NROW
              A(J,K) = S*A(J,K)
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
!
!     ...EXIT
     if (IC  ==  0) go to 130
!
     DO 120 K = 1, NROW
        ROWS(K) = DDOT(NCOL,A(K,1),NRDA,A(K,1),NRDA)
        ROWSAV(K) = ROWS(K)
        ANORM = ANORM + ROWS(K)
  120    CONTINUE
  130 CONTINUE
  return
end
