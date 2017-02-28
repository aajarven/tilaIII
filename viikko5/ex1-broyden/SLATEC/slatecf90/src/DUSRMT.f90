subroutine DUSRMT (I, J, AIJ, INDCAT, PRGOPT, DATTRV, IFLAG)
!
!! DUSRMT is subsidiary to DSPLP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (USRMAT-S, DUSRMT-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   The user may supply this code
!
!***SEE ALSO  DSPLP
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   811215  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DUSRMT
  DOUBLE PRECISION PRGOPT(*),DATTRV(*),AIJ
  INTEGER IFLAG(*)
!
!***FIRST EXECUTABLE STATEMENT  DUSRMT
  if ( IFLAG(1) == 1) THEN
!
!     THIS IS THE INITIALIZATION STEP.  THE VALUES OF IFLAG(K),K=2,3,4,
!     ARE RESPECTIVELY THE COLUMN INDEX, THE ROW INDEX (OR THE NEXT COL.
!     INDEX), AND THE POINTER TO THE MATRIX ENTRY'S VALUE WITHIN
!     DATTRV(*).  ALSO CHECK (DATTRV(1)=0.) SIGNIFYING NO DATA.
       if ( DATTRV(1) == 0.D0) THEN
       I = 0
       J = 0
       IFLAG(1) = 3
       ELSE
       IFLAG(2)=-DATTRV(1)
       IFLAG(3)= DATTRV(2)
       IFLAG(4)= 3
       ENDIF
!
       return
  ELSE
       J=IFLAG(2)
       I=IFLAG(3)
       L=IFLAG(4)
       if ( I == 0) THEN
!
!     SIGNAL THAT ALL OF THE NONZERO ENTRIES HAVE BEEN DEFINED.
            IFLAG(1)=3
            return
       ELSE if ( I < 0) THEN
!
!     SIGNAL THAT A SWITCH IS MADE TO A NEW COLUMN.
            J=-I
            I=DATTRV(L)
            L=L+1
       ENDIF
!
       AIJ=DATTRV(L)
!
!     UPDATE THE INDICES AND POINTERS FOR THE NEXT ENTRY.
       IFLAG(2)=J
       IFLAG(3)=DATTRV(L+1)
       IFLAG(4)=L+2
!
!     INDCAT=0 DENOTES THAT ENTRIES OF THE MATRIX ARE ASSIGNED THE
!     VALUES FROM DATTRV(*).  NO ACCUMULATION IS PERFORMED.
       INDCAT=0
       return
  end if
end