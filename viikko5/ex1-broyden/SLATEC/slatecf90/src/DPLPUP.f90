subroutine DPLPUP (DUSRMT, MRELAS, NVARS, PRGOPT, DATTRV, BL, BU, &
     IND, INFO, AMAT, IMAT, SIZEUP, ASMALL, ABIG)
!
!! DPLPUP is subsidiary to DSPLP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SPLPUP-S, DPLPUP-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
!     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
!
!     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
!     /REAL (12 BLANKS)/DOUBLE PRECISION/.
!
!     REVISED 810613-1130
!     REVISED YYMMDD-HHMM
!
!     THIS SUBROUTINE COLLECTS INFORMATION ABOUT THE BOUNDS AND MATRIX
!     FROM THE USER.  IT IS PART OF THE DSPLP( ) PACKAGE.
!
!***SEE ALSO  DSPLP
!***ROUTINES CALLED  DPCHNG, DPNNZR, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811215  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890605  Corrected references to XERRWV.  (WRB)
!   890605  Removed unreferenced labels.  (WRB)
!   891009  Removed unreferenced variables.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls, changed do-it-yourself
!           DO loops to DO loops.  (RWC)
!   900602  Get rid of ASSIGNed GOTOs.  (RWC)
!***END PROLOGUE  DPLPUP
  DOUBLE PRECISION ABIG,AIJ,AMAT(*),AMN,AMX,ASMALL,BL(*), &
   BU(*),DATTRV(*),PRGOPT(*),XVAL,ZERO
  INTEGER IFLAG(10),IMAT(*),IND(*)
  LOGICAL SIZEUP,FIRST
  CHARACTER*8 XERN1, XERN2
  CHARACTER*16 XERN3, XERN4
!
!***FIRST EXECUTABLE STATEMENT  DPLPUP
  ZERO = 0.D0
!
!     CHECK USER-SUPPLIED BOUNDS
!
!     CHECK THAT IND(*) VALUES ARE 1,2,3 OR 4.
!     ALSO CHECK CONSISTENCY OF UPPER AND LOWER BOUNDS.
!
  DO 10 J=1,NVARS
     if (IND(J) < 1 .OR. IND(J) > 4) THEN
        WRITE (XERN1, '(I8)') J
        call XERMSG ('SLATEC', 'DPLPUP', &
           'IN DSPLP, INDEPENDENT VARIABLE = ' // XERN1 // &
           ' IS NOT DEFINED.', 10, 1)
        INFO = -10
        return
     ENDIF
!
     if (IND(J) == 3) THEN
        if (BL(J) > BU(J)) THEN
           WRITE (XERN1, '(I8)') J
           WRITE (XERN3, '(1PE15.6)') BL(J)
           WRITE (XERN4, '(1PE15.6)') BU(J)
           call XERMSG ('SLATEC', 'DPLPUP', &
              'IN DSPLP, LOWER BOUND = ' // XERN3 // &
              ' AND UPPER BOUND = ' // XERN4 // &
              ' FOR INDEPENDENT VARIABLE = ' // XERN1 // &
              ' ARE NOT CONSISTENT.', 11, 1)
           return
        ENDIF
     ENDIF
   10 CONTINUE
!
  DO 20 I=NVARS+1,NVARS+MRELAS
     if (IND(I) < 1 .OR. IND(I) > 4) THEN
        WRITE (XERN1, '(I8)') I-NVARS
        call XERMSG ('SLATEC', 'DPLPUP', &
           'IN DSPLP, DEPENDENT VARIABLE = ' // XERN1 // &
           ' IS NOT DEFINED.', 12, 1)
        INFO = -12
        return
     ENDIF
!
     if (IND(I) == 3) THEN
        if (BL(I) > BU(I)) THEN
           WRITE (XERN1, '(I8)') I
           WRITE (XERN3, '(1PE15.6)') BL(I)
           WRITE (XERN4, '(1PE15.6)') BU(I)
           call XERMSG ('SLATEC', 'DPLPUP', &
              'IN DSPLP, LOWER BOUND = ' // XERN3 // &
              ' AND UPPER BOUND = ' // XERN4 // &
              ' FOR DEPENDANT VARIABLE = ' // XERN1 // &
              ' ARE NOT CONSISTENT.',13,1)
           INFO = -13
           return
        ENDIF
     ENDIF
   20 CONTINUE
!
!     GET UPDATES OR DATA FOR MATRIX FROM THE USER
!
!     GET THE ELEMENTS OF THE MATRIX FROM THE USER.  IT WILL BE STORED
!     BY COLUMNS USING THE SPARSE STORAGE CODES OF RJ HANSON AND
!     JA WISNIEWSKI.
!
  IFLAG(1) = 1
!
!     KEEP ACCEPTING ELEMENTS UNTIL THE USER IS FINISHED GIVING THEM.
!     LIMIT THIS LOOP TO 2*NVARS*MRELAS ITERATIONS.
!
  ITMAX = 2*NVARS*MRELAS+1
  ITCNT = 0
  FIRST = .TRUE.
!
!     CHECK ON THE ITERATION COUNT.
!
   30 ITCNT = ITCNT+1
  if (ITCNT > ITMAX) THEN
     call XERMSG ('SLATEC', 'DPLPUP', &
        'IN DSPLP, MORE THAN 2*NVARS*MRELAS ITERATIONS DEFINING ' // &
        'OR UPDATING MATRIX DATA.', 7, 1)
     INFO = -7
     return
  end if
!
  AIJ = ZERO
  call DUSRMT(I,J,AIJ,INDCAT,PRGOPT,DATTRV,IFLAG)
  if (IFLAG(1) == 1) THEN
     IFLAG(1) = 2
     go to 30
  end if
!
!     CHECK TO SEE THAT THE SUBSCRIPTS I AND J ARE VALID.
!
  if (I < 1 .OR. I > MRELAS .OR. J < 1 .OR. J > NVARS) THEN
!
!        CHECK ON SIZE OF MATRIX DATA
!        RECORD THE LARGEST AND SMALLEST(IN MAGNITUDE) NONZERO ELEMENTS.
!
     if (IFLAG(1) == 3) THEN
        if (SIZEUP .AND. ABS(AIJ) /= ZERO) THEN
           if (FIRST) THEN
              AMX = ABS(AIJ)
              AMN = ABS(AIJ)
              FIRST = .FALSE.
           ELSEIF (ABS(AIJ) > AMX) THEN
              AMX = ABS(AIJ)
           ELSEIF (ABS(AIJ) < AMN) THEN
              AMN = ABS(AIJ)
           ENDIF
        ENDIF
        go to 40
     ENDIF
!
     WRITE (XERN1, '(I8)') I
     WRITE (XERN2, '(I8)') J
     call XERMSG ('SLATEC', 'DPLPUP', &
        'IN DSPLP, ROW INDEX = ' // XERN1 // ' OR COLUMN INDEX = ' &
        // XERN2 // ' IS OUT OF RANGE.', 8, 1)
     INFO = -8
     return
  end if
!
!     if INDCAT=0 THEN SET A(I,J)=AIJ.
!     if INDCAT=1 THEN ACCUMULATE ELEMENT, A(I,J)=A(I,J)+AIJ.
!
  if (INDCAT == 0) THEN
     call DPCHNG(I,AIJ,IPLACE,AMAT,IMAT,J)
  ELSEIF (INDCAT == 1) THEN
     INDEX = -(I-1)
     call DPNNZR(INDEX,XVAL,IPLACE,AMAT,IMAT,J)
     if (INDEX == I) AIJ=AIJ+XVAL
     call DPCHNG(I,AIJ,IPLACE,AMAT,IMAT,J)
  ELSE
     WRITE (XERN1, '(I8)') INDCAT
     call XERMSG ('SLATEC', 'DPLPUP', &
        'IN DSPLP, INDICATION FLAG = ' // XERN1 // &
        ' FOR MATRIX DATA MUST BE EITHER 0 OR 1.', 9, 1)
     INFO = -9
     return
  end if
!
!     CHECK ON SIZE OF MATRIX DATA
!     RECORD THE LARGEST AND SMALLEST(IN MAGNITUDE) NONZERO ELEMENTS.
!
  if (SIZEUP .AND. ABS(AIJ) /= ZERO) THEN
     if (FIRST) THEN
        AMX = ABS(AIJ)
        AMN = ABS(AIJ)
        FIRST = .FALSE.
     ELSEIF (ABS(AIJ) > AMX) THEN
        AMX = ABS(AIJ)
     ELSEIF (ABS(AIJ) < AMN) THEN
        AMN = ABS(AIJ)
     ENDIF
  end if
  if (IFLAG(1) /= 3) go to 30
!
   40 if (SIZEUP .AND. .NOT. FIRST) THEN
     if (AMN < ASMALL .OR. AMX > ABIG) THEN
        call XERMSG ('SLATEC', 'DPLPUP', &
           'IN DSPLP, A MATRIX ELEMENT''S SIZE IS OUT OF THE ' // &
           'SPECIFIED RANGE.', 22, 1)
        INFO = -22
        return
     ENDIF
  end if
  return
end
