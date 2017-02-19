subroutine DBOLS (W, MDW, MROWS, NCOLS, BL, BU, IND, IOPT, X, &
     RNORM, MODE, RW, IW)
!
!! DBOLS solves the problem E*X = F (in the least squares sense) ...
!  with bounds on selected X values.
!
!***LIBRARY   SLATEC
!***CATEGORY  K1A2A, G2E, G2H1, G2H2
!***TYPE      DOUBLE PRECISION (SBOLS-S, DBOLS-D)
!***KEYWORDS  BOUNDS, CONSTRAINTS, INEQUALITY, LEAST SQUARES, LINEAR
!***AUTHOR  Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!   **** All INPUT and OUTPUT real variables are DOUBLE PRECISION ****
!
!     The user must have dimension statements of the form:
!
!       DIMENSION W(MDW,NCOLS+1), BL(NCOLS), BU(NCOLS),
!      * X(NCOLS+NX), RW(5*NCOLS)
!       INTEGER IND(NCOLS), IOPT(1+NI), IW(2*NCOLS)
!
!     (Here NX=number of extra locations required for option 4; NX=0
!     for no options; NX=NCOLS if this option is in use. Here NI=number
!     of extra locations required for options 1-6; NI=0 for no
!     options.)
!
!   INPUT
!   -----
!
!    --------------------
!    W(MDW,*),MROWS,NCOLS
!    --------------------
!     The array W(*,*) contains the matrix [E:F] on entry. The matrix
!     [E:F] has MROWS rows and NCOLS+1 columns. This data is placed in
!     the array W(*,*) with E occupying the first NCOLS columns and the
!     right side vector F in column NCOLS+1. The row dimension, MDW, of
!     the array W(*,*) must satisfy the inequality MDW .ge. MROWS.
!     Other values of MDW are errors. The values of MROWS and NCOLS
!     must be positive. Other values are errors. There is an exception
!     to this when using option 1 for accumulation of blocks of
!     equations. In that case MROWS is an OUTPUT variable ONLY, and the
!     matrix data for [E:F] is placed in W(*,*), one block of rows at a
!     time.  MROWS contains the number of rows in the matrix after
!     triangularizing several blocks of equations. This is an OUTPUT
!     parameter ONLY when option 1 is used. See IOPT(*) CONTENTS
!     for details about option 1.
!
!    ------------------
!    BL(*),BU(*),IND(*)
!    ------------------
!     These arrays contain the information about the bounds that the
!     solution values are to satisfy. The value of IND(J) tells the
!     type of bound and BL(J) and BU(J) give the explicit values for
!     the respective upper and lower bounds.
!
!    1.    For IND(J)=1, require X(J) .ge. BL(J).
!          (the value of BU(J) is not used.)
!    2.    For IND(J)=2, require X(J) .le. BU(J).
!          (the value of BL(J) is not used.)
!    3.    For IND(J)=3, require X(J) .ge. BL(J) and
!                                X(J) .le. BU(J).
!    4.    For IND(J)=4, no bounds on X(J) are required.
!          (the values of BL(J) and BU(J) are not used.)
!
!     Values other than 1,2,3 or 4 for IND(J) are errors. In the case
!     IND(J)=3 (upper and lower bounds) the condition BL(J) .gt. BU(J)
!     is an error.
!
!    -------
!    IOPT(*)
!    -------
!     This is the array where the user can specify nonstandard options
!     for DBOLSM( ). Most of the time this feature can be ignored by
!     setting the input value IOPT(1)=99. Occasionally users may have
!     needs that require use of the following subprogram options. For
!     details about how to use the options see below: IOPT(*) CONTENTS.
!
!     Option Number   Brief Statement of Purpose
!     ------ ------   ----- --------- -- -------
!           1         Return to user for accumulation of blocks
!                     of least squares equations.
!           2         Check lengths of all arrays used in the
!                     subprogram.
!           3         Standard scaling of the data matrix, E.
!           4         User provides column scaling for matrix E.
!           5         Provide option array to the low-level
!                     subprogram DBOLSM( ).
!           6         Move the IOPT(*) processing pointer.
!          99         No more options to change.
!
!    ----
!    X(*)
!    ----
!     This array is used to pass data associated with option 4. Ignore
!     this parameter if this option is not used. Otherwise see below:
!     IOPT(*) CONTENTS.
!
!    OUTPUT
!    ------
!
!    ----------
!    X(*),RNORM
!    ----------
!     The array X(*) contains a solution (if MODE .ge.0 or .eq.-22) for
!     the constrained least squares problem. The value RNORM is the
!     minimum residual vector length.
!
!    ----
!    MODE
!    ----
!     The sign of MODE determines whether the subprogram has completed
!     normally, or encountered an error condition or abnormal status. A
!     value of MODE .ge. 0 signifies that the subprogram has completed
!     normally. The value of MODE ( >=  0) is the number of variables
!     in an active status: not at a bound nor at the value ZERO, for
!     the case of free variables. A negative value of MODE will be one
!     of the cases -37,-36,...,-22, or -17,...,-2. Values .lt. -1
!     correspond to an abnormal completion of the subprogram. To
!     understand the abnormal completion codes see below: ERROR
!     MESSAGES for DBOLS( ). AN approximate solution will be returned
!     to the user only when max. iterations is reached, MODE=-22.
!     Values for MODE=-37,...,-22 come from the low-level subprogram
!     DBOLSM(). See the section ERROR MESSAGES for DBOLSM() in the
!     documentation for DBOLSM().
!
!    -----------
!    RW(*),IW(*)
!    -----------
!     These are working arrays with 5*NCOLS and 2*NCOLS entries.
!     (normally the user can ignore the contents of these arrays,
!     but they must be dimensioned properly.)
!
!    IOPT(*) CONTENTS
!    ------- --------
!     The option array allows a user to modify internal variables in
!     the subprogram without recompiling the source code. A central
!     goal of the initial software design was to do a good job for most
!     people. Thus the use of options will be restricted to a select
!     group of users. The processing of the option array proceeds as
!     follows: a pointer, here called LP, is initially set to the value
!     1. This value is updated as each option is processed. At the
!     pointer position the option number is extracted and used for
!     locating other information that allows for options to be changed.
!     The portion of the array IOPT(*) that is used for each option is
!     fixed; the user and the subprogram both know how many locations
!     are needed for each option. A great deal of error checking is
!     done by the subprogram on the contents of the option array.
!     Nevertheless it is still possible to give the subprogram optional
!     input that is meaningless. For example option 4 uses the
!     locations X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1) for passing
!     scaling data. The user must manage the allocation of these
!     locations.
!
!   1
!   -
!     This option allows the user to solve problems with a large number
!     of rows compared to the number of variables. The idea is that the
!     subprogram returns to the user (perhaps many times) and receives
!     new least squares equations from the calling program unit.
!     Eventually the user signals "that's all" and then computes the
!     solution with one final call to subprogram DBOLS( ). The value of
!     MROWS is an OUTPUT variable when this option is used. Its value
!     is always in the range 0 .le. MROWS .le. NCOLS+1. It is equal to
!     the number of rows after the triangularization of the entire set
!     of equations. If LP is the processing pointer for IOPT(*), the
!     usage for the sequential processing of blocks of equations is
!
!        IOPT(LP)=1
!        Move block of equations to W(*,*) starting at
!        the first row of W(*,*).
!        IOPT(LP+3)=# of rows in the block; user defined
!
!     The user now calls DBOLS( ) in a loop. The value of IOPT(LP+1)
!     directs the user's action. The value of IOPT(LP+2) points to
!     where the subsequent rows are to be placed in W(*,*).
!
!      .<LOOP
!      . call DBOLS()
!      . if ( IOPT(LP+1)  ==  1) THEN
!      .    IOPT(LP+3)=# OF ROWS IN THE NEW BLOCK; USER DEFINED
!      .    PLACE NEW BLOCK OF IOPT(LP+3) ROWS IN
!      .    W(*,*) STARTING AT ROW IOPT(LP+2).
!      .
!      .    if (  THIS IS THE LAST BLOCK OF EQUATIONS ) THEN
!      .       IOPT(LP+1)=2
!      .<------CYCLE LOOP
!      .    ELSE if (IOPT(LP+1)  ==  2) THEN
!      <-------EXIT LOOP SOLUTION COMPUTED if MODE  >=  0
!      . ELSE
!      . ERROR CONDITION; SHOULD NOT HAPPEN.
!      .<END LOOP
!
!     Use of this option adds 4 to the required length of IOPT(*).
!
!
!   2
!   -
!     This option is useful for checking the lengths of all arrays used
!     by DBOLS() against their actual requirements for this problem.
!     The idea is simple: the user's program unit passes the declared
!     dimension information of the arrays. These values are compared
!     against the problem-dependent needs within the subprogram. If any
!     of the dimensions are too small an error message is printed and a
!     negative value of MODE is returned, -11 to -17. The printed error
!     message tells how long the dimension should be. If LP is the
!     processing pointer for IOPT(*),
!
!        IOPT(LP)=2
!        IOPT(LP+1)=Row dimension of W(*,*)
!        IOPT(LP+2)=Col. dimension of W(*,*)
!        IOPT(LP+3)=Dimensions of BL(*),BU(*),IND(*)
!        IOPT(LP+4)=Dimension of X(*)
!        IOPT(LP+5)=Dimension of RW(*)
!        IOPT(LP+6)=Dimension of IW(*)
!        IOPT(LP+7)=Dimension of IOPT(*)
!         .
!        call DBOLS()
!
!     Use of this option adds 8 to the required length of IOPT(*).
!
!   3
!   -
!     This option changes the type of scaling for the data matrix E.
!     Nominally each nonzero column of E is scaled so that the
!     magnitude of its largest entry is equal to the value ONE. If LP
!     is the processing pointer for IOPT(*),
!
!        IOPT(LP)=3
!        IOPT(LP+1)=1,2 or 3
!            1= Nominal scaling as noted;
!            2= Each nonzero column scaled to have length ONE;
!            3= Identity scaling; scaling effectively suppressed.
!         .
!        call DBOLS()
!
!     Use of this option adds 2 to the required length of IOPT(*).
!
!   4
!   -
!     This option allows the user to provide arbitrary (positive)
!     column scaling for the matrix E. If LP is the processing pointer
!     for IOPT(*),
!
!        IOPT(LP)=4
!        IOPT(LP+1)=IOFF
!        X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1)
!        = Positive scale factors for cols. of E.
!         .
!        call DBOLS()
!
!     Use of this option adds 2 to the required length of IOPT(*) and
!     NCOLS to the required length of X(*).
!
!   5
!   -
!     This option allows the user to provide an option array to the
!     low-level subprogram DBOLSM(). If LP is the processing pointer
!     for IOPT(*),
!
!        IOPT(LP)=5
!        IOPT(LP+1)= Position in IOPT(*) where option array
!                    data for DBOLSM() begins.
!         .
!        call DBOLS()
!
!     Use of this option adds 2 to the required length of IOPT(*).
!
!   6
!   -
!     Move the processing pointer (either forward or backward) to the
!     location IOPT(LP+1). The processing point is moved to entry
!     LP+2 of IOPT(*) if the option is left with -6 in IOPT(LP).  For
!     example to skip over locations 3,...,NCOLS+2 of IOPT(*),
!
!       IOPT(1)=6
!       IOPT(2)=NCOLS+3
!       (IOPT(I), I=3,...,NCOLS+2 are not defined here.)
!       IOPT(NCOLS+3)=99
!       call DBOLS()
!
!     CAUTION: Misuse of this option can yield some very hard
!     -to-find bugs.  Use it with care.
!
!   99
!   --
!     There are no more options to change.
!
!     Only option numbers -99, -6,-5,...,-1, 1,2,...,6, and 99 are
!     permitted. Other values are errors. Options -99,-1,...,-6 mean
!     that the respective options 99,1,...,6 are left at their default
!     values. An example is the option to modify the (rank) tolerance:
!
!       IOPT(1)=-3 Option is recognized but not changed
!       IOPT(2)=2  Scale nonzero cols. to have length ONE
!       IOPT(3)=99
!
!    ERROR MESSAGES for DBOLS()
!    ----- -------- --- -------
!
! WARNING IN...
! DBOLS(). MDW=(I1) MUST BE POSITIVE.
!           IN ABOVE MESSAGE, I1=         0
! ERROR NUMBER =         2
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
! WARNING IN...
! DBOLS(). NCOLS=(I1) THE NO. OF VARIABLES MUST BE POSITIVE.
!           IN ABOVE MESSAGE, I1=         0
! ERROR NUMBER =         3
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
! WARNING IN...
! DBOLS(). FOR J=(I1), IND(J)=(I2) MUST BE 1-4.
!           IN ABOVE MESSAGE, I1=         1
!           IN ABOVE MESSAGE, I2=         0
! ERROR NUMBER =         4
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
! WARNING IN...
! DBOLS(). FOR J=(I1), BOUND BL(J)=(R1) IS  >  BU(J)=(R2).
!           IN ABOVE MESSAGE, I1=         1
!           IN ABOVE MESSAGE, R1=    0.
!           IN ABOVE MESSAGE, R2=    ABOVE MESSAGE, I1=         0
! ERROR NUMBER =         6
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
! WARNING IN...
! DBOLS(). ISCALE OPTION=(I1) MUST BE 1-3.
!           IN ABOVE MESSAGE, I1=         0
! ERROR NUMBER =         7
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
! WARNING IN...
! DBOLS(). OFFSET PAST X(NCOLS) (I1) FOR USER-PROVIDED  COLUMN SCALING
! MUST BE POSITIVE.
!           IN ABOVE MESSAGE, I1=         0
! ERROR NUMBER =         8
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
! WARNING IN...
! DBOLS(). EACH PROVIDED COL. SCALE FACTOR MUST BE POSITIVE.
! COMPONENT (I1) NOW = (R1).
!           IN ABOVE MESSAGE, I1=        ND.  <=  MDW=(I2).
!           IN ABOVE MESSAGE, I1=         1
!           IN ABOVE MESSAGE, I2=         0
! ERROR NUMBER =        10
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
! WARNING IN...
! DBOLS().THE ROW DIMENSION OF W(,)=(I1) MUST BE  >= THE NUMBER OF ROWS=
! (I2).
!           IN ABOVE MESSAGE, I1=         0
!           IN ABOVE MESSAGE, I2=         1
! ERROR NUMBER =        11
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
! WARNING IN...
! DBOLS(). THE COLUMN DIMENSION OF W(,)=(I1) MUST BE  >=  NCOLS+1=(I2).
!           IN ABOVE MESSAGE, I1=         0
!           IN ABOVE MESSAGE, I2=         2
! ERROR NUMBER =        12
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
! WARNING IN...
! DBOLS().THE DIMENSIONS OF THE ARRAYS BL(),BU(), AND IND()=(I1) MUST BE
!  >=  NCOLS=(I2).
!           IN ABOVE MESSAGE, I1=         0
!           IN ABOVE MESSAGE, I2=         1
! ERROR NUMBER =        13
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
! WARNING IN...
! DBOLS(). THE DIMENSION OF X()=(I1) MUST BE  >=  THE REQD. LENGTH=(I2).
!           IN ABOVE MESSAGE, I1=         0
!           IN ABOVE MESSAGE, I2=         2
! ERROR NUMBER =        14
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
! WARNING IN...
! DBOLS(). THE DIMENSION OF RW()=(I1) MUST BE  >=  5*NCOLS=(I2).
!           IN ABOVE MESSAGE, I1=         0
!           IN ABOVE MESSAGE, I2=         3
! ERROR NUMBER =        15
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
! WARNING IN...
! DBOLS() THE DIMENSION OF IW()=(I1) MUST BE  >=  2*NCOLS=(I2).
!           IN ABOVE MESSAGE, I1=         0
!           IN ABOVE MESSAGE, I2=         2
! ERROR NUMBER =        16
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
! WARNING IN...
! DBOLS() THE DIMENSION OF IOPT()=(I1) MUST BE  >=  THE REQD. LEN.=(I2).
!           IN ABOVE MESSAGE, I1=         0
!           IN ABOVE MESSAGE, I2=         1
! ERROR NUMBER =        17
! (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
!
!***REFERENCES  R. J. Hanson, Linear least squares with bounds and
!                 linear constraints, Report SAND82-1517, Sandia
!                 Laboratories, August 1982.
!***ROUTINES CALLED  DBOLSM, DCOPY, DNRM2, DROT, DROTG, IDAMAX, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   821220  DATE WRITTEN
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBOLS
!
!     SOLVE LINEAR LEAST SQUARES SYSTEM WITH BOUNDS ON
!     SELECTED VARIABLES.
!     REVISED 850329-1400
!     REVISED YYMMDD-HHMM
!     TO CHANGE THIS SUBPROGRAM FROM SINGLE TO DOUBLE PRECISION BEGIN
!     EDITING AT THE CARD 'C++'.
!     CHANGE THIS SUBPROGRAM NAME TO DBOLS AND THE STRINGS
!     /SCOPY/ TO /DCOPY/, /SBOL/ TO /DBOL/,
!     /SNRM2/ TO /DNRM2/, /ISAMAX/ TO /IDAMAX/,
!     /SROTG/ TO /DROTG/, /SROT/ TO /DROT/, /E0/ TO /D0/,
!     /REAL            / TO /DOUBLE PRECISION/.
! ++
  DOUBLE PRECISION W(MDW,*),BL(*),BU(*),X(*),RW(*)
  DOUBLE PRECISION SC, SS, ONE, DNRM2, RNORM, ZERO
!
!     THIS VARIABLE SHOULD REMAIN TYPE REAL.
  INTEGER IND(*),IOPT(*),IW(*)
  LOGICAL CHECKL
  CHARACTER*8 XERN1, XERN2
  CHARACTER*16 XERN3, XERN4
  SAVE IGO,LOCACC,LOPT,ISCALE
  DATA IGO/0/
!***FIRST EXECUTABLE STATEMENT  DBOLS
  NERR = 0
  MODE = 0
  if (IGO == 0) THEN
!     DO(CHECK VALIDITY OF INPUT DATA)
!     PROCEDURE(CHECK VALIDITY OF INPUT DATA)
!
!     SEE THAT MDW IS  > 0. GROSS CHECK ONLY.
      if (MDW <= 0) THEN
          WRITE (XERN1, '(I8)') MDW
          call XERMSG ('SLATEC', 'DBOLS', 'MDW = ' // XERN1 // &
             ' MUST BE POSITIVE.', 2, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
          go to 190
      ENDIF
!
!     SEE THAT NUMBER OF UNKNOWNS IS POSITIVE.
      if (NCOLS <= 0) THEN
          WRITE (XERN1, '(I8)') NCOLS
          call XERMSG ('SLATEC', 'DBOLS', 'NCOLS = ' // XERN1 // &
             ' THE NO. OF VARIABLES MUST BE POSITIVE.', 3, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
          go to 190
      ENDIF
!
!     SEE THAT CONSTRAINT INDICATORS ARE ALL WELL-DEFINED.
      DO 10 J = 1,NCOLS
          if (IND(J) < 1 .OR. IND(J) > 4) THEN
              WRITE (XERN1, '(I8)') J
              WRITE (XERN2, '(I8)') IND(J)
              call XERMSG ('SLATEC', 'DBOLS', 'IND(' // XERN1 // &
                 ') = ' // XERN2 // ' MUST BE 1-4.', 4, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 190
          ENDIF
   10     CONTINUE
!
!     SEE THAT BOUNDS ARE CONSISTENT.
      DO 20 J = 1,NCOLS
          if (IND(J) == 3) THEN
              if (BL(J) > BU(J)) THEN
                  WRITE (XERN1, '(I8)') J
                  WRITE (XERN3, '(1PE15.6)') BL(J)
                  WRITE (XERN4, '(1PE15.6)') BU(J)
                  call XERMSG ('SLATEC', 'DBOLS', 'BOUND BL(' // &
                     XERN1 // ') = ' // XERN3 // ' IS  >  BU(' // &
                     XERN1 // ') = ' // XERN4, 5, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
                  go to 190
              ENDIF
          ENDIF
   20     CONTINUE
!     END PROCEDURE
!     DO(PROCESS OPTION ARRAY)
!     PROCEDURE(PROCESS OPTION ARRAY)
      ZERO = 0.D0
      ONE = 1.D0
      CHECKL = .FALSE.
      LENX = NCOLS
      ISCALE = 1
      IGO = 2
      LOPT = 0
      LP = 0
      LDS = 0
   30     CONTINUE
      LP = LP + LDS
      IP = IOPT(LP+1)
      JP = ABS(IP)
!
!     TEST FOR NO MORE OPTIONS.
      if (IP == 99) THEN
          if (LOPT == 0) LOPT = LP + 1
          go to 50
      ELSE if (JP == 99) THEN
          LDS = 1
          go to 30
      ELSE if (JP == 1) THEN
          if (IP > 0) THEN
!
!     SET UP DIRECTION FLAG, ROW STACKING POINTER
!     LOCATION, AND LOCATION FOR NUMBER OF NEW ROWS.
              LOCACC = LP + 2
!
!                  IOPT(LOCACC-1)=OPTION NUMBER FOR SEQ. ACCUMULATION.
!     CONTENTS..   IOPT(LOCACC  )=USER DIRECTION FLAG, 1 OR 2.
!                  IOPT(LOCACC+1)=ROW STACKING POINTER.
!                  IOPT(LOCACC+2)=NUMBER OF NEW ROWS TO PROCESS.
!     USER ACTION WITH THIS OPTION..
!      (SET UP OPTION DATA FOR SEQ. ACCUMULATION IN IOPT(*).
!      MUST ALSO START PROCESS WITH IOPT(LOCACC)=1.)
!      (MOVE BLOCK OF EQUATIONS INTO W(*,*)  STARTING AT FIRST
!       ROW OF W(*,*).  SET IOPT(LOCACC+2)=NO. OF ROWS IN BLOCK.)
!              LOOP
!              call DBOLS()
!
!                  if ( IOPT(LOCACC)  ==  1) THEN
!                      STACK EQUAS., STARTING AT ROW IOPT(LOCACC+1),
!                       INTO W(*,*).
!                       SET IOPT(LOCACC+2)=NO. OF EQUAS.
!                      if LAST BLOCK OF EQUAS., SET IOPT(LOCACC)=2.
!                  ELSE if IOPT(LOCACC)  ==  2) THEN
!                      (PROCESS IS OVER. EXIT LOOP.)
!                  ELSE
!                      (ERROR CONDITION. SHOULD NOT HAPPEN.)
!                  end if
!              END LOOP
!              SET IOPT(LOCACC-1)=-OPTION NUMBER FOR SEQ. ACCUMULATION.
!              call DBOLS( )
              IOPT(LOCACC+1) = 1
              IGO = 1
          ENDIF
          LDS = 4
          go to 30
      ELSE if (JP == 2) THEN
          if (IP > 0) THEN
!
!     GET ACTUAL LENGTHS OF ARRAYS FOR CHECKING AGAINST NEEDS.
              LOCDIM = LP + 2
!
!     LMDW >= MROWS
!     LNDW >= NCOLS+1
!     LLB  >= NCOLS
!     LLX  >= NCOLS+EXTRA REQD. IN OPTIONS.
!     LLRW >= 5*NCOLS
!     LLIW >= 2*NCOLS
!     LIOP >=  AMOUNT REQD. FOR IOPTION ARRAY.
              LMDW = IOPT(LOCDIM)
              LNDW = IOPT(LOCDIM+1)
              LLB = IOPT(LOCDIM+2)
              LLX = IOPT(LOCDIM+3)
              LLRW = IOPT(LOCDIM+4)
              LLIW = IOPT(LOCDIM+5)
              LIOPT = IOPT(LOCDIM+6)
              CHECKL = .TRUE.
          ENDIF
          LDS = 8
          go to 30
!
!     OPTION TO MODIFY THE COLUMN SCALING.
      ELSE if (JP == 3) THEN
          if (IP > 0) THEN
              ISCALE = IOPT(LP+2)
!
!     SEE THAT ISCALE IS 1 THRU 3.
              if (ISCALE < 1 .OR. ISCALE > 3) THEN
                  WRITE (XERN1, '(I8)') ISCALE
                  call XERMSG ('SLATEC', 'DBOLS', 'ISCALE OPTION = ' &
                     // XERN1 // ' MUST BE 1-3', 7, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
                  go to 190
              ENDIF
          ENDIF
          LDS = 2
!     CYCLE FOREVER
          go to 30
!
!     IN THIS OPTION THE USER HAS PROVIDED SCALING.  THE
!     SCALE FACTORS FOR THE COLUMNS BEGIN IN X(NCOLS+IOPT(LP+2)).
      ELSE if (JP == 4) THEN
          if (IP > 0) THEN
              ISCALE = 4
              if (IOPT(LP+2) <= 0) THEN
                  WRITE (XERN1, '(I8)') IOPT(LP+2)
                  call XERMSG ('SLATEC', 'DBOLS', &
                     'OFFSET PAST X(NCOLS) (' // XERN1 // &
                     ') FOR USER-PROVIDED COLUMN SCALING MUST ' // &
                     'BE POSITIVE.',  8, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
                  go to 190
              ENDIF
              call DCOPY(NCOLS,X(NCOLS+IOPT(LP+2)),1,RW,1)
              LENX = LENX + NCOLS
              DO 40 J = 1,NCOLS
                  if (RW(J) <= ZERO) THEN
                      WRITE (XERN1, '(I8)') J
                      WRITE (XERN3, '(1PE15.6)') RW(J)
                      call XERMSG ('SLATEC', 'DBOLS', &
                         'EACH PROVIDED COLUMN SCALE FACTOR ' // &
                         'MUST BE POSITIVE.$$COMPONENT ' // XERN1 // &
                         ' NOW = ' // XERN3, 9, 1)
                      go to 190
                  ENDIF
   40             CONTINUE
          ENDIF
          LDS = 2
!     CYCLE FOREVER
          go to 30
!
!     IN THIS OPTION AN OPTION ARRAY IS PROVIDED TO DBOLSM().
      ELSE if (JP == 5) THEN
          if (IP > 0) THEN
              LOPT = IOPT(LP+2)
          ENDIF
          LDS = 2
!     CYCLE FOREVER
          go to 30
!
!     THIS OPTION USES THE NEXT LOC OF IOPT(*) AS AN
!     INCREMENT TO SKIP.
      ELSE if (JP == 6) THEN
          if (IP > 0) THEN
              LP = IOPT(LP+2) - 1
              LDS = 0
          ELSE
              LDS = 2
          ENDIF
!     CYCLE FOREVER
          go to 30
!
!     NO VALID OPTION NUMBER WAS NOTED. THIS IS AN ERROR CONDITION.
      ELSE
          WRITE (XERN1, '(I8)') JP
          call XERMSG ('SLATEC', 'DBOLS', 'THE OPTION NUMBER = ' // &
             XERN1 // ' IS NOT DEFINED.', 6, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
          go to 190
      ENDIF
   50     CONTINUE
!     END PROCEDURE
      if (CHECKL) THEN
!     DO(CHECK LENGTHS OF ARRAYS)
!     PROCEDURE(CHECK LENGTHS OF ARRAYS)
!
!     THIS FEATURE ALLOWS THE USER TO MAKE SURE THAT THE
!     ARRAYS ARE LONG ENOUGH FOR THE INTENDED PROBLEM SIZE AND USE.
          if (LMDW < MROWS) THEN
              WRITE (XERN1, '(I8)') LMDW
              WRITE (XERN2, '(I8)') MROWS
              call XERMSG ('SLATEC', 'DBOLS', &
                 'THE ROW DIMENSION OF W(,) = ' // XERN1 // &
                 ' MUST BE  >=  THE NUMBER OF ROWS = ' // XERN2, &
                 11, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 190
          ENDIF
          if (LNDW < NCOLS+1) THEN
              WRITE (XERN1, '(I8)') LNDW
              WRITE (XERN2, '(I8)') NCOLS+1
              call XERMSG ('SLATEC', 'DBOLS', &
                 'THE COLUMN DIMENSION OF W(,) = ' // XERN1 // &
                 ' MUST BE  >=  NCOLS+1 = ' // XERN2, 12, 1)
              go to 190
          ENDIF
          if (LLB < NCOLS) THEN
              WRITE (XERN1, '(I8)') LLB
              WRITE (XERN2, '(I8)') NCOLS
              call XERMSG ('SLATEC', 'DBOLS', &
             'THE DIMENSIONS OF THE ARRAYS BL(), BU(), AND IND() = ' &
                 // XERN1 // ' MUST BE  >=  NCOLS = ' // XERN2, &
                 13, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 190
          ENDIF
          if (LLX < LENX) THEN
              WRITE (XERN1, '(I8)') LLX
              WRITE (XERN2, '(I8)') LENX
              call XERMSG ('SLATEC', 'DBOLS', &
                 'THE DIMENSION OF X() = ' // XERN1 // &
                 ' MUST BE  >=  THE REQUIRED LENGTH = ' // XERN2, &
                 14, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 190
          ENDIF
          if (LLRW < 5*NCOLS) THEN
              WRITE (XERN1, '(I8)') LLRW
              WRITE (XERN2, '(I8)') 5*NCOLS
              call XERMSG ('SLATEC', 'DBOLS', &
                 'THE DIMENSION OF RW() = ' // XERN1 // &
                 ' MUST BE  >=  5*NCOLS = ' // XERN2, 15, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 190
          ENDIF
          if (LLIW < 2*NCOLS) THEN
              WRITE (XERN1, '(I8)') LLIW
              WRITE (XERN2, '(I8)') 2*NCOLS
              call XERMSG ('SLATEC', 'DBOLS', &
                 'THE DIMENSION OF IW() = ' // XERN1 // &
                 ' MUST BE  >=  2*NCOLS = ' // XERN2, 16, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 190
          ENDIF
          if (LIOPT < LP+1) THEN
              WRITE (XERN1, '(I8)') LIOPT
              WRITE (XERN2, '(I8)') LP+1
              call XERMSG ('SLATEC', 'DBOLS', &
                 'THE DIMENSION OF IOPT() = ' // XERN1 // &
                 ' MUST BE  >=  THE REQUIRED LEN = ' // XERN2, 17,1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 190
          ENDIF
!     END PROCEDURE
      ENDIF
  end if
  go to (60,90),IGO
  go to 180
!
!     GO BACK TO THE USER FOR ACCUMULATION OF LEAST SQUARES
!     EQUATIONS AND DIRECTIONS TO QUIT PROCESSING.
!     CASE 1
   60 CONTINUE
!     DO(ACCUMULATE LEAST SQUARES EQUATIONS)
!     PROCEDURE(ACCUMULATE LEAST SQUARES EQUATIONS)
  MROWS = IOPT(LOCACC+1) - 1
  INROWS = IOPT(LOCACC+2)
  MNEW = MROWS + INROWS
  if (MNEW < 0 .OR. MNEW > MDW) THEN
      WRITE (XERN1, '(I8)') MNEW
      WRITE (XERN2, '(I8)') MDW
      call XERMSG ('SLATEC', 'DBOLS', 'NO. OF ROWS = ' // XERN1 // &
         ' MUST BE  >=  0 .AND.  <=  MDW = ' // XERN2, 10, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
      go to 190
  end if
  DO 80 J = 1,MIN(NCOLS+1,MNEW)
      DO 70 I = MNEW,MAX(MROWS,J) + 1,-1
          IBIG = IDAMAX(I-J,W(J,J),1) + J - 1
!
!     PIVOT FOR INCREASED STABILITY.
          call DROTG(W(IBIG,J),W(I,J),SC,SS)
          call DROT(NCOLS+1-J,W(IBIG,J+1),MDW,W(I,J+1),MDW,SC,SS)
          W(I,J) = ZERO
   70     CONTINUE
   80 CONTINUE
  MROWS = MIN(NCOLS+1,MNEW)
  IOPT(LOCACC+1) = MROWS + 1
  IGO = IOPT(LOCACC)
!     END PROCEDURE
  if (IGO == 2) THEN
      IGO = 0
  end if
  go to 180
!     CASE 2
   90 CONTINUE
!     DO(INITIALIZE VARIABLES AND DATA VALUES)
!     PROCEDURE(INITIALIZE VARIABLES AND DATA VALUES)
  DO 150 J = 1,NCOLS
      go to (100,110,120,130),ISCALE
      go to 140
  100     CONTINUE
!     CASE 1
!
!     THIS IS THE NOMINAL SCALING. EACH NONZERO
!     COL. HAS MAX. NORM EQUAL TO ONE.
      IBIG = IDAMAX(MROWS,W(1,J),1)
      RW(J) = ABS(W(IBIG,J))
      if (RW(J) == ZERO) THEN
          RW(J) = ONE
      ELSE
          RW(J) = ONE/RW(J)
      ENDIF
      go to 140
  110     CONTINUE
!     CASE 2
!
!     THIS CHOICE OF SCALING MAKES EACH NONZERO COLUMN
!     HAVE EUCLIDEAN LENGTH EQUAL TO ONE.
      RW(J) = DNRM2(MROWS,W(1,J),1)
      if (RW(J) == ZERO) THEN
          RW(J) = ONE
      ELSE
          RW(J) = ONE/RW(J)
      ENDIF
      go to 140
  120     CONTINUE
!     CASE 3
!
!     THIS CASE EFFECTIVELY SUPPRESSES SCALING BY SETTING
!     THE SCALING MATRIX TO THE IDENTITY MATRIX.
      RW(1) = ONE
      call DCOPY(NCOLS,RW,0,RW,1)
      go to 160
  130     CONTINUE
!     CASE 4
      go to 160
  140     CONTINUE
  150 CONTINUE
  160 CONTINUE
!     END PROCEDURE
!     DO(SOLVE BOUNDED LEAST SQUARES PROBLEM)
!     PROCEDURE(SOLVE BOUNDED LEAST SQUARES PROBLEM)
!
!     INITIALIZE IBASIS(*), J=1,NCOLS, AND IBB(*), J=1,NCOLS,
!     TO =J,AND =1, FOR USE IN DBOLSM( ).
  DO 170 J = 1,NCOLS
      IW(J) = J
      IW(J+NCOLS) = 1
      RW(3*NCOLS+J) = BL(J)
      RW(4*NCOLS+J) = BU(J)
  170 CONTINUE
  call DBOLSM(W,MDW,MROWS,NCOLS,RW(3*NCOLS+1),RW(4*NCOLS+1),IND, &
              IOPT(LOPT),X,RNORM,MODE,RW(NCOLS+1),RW(2*NCOLS+1),RW, &
              IW,IW(NCOLS+1))
!     END PROCEDURE
  IGO = 0
  180 CONTINUE
  return
!     PROCEDURE(RETURN TO USER PROGRAM UNIT)
  190 if ( MODE >= 0)MODE = -NERR
  IGO = 0
  return
!     END PROCEDURE
end
