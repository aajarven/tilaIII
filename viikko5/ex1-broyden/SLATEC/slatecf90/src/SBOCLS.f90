subroutine SBOCLS (W, MDW, MCON, MROWS, NCOLS, BL, BU, IND, IOPT, &
     X, RNORMC, RNORM, MODE, RW, IW)
!
!! SBOCLS solves the bounded and constrained least squares problem ...
!  consisting of solving the equation
!                      E*X = F  (in the least squares sense)
!             subject to the linear constraints
!                            C*X = Y.
!***LIBRARY   SLATEC
!***CATEGORY  K1A2A, G2E, G2H1, G2H2
!***TYPE      SINGLE PRECISION (SBOCLS-S, DBOCLS-D)
!***KEYWORDS  BOUNDS, CONSTRAINTS, INEQUALITY, LEAST SQUARES, LINEAR
!***AUTHOR  Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!     This subprogram solves the bounded and constrained least squares
!     problem. The problem statement is:
!
!     Solve E*X = F (least squares sense), subject to constraints
!     C*X=Y.
!
!     In this formulation both X and Y are unknowns, and both may
!     have bounds on any of their components.  This formulation
!     of the problem allows the user to have equality and inequality
!     constraints as well as simple bounds on the solution components.
!
!     This constrained linear least squares subprogram solves E*X=F
!     subject to C*X=Y, where E is MROWS by NCOLS, C is MCON by NCOLS.
!
!      The user must have dimension statements of the form
!
!      DIMENSION W(MDW,NCOLS+MCON+1), BL(NCOLS+MCON), BU(NCOLS+MCON),
!     * X(2*(NCOLS+MCON)+2+NX), RW(6*NCOLS+5*MCON)
!       INTEGER IND(NCOLS+MCON), IOPT(17+NI), IW(2*(NCOLS+MCON))
!
!     (here NX=number of extra locations required for the options; NX=0
!     if no options are in use. Also NI=number of extra locations
!     for options 1-9.)
!
!    INPUT
!    -----
!
!    -------------------------
!    W(MDW,*),MCON,MROWS,NCOLS
!    -------------------------
!     The array W contains the (possibly null) matrix [C:*] followed by
!     [E:F].  This must be placed in W as follows:
!          [C  :  *]
!     W  = [       ]
!          [E  :  F]
!     The (*) after C indicates that this data can be undefined. The
!     matrix [E:F] has MROWS rows and NCOLS+1 columns. The matrix C is
!     placed in the first MCON rows of W(*,*) while [E:F]
!     follows in rows MCON+1 through MCON+MROWS of W(*,*). The vector F
!     is placed in rows MCON+1 through MCON+MROWS, column NCOLS+1. The
!     values of MDW and NCOLS must be positive; the value of MCON must
!     be nonnegative. An exception to this occurs when using option 1
!     for accumulation of blocks of equations. In that case MROWS is an
!     OUTPUT variable only, and the matrix data for [E:F] is placed in
!     W(*,*), one block of rows at a time. See IOPT(*) contents, option
!     number 1, for further details. The row dimension, MDW, of the
!     array W(*,*) must satisfy the inequality:
!
!     If using option 1,
!                     MDW .ge. MCON + max(max. number of
!                     rows accumulated, NCOLS) + 1.
!     If using option 8,
!                     MDW .ge. MCON + MROWS.
!     Else
!                     MDW .ge. MCON + max(MROWS, NCOLS).
!
!     Other values are errors, but this is checked only when using
!     option=2.  The value of MROWS is an output parameter when
!     using option number 1 for accumulating large blocks of least
!     squares equations before solving the problem.
!     See IOPT(*) contents for details about option 1.
!
!    ------------------
!    BL(*),BU(*),IND(*)
!    ------------------
!     These arrays contain the information about the bounds that the
!     solution values are to satisfy. The value of IND(J) tells the
!     type of bound and BL(J) and BU(J) give the explicit values for
!     the respective upper and lower bounds on the unknowns X and Y.
!     The first NVARS entries of IND(*), BL(*) and BU(*) specify
!     bounds on X; the next MCON entries specify bounds on Y.
!
!    1.    For IND(J)=1, require X(J) .ge. BL(J);
!          if J.gt.NCOLS,        Y(J-NCOLS) .ge. BL(J).
!          (the value of BU(J) is not used.)
!    2.    For IND(J)=2, require X(J) .le. BU(J);
!          if J.gt.NCOLS,        Y(J-NCOLS) .le. BU(J).
!          (the value of BL(J) is not used.)
!    3.    For IND(J)=3, require X(J) .ge. BL(J) and
!                                X(J) .le. BU(J);
!          if J.gt.NCOLS,        Y(J-NCOLS) .ge. BL(J) and
!                                Y(J-NCOLS) .le. BU(J).
!          (to impose equality constraints have BL(J)=BU(J)=
!          constraining value.)
!    4.    For IND(J)=4, no bounds on X(J) or Y(J-NCOLS) are required.
!          (the values of BL(J) and BU(J) are not used.)
!
!     Values other than 1,2,3 or 4 for IND(J) are errors. In the case
!     IND(J)=3 (upper and lower bounds) the condition BL(J) .gt. BU(J)
!     is  an  error.   The values BL(J), BU(J), J .gt. NCOLS, will be
!     changed.  Significant changes mean that the constraints are
!     infeasible.  (Users must make this decision themselves.)
!     The new values for BL(J), BU(J), J .gt. NCOLS, define a
!     region such that the perturbed problem is feasible.  If users
!     know that their problem is feasible, this step can be skipped
!     by using option number 8 described below.
!
!     See IOPT(*) description.
!
!
!    -------
!    IOPT(*)
!    -------
!     This is the array where the user can specify nonstandard options
!     for SBOCLS( ). Most of the time this feature can be ignored by
!     setting the input value IOPT(1)=99. Occasionally users may have
!     needs that require use of the following subprogram options. For
!     details about how to use the options see below: IOPT(*) CONTENTS.
!
!     Option Number   Brief Statement of Purpose
!     ------ ------   ----- --------- -- -------
!           1         Return to user for accumulation of blocks
!                     of least squares equations.  The values
!                     of IOPT(*) are changed with this option.
!                     The changes are updates to pointers for
!                     placing the rows of equations into position
!                     for processing.
!           2         Check lengths of all arrays used in the
!                     subprogram.
!           3         Column scaling of the data matrix, [C].
!                                                        [E]
!           4         User provides column scaling for matrix [C].
!                                                             [E]
!           5         Provide option array to the low-level
!                     subprogram SBOLS( ).
!           6         Provide option array to the low-level
!                     subprogram SBOLSM( ).
!           7         Move the IOPT(*) processing pointer.
!           8         Do not preprocess the constraints to
!                     resolve infeasibilities.
!           9         Do not pretriangularize the least squares matrix.
!          99         No more options to change.
!
!    ----
!    X(*)
!    ----
!     This array is used to pass data associated with options 4,5 and
!     6. Ignore this parameter (on input) if no options are used.
!     Otherwise see below: IOPT(*) CONTENTS.
!
!
!    OUTPUT
!    ------
!
!    -----------------
!    X(*),RNORMC,RNORM
!    -----------------
!     The array X(*) contains a solution (if MODE .ge.0 or .eq.-22) for
!     the constrained least squares problem. The value RNORMC is the
!     minimum residual vector length for the constraints C*X - Y = 0.
!     The value RNORM is the minimum residual vector length for the
!     least squares equations. Normally RNORMC=0, but in the case of
!     inconsistent constraints this value will be nonzero.
!     The values of X are returned in the first NVARS entries of X(*).
!     The values of Y are returned in the last MCON entries of X(*).
!
!    ----
!    MODE
!    ----
!     The sign of MODE determines whether the subprogram has completed
!     normally, or encountered an error condition or abnormal status. A
!     value of MODE .ge. 0 signifies that the subprogram has completed
!     normally. The value of mode (.ge. 0) is the number of variables
!     in an active status: not at a bound nor at the value zero, for
!     the case of free variables. A negative value of MODE will be one
!     of the cases (-57)-(-41), (-37)-(-22), (-19)-(-2). Values .lt. -1
!     correspond to an abnormal completion of the subprogram. These
!     error messages are in groups for the subprograms SBOCLS(),
!     SBOLSM(), and SBOLS().  An approximate solution will be returned
!     to the user only when max. iterations is reached, MODE=-22.
!
!    -----------
!    RW(*),IW(*)
!    -----------
!     These are working arrays.  (normally the user can ignore the
!     contents of these arrays.)
!
!    IOPT(*) CONTENTS
!    ------- --------
!     The option array allows a user to modify some internal variables
!     in the subprogram without recompiling the source code. A central
!     goal of the initial software design was to do a good job for most
!     people. Thus the use of options will be restricted to a select
!     group of users. The processing of the option array proceeds as
!     follows: a pointer, here called LP, is initially set to the value
!     1. At the pointer position the option number is extracted and
!     used for locating other information that allows for options to be
!     changed. The portion of the array IOPT(*) that is used for each
!     option is fixed; the user and the subprogram both know how many
!     locations are needed for each option. The value of LP is updated
!     for each option based on the amount of storage in IOPT(*) that is
!     required. A great deal of error checking is done by the
!     subprogram on the contents of the option array. Nevertheless it
!     is still possible to give the subprogram optional input that is
!     meaningless. For example option 4 uses the locations
!     X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1) for passing scaling data.
!     The user must manage the allocation of these locations.
!
!   1
!   -
!     This option allows the user to solve problems with a large number
!     of rows compared to the number of variables. The idea is that the
!     subprogram returns to the user (perhaps many times) and receives
!     new least squares equations from the calling program unit.
!     Eventually the user signals "that's all" and a solution is then
!     computed. The value of MROWS is an output variable when this
!     option is used. Its value is always in the range 0 .le. MROWS
!     .le. NCOLS+1. It is the number of rows after the
!     triangularization of the entire set of equations. If LP is the
!     processing pointer for IOPT(*), the usage for the sequential
!     processing of blocks of equations is
!
!
!        IOPT(LP)=1
!         Move block of equations to W(*,*) starting at
!         the first row of W(*,*).
!        IOPT(LP+3)=# of rows in the block; user defined
!
!     The user now calls SBOCLS( ) in a loop. The value of IOPT(LP+1)
!     directs the user's action. The value of IOPT(LP+2) points to
!     where the subsequent rows are to be placed in W(*,*). Both of
!     these values are first defined in the subprogram. The user
!     changes the value of IOPT(LP+1) (to 2) as a signal that all of
!     the rows have been processed.
!
!
!      .<LOOP
!      . call SBOCLS( )
!      . if ( IOPT(LP+1)  ==  1) THEN
!      .    IOPT(LP+3)=# OF ROWS IN THE NEW BLOCK; USER DEFINED
!      .    PLACE NEW BLOCK OF IOPT(LP+3) ROWS IN
!      .    W(*,*) STARTING AT ROW MCON + IOPT(LP+2).
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
!   2
!   -
!     This option is useful for checking the lengths of all arrays used
!     by SBOCLS( ) against their actual requirements for this problem.
!     The idea is simple: the user's program unit passes the declared
!     dimension information of the arrays. These values are compared
!     against the problem-dependent needs within the subprogram. If any
!     of the dimensions are too small an error message is printed and a
!     negative value of MODE is returned, -41 to -47. The printed error
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
!        call SBOCLS( )
!
!     Use of this option adds 8 to the required length of IOPT(*).
!
!   3
!   -
!     This option can change the type of scaling for the data matrix.
!     Nominally each nonzero column of the matrix is scaled so that the
!     magnitude of its largest entry is equal to the value ONE. If LP
!     is the processing pointer for IOPT(*),
!
!        IOPT(LP)=3
!        IOPT(LP+1)=1,2 or 3
!            1= Nominal scaling as noted;
!            2= Each nonzero column scaled to have length ONE;
!            3= Identity scaling; scaling effectively suppressed.
!         .
!        call SBOCLS( )
!
!     Use of this option adds 2 to the required length of IOPT(*).
!
!   4
!   -
!     This options allows the user to provide arbitrary (positive)
!     column scaling for the matrix. If LP is the processing pointer
!     for IOPT(*),
!
!        IOPT(LP)=4
!        IOPT(LP+1)=IOFF
!        X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1)
!        = Positive scale factors for cols. of E.
!         .
!        call SBOCLS( )
!
!     Use of this option adds 2 to the required length of IOPT(*)
!     and NCOLS to the required length of X(*).
!
!   5
!   -
!     This option allows the user to provide an option array to the
!     low-level subprogram SBOLS( ). If LP is the processing pointer
!     for IOPT(*),
!
!        IOPT(LP)=5
!        IOPT(LP+1)= Position in IOPT(*) where option array
!                    data for SBOLS( ) begins.
!         .
!        call SBOCLS( )
!
!     Use of this option adds 2 to the required length of IOPT(*).
!
!   6
!   -
!     This option allows the user to provide an option array to the
!     low-level subprogram SBOLSM( ). If LP is the processing pointer
!     for IOPT(*),
!
!        IOPT(LP)=6
!        IOPT(LP+1)= Position in IOPT(*) where option array
!                    data for SBOLSM( ) begins.
!         .
!        call SBOCLS( )
!
!     Use of this option adds 2 to the required length of IOPT(*).
!
!   7
!   -
!     Move the processing pointer (either forward or backward) to the
!     location IOPT(LP+1). The processing pointer moves to locations
!     LP+2 if option number 7 is used with the value -7.  For
!     example to skip over locations 3,...,NCOLS+2,
!
!       IOPT(1)=7
!       IOPT(2)=NCOLS+3
!       (IOPT(I), I=3,...,NCOLS+2 are not defined here.)
!       IOPT(NCOLS+3)=99
!       call SBOCLS( )
!
!     CAUTION: Misuse of this option can yield some very hard-to-find
!     bugs. Use it with care. It is intended to be used for passing
!     option arrays to other subprograms.
!
!   8
!   -
!     This option allows the user to suppress the algorithmic feature
!     of SBOCLS( ) that processes the constraint equations C*X = Y and
!     resolves infeasibilities. The steps normally done are to solve
!     C*X - Y = 0 in a least squares sense using the stated bounds on
!     both X and Y. Then the "reachable" vector Y = C*X is computed
!     using the solution X obtained. Finally the stated bounds for Y are
!     enlarged to include C*X. To suppress the feature:
!
!
!       IOPT(LP)=8
!         .
!       call SBOCLS( )
!
!     Use of this option adds 1 to the required length of IOPT(*).
!
!   9
!   -
!     This option allows the user to suppress the pretriangularizing
!     step of the least squares matrix that is done within SBOCLS( ).
!     This is primarily a means of enhancing the subprogram efficiency
!     and has little effect on accuracy. To suppress the step, set:
!
!       IOPT(LP)=9
!         .
!       call SBOCLS( )
!
!     Use of this option adds 1 to the required length of IOPT(*).
!
!   99
!   --
!     There are no more options to change.
!
!     Only option numbers -99, -9,-8,...,-1, 1,2,...,9, and 99 are
!     permitted. Other values are errors. Options -99,-1,...,-9 mean
!     that the respective options 99,1,...,9 are left at their default
!     values. An example is the option to suppress the preprocessing of
!     constraints:
!
!       IOPT(1)=-8 Option is recognized but not changed
!       IOPT(2)=99
!       call SBOCLS( )
!
!    Error Messages for SBOCLS()
!    ----- -------- --- --------
!
! WARNING in...
! SBOCLS(). THE ROW DIMENSION OF W(,)=(I1) MUST BE  >=  THE NUMBER
! OF EFFECTIVE ROWS=(I2).
!           IN ABOVE MESSAGE, I1=         1
!           IN ABOVE MESSAGE, I2=         2
! ERROR NUMBER =        41
!
! WARNING IN...
! SBOCLS(). THE COLUMN DIMENSION OF W(,)=(I1) MUST BE  >=  NCOLS+
! MCON+1=(I2).
!           IN ABOVE MESSAGE, I1=         2
!           IN ABOVE MESSAGE, I2=         3
! ERROR NUMBER =        42
!
! WARNING IN...
! SBOCLS(). THE DIMENSIONS OF THE ARRAYS BL(),BU(), AND IND()=(I1)
! MUST BE  >=  NCOLS+MCON=(I2).
!           IN ABOVE MESSAGE, I1=         1
!           IN ABOVE MESSAGE, I2=         2
! ERROR NUMBER =        43
!
! WARNING IN...
! SBOCLS(). THE DIMENSION OF X()=(I1) MUST BE
!  >=  THE REQD.LENGTH=(I2).
!           IN ABOVE MESSAGE, I1=         1
!           IN ABOVE MESSAGE, I2=         2
! ERROR NUMBER =        44
!
! WARNING IN...
! SBOCLS(). THE .
! SBOCLS() THE DIMENSION OF IW()=(I1) MUST BE  >=  2*NCOLS+2*MCON=(I2).
!           IN ABOVE MESSAGE, I1=         1
!           IN ABOVE MESSAGE, I2=         4
! ERROR NUMBER =        46
!
! WARNING IN...
! SBOCLS(). THE DIMENSION OF IOPT()=(I1) MUST BE  >=  THE REQD.
! LEN.=(I2).
!           IN ABOVE MESSAGE, I1=        16
!           IN ABOVE MESSAGE, I2=        18
! ERROR NUMBER =        47
!
! WARNING IN...
! SBOCLS(). ISCALE OPTION=(I1) MUST BE 1-3.
!           IN ABOVE MESSAGE, I1=         0
! ERROR NUMBER =        48
!
! WARNING IN...
! SBOCLS(). OFFSET PAST X(NCOLS) (I1) FOR USER-PROVIDED COLUMN SCALING
! MUST BE POSITIVE.
!           IN ABOVE MESSAGE, I1=         0
! ERROR NUMBER =        49
!
! WARNING IN...
! SBOCLS(). EACH PROVIDED COL. SCALE FACTOR MUST BE POSITIVE.
!  COMPONENT (I1) NOW = (R1).
!           IN ABOVE MESSAGE, I1=         1
!           IN ABOVE MESSAGE, R1=    0.
! ERROR NUMBER =        50
!
! WARNING IN...
! SBOCLS(). THE OPTION NUMBER=(I1) IS NOT DEFINED.
!           IN ABOVE MESSAGE, I1=      1001
! ERROR NUMBER =        51
!
! WARNING IN...
! SBOCLS(). NO. OF ROWS=(I1) MUST BE  >=  0 .AND.  <=  MDW-MCON=(I2).
!           IN ABOVE MESSAGE, I1=         2
!           IN ABOVE MESSAGE, I2=         1
! ERROR NUMBER =        52
!
! WARNING IN...
! SBOCLS(). MDW=(I1) MUST BE POSITIVE.
!           IN ABOVE MESSAGE, I1=         0
! ERROR NUMBER =        53
!
! WARNING IN...
! SBOCLS(). MCON=(I1) MUST BE NONNEGATIVE.
!           IN ABOVE MESSAGE, I1=        -1
! ERROR NUMBER =        54
!
! WARNING IN...
! SBOCLS(). NCOLS=(I1) THE NO. OF VARIABLES MUST BE POSITIVE.
!           IN ABOVE MESSAGE, I1=         0
! ERROR NUMBER =        55
!
! WARNING IN...
! SBOCLS(). FOR J=(I1), IND(J)=(I2) MUST BE 1-4.
!           IN ABOVE MESSAGE, I1=         1
!           IN ABOVE MESSAGE, I2=         0
! ERROR NUMBER =        56
!
! WARNING IN...
! SBOCLS(). FOR J=(I1), BOUND BL(J)=(R1) IS  >  BU(J)=(R2).
!           IN ABOVE MESSAGE, I1=         1
!           IN ABOVE MESSAGE, R1=     .1000000000E+01
!           IN ABOVE MESSAGE, R2=    0.
! ERROR NUMBER =        57
!           LINEAR CONSTRAINTS, SNLA REPT. SAND82-1517, AUG. (1982).
!
!***REFERENCES  R. J. Hanson, Linear least squares with bounds and
!                 linear constraints, Report SAND82-1517, Sandia
!                 Laboratories, August 1982.
!***ROUTINES CALLED  R1MACH, SASUM, SBOLS, SCOPY, SDOT, SNRM2, SSCAL,
!                    XERMSG
!***REVISION HISTORY  (YYMMDD)
!   821220  DATE WRITTEN
!   870803  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!   910819  Added variable M for MOUT+MCON in reference to SBOLS.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SBOCLS
!     REVISED 850604-0900
!     REVISED YYMMDD-HHMM
!
!    PURPOSE
!    -------
!     THIS IS THE MAIN SUBPROGRAM THAT SOLVES THE LEAST SQUARES
!     PROBLEM CONSISTING OF LINEAR CONSTRAINTS
!
!              C*X = Y
!
!     AND LEAST SQUARES EQUATIONS
!
!              E*X = F
!
!     IN THIS FORMULATION THE VECTORS X AND Y ARE BOTH UNKNOWNS.
!     FURTHER, X AND Y MAY BOTH HAVE USER-SPECIFIED BOUNDS ON EACH
!     COMPONENT.  THE USER MUST HAVE DIMENSION STATEMENTS OF THE
!     FORM
!
!     DIMENSION W(MDW,NCOLS+MCON+1), BL(NCOLS+MCON),BU(NCOLS+MCON),
!               X(2*(NCOLS+MCON)+2+NX), RW(6*NCOLS+5*MCON)
!
!     INTEGER IND(NCOLS+MCON), IOPT(16+NI), IW(2*(NCOLS+MCON))
!
!     TO CHANGE THIS SUBPROGRAM FROM SINGLE TO DOUBLE PRECISION BEGIN
!     EDITING AT THE CARD 'C++'.
!     CHANGE THIS SUBPROGRAM TO SBOCLS AND THE STRINGS
!     /SDOT/ TO /DDOT/, /SNRM2/ TO /DNRM2/, /SRELPR/ TO /DRELPR/,
!     /R1MACH/ TO /D1MACH/, /E0/ TO /D0/, /SCOPY/ TO /DCOPY/,
!     /SSCAL/ TO /DSCAL/, /SASUM/ TO /DASUM/, /SBOLS/ TO /DBOLS/,
!     /REAL            / TO /DOUBLE PRECISION/.
! ++
  REAL             W(MDW,*),BL(*),BU(*),X(*),RW(*)
  REAL             ANORM, CNORM, ONE, RNORM, RNORMC, SRELPR
  REAL             T, T1, T2, SDOT, SNRM2, WT, ZERO
  REAL             SASUM, R1MACH
!     THIS VARIABLE REMAINS TYPED REAL.
  INTEGER IND(*),IOPT(*),IW(*),JOPT(05)
  LOGICAL CHECKL,FILTER,ACCUM,PRETRI
  CHARACTER*8 XERN1, XERN2
  CHARACTER*16 XERN3, XERN4
  SAVE IGO,ACCUM,CHECKL
  DATA IGO/0/
!***FIRST EXECUTABLE STATEMENT  SBOCLS
  NERR = 0
  MODE = 0
  if (IGO == 0) THEN
!     DO(CHECK VALIDITY OF INPUT DATA)
!     PROCEDURE(CHECK VALIDITY OF INPUT DATA)
!
!     SEE THAT MDW IS  > 0. GROSS CHECK ONLY.
      if (MDW <= 0) THEN
          WRITE (XERN1, '(I8)') MDW
          call XERMSG ('SLATEC', 'SBOCLS', 'MDW = ' // XERN1 // &
             ' MUST BE POSITIVE.', 53, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
          go to 260
      ENDIF
!
!     SEE THAT NUMBER OF CONSTRAINTS IS NONNEGATIVE.
      if (MCON < 0) THEN
          WRITE (XERN1, '(I8)') MCON
          call XERMSG ('SLATEC', 'SBOCLS', 'MCON = ' // XERN1 // &
             ' MUST BE NON-NEGATIVE', 54, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
          go to 260
      ENDIF
!
!     SEE THAT NUMBER OF UNKNOWNS IS POSITIVE.
      if (NCOLS <= 0) THEN
          WRITE (XERN1, '(I8)') NCOLS
          call XERMSG ('SLATEC', 'SBOCLS', 'NCOLS = ' // XERN1 // &
             ' THE NO. OF VARIABLES, MUST BE POSITIVE.', 55, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
          go to 260
      ENDIF
!
!     SEE THAT CONSTRAINT INDICATORS ARE ALL WELL-DEFINED.
      DO 10 J = 1,NCOLS + MCON
          if (IND(J) < 1 .OR. IND(J) > 4) THEN
              WRITE (XERN1, '(I8)') J
              WRITE (XERN2, '(I8)') IND(J)
              call XERMSG ('SLATEC', 'SBOCLS', &
                'IND(' // XERN1 // ') = ' // XERN2 // &
                ' MUST BE 1-4.', 56, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 260
          ENDIF
   10     CONTINUE
!
!     SEE THAT BOUNDS ARE CONSISTENT.
      DO 20 J = 1,NCOLS + MCON
          if (IND(J) == 3) THEN
              if (BL(J) > BU(J)) THEN
                 WRITE (XERN1, '(I8)') J
                 WRITE (XERN3, '(1PE15.6)') BL(J)
                 WRITE (XERN4, '(1PE15.6)') BU(J)
                 call XERMSG ('SLATEC', 'SBOCLS', &
                    'BOUND BL(' // XERN1 // ') = ' // XERN3 // &
                    ' IS  >  BU(' // XERN1 // ') = ' // XERN4, &
                    57, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
                  go to 260
              ENDIF
          ENDIF
   20     CONTINUE
!     END PROCEDURE
!     DO(PROCESS OPTION ARRAY)
!     PROCEDURE(PROCESS OPTION ARRAY)
      ZERO = 0.E0
      ONE = 1.E0
      SRELPR = R1MACH(4)
      CHECKL = .FALSE.
      FILTER = .TRUE.
      LENX = 2* (NCOLS+MCON) + 2
      ISCALE = 1
      IGO = 1
      ACCUM = .FALSE.
      PRETRI = .TRUE.
      LOPT = 0
      MOPT = 0
      LP = 0
      LDS = 0
!     DO FOREVER
   30     CONTINUE
      LP = LP + LDS
      IP = IOPT(LP+1)
      JP = ABS(IP)
!
!     TEST FOR NO MORE OPTIONS TO CHANGE.
      if (IP == 99) THEN
          if (LOPT == 0) LOPT = - (LP+2)
          if (MOPT == 0) MOPT = - (ABS(LOPT)+7)
          if (LOPT < 0) THEN
              LBOU = ABS(LOPT)
          ELSE
              LBOU = LOPT - 15
          ENDIF
!
!     SEND COL. SCALING TO SBOLS().
          IOPT(LBOU) = 4
          IOPT(LBOU+1) = 1
!
!     PASS AN OPTION ARRAY FOR SBOLSM().
          IOPT(LBOU+2) = 5
!
!     LOC. OF OPTION ARRAY FOR SBOLSM( ).
          IOPT(LBOU+3) = 8
!
!     SKIP TO START OF USER-GIVEN OPTION ARRAY FOR SBOLS().
          IOPT(LBOU+4) = 6
          IOPT(LBOU+6) = 99
          if (LOPT > 0) THEN
              IOPT(LBOU+5) = LOPT - LBOU + 1
          ELSE
              IOPT(LBOU+4) = -IOPT(LBOU+4)
          ENDIF
          if (MOPT < 0) THEN
              LBOUM = ABS(MOPT)
          ELSE
              LBOUM = MOPT - 8
          ENDIF
!
!     CHANGE PRETRIANGULARIZATION FACTOR IN SBOLSM().
          IOPT(LBOUM) = 5
          IOPT(LBOUM+1) = NCOLS + MCON + 1
!
!     PASS WEIGHT TO SBOLSM() FOR RANK TEST.
          IOPT(LBOUM+2) = 6
          IOPT(LBOUM+3) = NCOLS + MCON + 2
          IOPT(LBOUM+4) = MCON
!
!     SKIP TO USER-GIVEN OPTION ARRAY FOR SBOLSM( ).
          IOPT(LBOUM+5) = 1
          IOPT(LBOUM+7) = 99
          if (MOPT > 0) THEN
              IOPT(LBOUM+6) = MOPT - LBOUM + 1
          ELSE
              IOPT(LBOUM+5) = -IOPT(LBOUM+5)
          ENDIF
!     EXIT FOREVER
          go to 50
      ELSE if (JP == 99) THEN
          LDS = 1
!     CYCLE FOREVER
          go to 50
      ELSE if (JP == 1) THEN
          if (IP > 0) THEN
!
!     SET UP DIRECTION FLAG LOCATION, ROW STACKING POINTER
!     LOCATION, AND LOCATION FOR NUMBER OF NEW ROWS.
              LOCACC = LP + 2
!
!                  IOPT(LOCACC-1)=OPTION NUMBER FOR SEQ. ACCUMULATION.
!     CONTENTS..   IOPT(LOCACC  )=USER DIRECTION FLAG, 1 OR 2.
!                  IOPT(LOCACC+1)=ROW STACKING POINTER.
!                  IOPT(LOCACC+2)=NUMBER OF NEW ROWS TO PROCESS.
!     USER ACTION WITH THIS OPTION..
!      (SET UP OPTION DATA FOR SEQ. ACCUMULATION IN IOPT(*).)
!      (MOVE BLOCK OF EQUATIONS INTO W(*,*)  STARTING AT FIRST
!       ROW OF W(*,*) BELOW THE ROWS FOR THE CONSTRAINT MATRIX C.
!       SET IOPT(LOCACC+2)=NO. OF LEAST SQUARES EQUATIONS IN BLOCK.
!              LOOP
!              call SBOCLS()
!
!                  if ( IOPT(LOCACC)  ==  1) THEN
!                      STACK EQUAS. INTO W(*,*), STARTING AT
!                      ROW IOPT(LOCACC+1).
!                       INTO W(*,*).
!                       SET IOPT(LOCACC+2)=NO. OF EQUAS.
!                      if LAST BLOCK OF EQUAS., SET IOPT(LOCACC)=2.
!                  ELSE if IOPT(LOCACC)  ==  2) THEN
!                      (PROCESS IS OVER. EXIT LOOP.)
!                  ELSE
!                      (ERROR CONDITION. SHOULD NOT HAPPEN.)
!                  end if
!              END LOOP
              IOPT(LOCACC+1) = MCON + 1
              ACCUM = .TRUE.
              IOPT(LOCACC) = IGO
          ENDIF
          LDS = 4
!     CYCLE FOREVER
          go to 30
      ELSE if (JP == 2) THEN
          if (IP > 0) THEN
!
!     GET ACTUAL LENGTHS OF ARRAYS FOR CHECKING AGAINST NEEDS.
              LOCDIM = LP + 2
!
!     LMDW >= MCON+MAX(MOUT,NCOLS), if MCON > 0 .AND FILTER
!     LMDW >= MCON+MOUT, OTHERWISE
!
!     LNDW >= NCOLS+MCON+1
!     LLB  >= NCOLS+MCON
!     LLX  >= 2*(NCOLS+MCON)+2+EXTRA REQD. IN OPTIONS.
!     LLRW >= 6*NCOLS+5*MCON
!     LLIW >= 2*(NCOLS+MCON)
!     LIOP >=  AMOUNT REQD. FOR OPTION ARRAY.
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
!     CYCLE FOREVER
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
                  call XERMSG ('SLATEC', 'SBOCLS', &
                     'ISCALE OPTION = ' // XERN1 // ' MUST BE 1-3', &
                     48, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
                  go to 260
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
                  call XERMSG ('SLATEC', 'SBOCLS', &
                     'OFFSET PAST X(NCOLS) (' // XERN1 // &
             ') FOR USER-PROVIDED COLUMN SCALING MUST BE POSITIVE.', &
                     49, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
                  go to 260
              ENDIF
              call SCOPY(NCOLS,X(NCOLS+IOPT(LP+2)),1,RW,1)
              LENX = LENX + NCOLS
              DO 40 J = 1,NCOLS
                  if (RW(J) <= ZERO) THEN
                      WRITE (XERN1, '(I8)') J
                      WRITE (XERN3, '(1PE15.6)') RW(J)
                      call XERMSG ('SLATEC', 'SBOCLS', &
                         'EACH PROVIDED COLUMN SCALE FACTOR ' // &
                         'MUST BE POSITIVE.$$COMPONENT ' // XERN1 // &
                         ' NOW = ' // XERN3, 50, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
                      go to 260
                  ENDIF
   40             CONTINUE
          ENDIF
          LDS = 2
!     CYCLE FOREVER
          go to 30
!
!     IN THIS OPTION AN OPTION ARRAY IS PROVIDED TO SBOLS().
      ELSE if (JP == 5) THEN
          if (IP > 0) THEN
              LOPT = IOPT(LP+2)
          ENDIF
          LDS = 2
!     CYCLE FOREVER
          go to 30
!
!     IN THIS OPTION AN OPTION ARRAY IS PROVIDED TO SBOLSM().
      ELSE if (JP == 6) THEN
          if (IP > 0) THEN
              MOPT = IOPT(LP+2)
          ENDIF
          LDS = 2
!     CYCLE FOREVER
          go to 30
!
!     THIS OPTION USES THE NEXT LOC OF IOPT(*) AS A
!     POINTER VALUE TO SKIP TO NEXT.
      ELSE if (JP == 7) THEN
          if (IP > 0) THEN
              LP = IOPT(LP+2) - 1
              LDS = 0
          ELSE
              LDS = 2
          ENDIF
!     CYCLE FOREVER
          go to 30
!
!     THIS OPTION AVOIDS THE CONSTRAINT RESOLVING PHASE FOR
!     THE LINEAR CONSTRAINTS C*X=Y.
      ELSE if (JP == 8) THEN
          FILTER = .NOT. (IP > 0)
          LDS = 1
!     CYCLE FOREVER
          go to 30
!
!     THIS OPTION SUPPRESSES PRE-TRIANGULARIZATION OF THE LEAST
!     SQUARES EQUATIONS.
      ELSE if (JP == 9) THEN
          PRETRI = .NOT. (IP > 0)
          LDS = 1
!     CYCLE FOREVER
          go to 30
!
!     NO VALID OPTION NUMBER WAS NOTED. THIS IS AN ERROR CONDITION.
      ELSE
          WRITE (XERN1, '(I8)') JP
          call XERMSG ('SLATEC', 'SBOCLS', 'OPTION NUMBER = ' // &
             XERN1 // ' IS NOT DEFINED.', 51, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
          go to 260
      ENDIF
!     END FOREVER
!     END PROCEDURE
   50     CONTINUE
      if (CHECKL) THEN
!     DO(CHECK LENGTHS OF ARRAYS)
!     PROCEDURE(CHECK LENGTHS OF ARRAYS)
!
!     THIS FEATURE ALLOWS THE USER TO MAKE SURE THAT THE
!     ARRAYS ARE LONG ENOUGH FOR THE INTENDED PROBLEM SIZE AND USE.
       if ( FILTER .AND. .NOT.ACCUM) THEN
            MDWL=MCON+MAX(MROWS,NCOLS)
       ELSE
            MDWL=MCON+NCOLS+1
       ENDIF
          if (LMDW < MDWL) THEN
              WRITE (XERN1, '(I8)') LMDW
              WRITE (XERN2, '(I8)') MDWL
              call XERMSG ('SLATEC', 'SBOCLS', &
                 'THE ROW DIMENSION OF W(,) = ' // XERN1 // &
                 ' MUST BE  >=  THE NUMBER OF EFFECTIVE ROWS = ' // &
                 XERN2, 41, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 260
          ENDIF
          if (LNDW < NCOLS+MCON+1) THEN
              WRITE (XERN1, '(I8)') LNDW
              WRITE (XERN2, '(I8)') NCOLS+MCON+1
              call XERMSG ('SLATEC', 'SBOCLS', &
                 'THE COLUMN DIMENSION OF W(,) = ' // XERN1 // &
                 ' MUST BE  >=  NCOLS+MCON+1 = ' // XERN2, 42, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 260
          ENDIF
          if (LLB < NCOLS+MCON) THEN
              WRITE (XERN1, '(I8)') LLB
              WRITE (XERN2, '(I8)') NCOLS+MCON
              call XERMSG ('SLATEC', 'SBOCLS', &
             'THE DIMENSIONS OF THE ARRAYS BS(), BU(), AND IND() = ' &
                 // XERN1 // ' MUST BE  >=  NCOLS+MCON = ' // XERN2, &
                 43, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 260
          ENDIF
          if (LLX < LENX) THEN
              WRITE (XERN1, '(I8)') LLX
              WRITE (XERN2, '(I8)') LENX
              call XERMSG ('SLATEC', 'SBOCLS', &
                'THE DIMENSION OF X() = ' // XERN1 // &
                ' MUST BE  >=  THE REQD. LENGTH = ' // XERN2, 44, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 260
          ENDIF
          if (LLRW < 6*NCOLS+5*MCON) THEN
              WRITE (XERN1, '(I8)') LLRW
              WRITE (XERN2, '(I8)') 6*NCOLS+5*MCON
              call XERMSG ('SLATEC', 'SBOCLS', &
                 'THE DIMENSION OF RW() = ' // XERN1 // &
                 ' MUST BE  >=  6*NCOLS+5*MCON = ' // XERN2, 45, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 260
          ENDIF
          if (LLIW < 2*NCOLS+2*MCON) THEN
              WRITE (XERN1, '(I8)') LLIW
              WRITE (XERN2, '(I8)') 2*NCOLS+2*MCON
              call XERMSG ('SLATEC', 'SBOCLS', &
                 'THE DIMENSION OF IW() = ' // XERN1 // &
                 ' MUST BE  >=  2*NCOLS+2*MCON = ' // XERN2, 46, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 260
          ENDIF
          if (LIOPT < LP+17) THEN
              WRITE (XERN1, '(I8)') LIOPT
              WRITE (XERN2, '(I8)') LP+17
              call XERMSG ('SLATEC', 'SBOCLS', &
                 'THE DIMENSION OF IOPT() = ' // XERN1 // &
                 ' MUST BE  >=  THE REQD. LEN = ' // XERN2, 47, 1)
!     DO(RETURN TO USER PROGRAM UNIT)
              go to 260
          ENDIF
!     END PROCEDURE
      ENDIF
  end if
!
!     OPTIONALLY GO BACK TO THE USER FOR ACCUMULATION OF LEAST SQUARES
!     EQUATIONS AND DIRECTIONS FOR PROCESSING THESE EQUATIONS.
!     DO(ACCUMULATE LEAST SQUARES EQUATIONS)
!     PROCEDURE(ACCUMULATE LEAST SQUARES EQUATIONS)
  if (ACCUM) THEN
      MROWS = IOPT(LOCACC+1) - 1 - MCON
      INROWS = IOPT(LOCACC+2)
      MNEW = MROWS + INROWS
      if (MNEW < 0 .OR. MNEW+MCON > MDW) THEN
          WRITE (XERN1, '(I8)') MNEW
          WRITE (XERN2, '(I8)') MDW-MCON
          call XERMSG ('SLATEC', 'SBOCLS', 'NO. OF ROWS = ' // &
             XERN1 //  ' MUST BE  >=  0 .AND.  <=  MDW-MCON = ' // &
             XERN2, 52, 1)
!    (RETURN TO USER PROGRAM UNIT)
          go to 260
      ENDIF
  end if
!
!     USE THE SOFTWARE OF SBOLS( ) FOR THE TRIANGULARIZATION OF THE
!     LEAST SQUARES MATRIX.  THIS MAY INVOLVE A SYSTALTIC INTERCHANGE
!     OF PROCESSING POINTERS BETWEEN THE CALLING AND CALLED (SBOLS())
!     PROGRAM UNITS.
  JOPT(01) = 1
  JOPT(02) = 2
  JOPT(04) = MROWS
  JOPT(05) = 99
  IRW = NCOLS + 1
  IIW = 1
  if (ACCUM .OR. PRETRI) THEN
      call SBOLS(W(MCON+1,1),MDW,MOUT,NCOLS,BL,BU,IND,JOPT,X,RNORM, &
                 MODE,RW(IRW),IW(IIW))
  ELSE
      MOUT = MROWS
  end if
  if (ACCUM) THEN
      ACCUM = IOPT(LOCACC)  ==  1
      IOPT(LOCACC+1) = JOPT(03) + MCON
      MROWS = MIN(NCOLS+1,MNEW)
  end if
!     END PROCEDURE
  if (ACCUM) RETURN
!     DO(SOLVE CONSTRAINED AND BOUNDED LEAST SQUARES PROBLEM)
!     PROCEDURE(SOLVE CONSTRAINED AND BOUNDED LEAST SQUARES PROBLEM)
!
!     MOVE RIGHT HAND SIDE OF LEAST SQUARES EQUATIONS.
  call SCOPY(MOUT,W(MCON+1,NCOLS+1),1,W(MCON+1,NCOLS+MCON+1),1)
  if (MCON > 0 .AND. FILTER) THEN
!
!     PROJECT THE LINEAR CONSTRAINTS INTO A REACHABLE SET.
      DO 60 I = 1,MCON
          call SCOPY(NCOLS,W(I,1),MDW,W(MCON+1,NCOLS+I),1)
   60     CONTINUE
!
!      PLACE (-)IDENTITY MATRIX AFTER CONSTRAINT DATA.
      DO 70 J = NCOLS + 1,NCOLS + MCON + 1
          W(1,J) = ZERO
          call SCOPY(MCON,W(1,J),0,W(1,J),1)
   70     CONTINUE
      W(1,NCOLS+1) = -ONE
      call SCOPY(MCON,W(1,NCOLS+1),0,W(1,NCOLS+1),MDW+1)
!
!     OBTAIN A 'FEASIBLE POINT' FOR THE LINEAR CONSTRAINTS.
      JOPT(01) = 99
      IRW = NCOLS + 1
      IIW = 1
      call SBOLS(W,MDW,MCON,NCOLS+MCON,BL,BU,IND,JOPT,X,RNORMC, &
                 MODEC,RW(IRW),IW(IIW))
!
!     ENLARGE THE BOUNDS SET, if REQUIRED, TO INCLUDE POINTS THAT
!     CAN BE REACHED.
      DO 130 J = NCOLS + 1,NCOLS + MCON
          ICASE = IND(J)
          if (ICASE < 4) THEN
              T = SDOT(NCOLS,W(MCON+1,J),1,X,1)
          ENDIF
          go to (80,90,100,110),ICASE
          go to 120
!     CASE 1
   80         BL(J) = MIN(T,BL(J))
          go to 120
!     CASE 2
   90         BU(J) = MAX(T,BU(J))
          go to 120
!     CASE 3
  100         BL(J) = MIN(T,BL(J))
          BU(J) = MAX(T,BU(J))
          go to 120
!     CASE 4
  110         CONTINUE
  120         CONTINUE
  130     CONTINUE
!
!     MOVE CONSTRAINT DATA BACK TO THE ORIGINAL AREA.
      DO 140 J = NCOLS + 1,NCOLS + MCON
          call SCOPY(NCOLS,W(MCON+1,J),1,W(J-NCOLS,1),MDW)
  140     CONTINUE
  end if
  if (MCON > 0) THEN
      DO 150 J = NCOLS + 1,NCOLS + MCON
          W(MCON+1,J) = ZERO
          call SCOPY(MOUT,W(MCON+1,J),0,W(MCON+1,J),1)
  150     CONTINUE
!
!     PUT IN (-)IDENTITY MATRIX (POSSIBLY) ONCE AGAIN.
      DO 160 J = NCOLS + 1,NCOLS + MCON + 1
          W(1,J) = ZERO
          call SCOPY(MCON,W(1,J),0,W(1,J),1)
  160     CONTINUE
      W(1,NCOLS+1) = -ONE
      call SCOPY(MCON,W(1,NCOLS+1),0,W(1,NCOLS+1),MDW+1)
  end if
!
!     COMPUTE NOMINAL COLUMN SCALING FOR THE UNWEIGHTED MATRIX.
  CNORM = ZERO
  ANORM = ZERO
  DO 170 J = 1,NCOLS
      T1 = SASUM(MCON,W(1,J),1)
      T2 = SASUM(MOUT,W(MCON+1,1),1)
      T = T1 + T2
      if (T == ZERO) T = ONE
      CNORM = MAX(CNORM,T1)
      ANORM = MAX(ANORM,T2)
      X(NCOLS+MCON+J) = ONE/T
  170 CONTINUE
  go to (180,190,210,220),ISCALE
  go to 230
!     CASE 1
  180 CONTINUE
  go to 230
!     CASE 2
!
!     SCALE COLS. (BEFORE WEIGHTING) TO HAVE LENGTH ONE.
  190 DO 200 J = 1,NCOLS
      T = SNRM2(MCON+MOUT,W(1,J),1)
      if (T == ZERO) T = ONE
      X(NCOLS+MCON+J) = ONE/T
  200 CONTINUE
  go to 230
!     CASE 3
!
!     SUPPRESS SCALING (USE UNIT MATRIX).
  210 X(NCOLS+MCON+1) = ONE
  call SCOPY(NCOLS,X(NCOLS+MCON+1),0,X(NCOLS+MCON+1),1)
  go to 230
!     CASE 4
!
!     THE USER HAS PROVIDED SCALING.
  220 call SCOPY(NCOLS,RW,1,X(NCOLS+MCON+1),1)
  230 CONTINUE
  DO 240 J = NCOLS + 1,NCOLS + MCON
      X(NCOLS+MCON+J) = ONE
  240 CONTINUE
!
!     WEIGHT THE LEAST SQUARES EQUATIONS.
  WT = SRELPR
  if (ANORM > ZERO) WT = WT/ANORM
  if (CNORM > ZERO) WT = WT*CNORM
  DO 250 I = 1,MOUT
      call SSCAL(NCOLS,WT,W(I+MCON,1),MDW)
  250 CONTINUE
  call SSCAL(MOUT,WT,W(MCON+1,MCON+NCOLS+1),1)
  LRW = 1
  LIW = 1
!
!     SET THE NEW TRIANGULARIZATION FACTOR.
  X(2* (NCOLS+MCON)+1) = ZERO
!
!     SET THE WEIGHT TO USE IN COMPONENTS  >  MCON,
!     WHEN MAKING LINEAR INDEPENDENCE TEST.
  X(2* (NCOLS+MCON)+2) = ONE/WT
  M=MOUT+MCON
  call SBOLS(W,MDW,M,NCOLS+MCON,BL,BU,IND,IOPT(LBOU),X, &
             RNORM,MODE,RW(LRW),IW(LIW))
  RNORM = RNORM/WT
!     END PROCEDURE
!     PROCEDURE(RETURN TO USER PROGRAM UNIT)
  260 if ( MODE >= 0)MODE = -NERR
  IGO = 0
  return
!     END PROGRAM
end
