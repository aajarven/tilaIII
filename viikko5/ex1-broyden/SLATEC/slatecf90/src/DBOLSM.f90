subroutine DBOLSM (W, MDW, MINPUT, NCOLS, BL, BU, IND, IOPT, X, &
     RNORM, MODE, RW, WW, SCL, IBASIS, IBB)
!
!! DBOLSM is subsidiary to DBOCLS and DBOLS.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SBOLSM-S, DBOLSM-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!            **** Double Precision Version of SBOLSM ****
!   **** All INPUT and OUTPUT real variables are DOUBLE PRECISION ****
!
!          Solve E*X = F (least squares sense) with bounds on
!            selected X values.
!     The user must have DIMENSION statements of the form:
!
!       DIMENSION W(MDW,NCOLS+1), BL(NCOLS), BU(NCOLS),
!      * X(NCOLS+NX), RW(NCOLS), WW(NCOLS), SCL(NCOLS)
!       INTEGER IND(NCOLS), IOPT(1+NI), IBASIS(NCOLS), IBB(NCOLS)
!
!     (Here NX=number of extra locations required for options 1,...,7;
!     NX=0 for no options; here NI=number of extra locations possibly
!     required for options 1-7; NI=0 for no options; NI=14 if all the
!     options are simultaneously in use.)
!
!    INPUT
!    -----
!
!    --------------------
!    W(MDW,*),MINPUT,NCOLS
!    --------------------
!     The array W(*,*) contains the matrix [E:F] on entry. The matrix
!     [E:F] has MINPUT rows and NCOLS+1 columns. This data is placed in
!     the array W(*,*) with E occupying the first NCOLS columns and the
!     right side vector F in column NCOLS+1. The row dimension, MDW, of
!     the array W(*,*) must satisfy the inequality MDW .ge. MINPUT.
!     Other values of MDW are errors. The values of MINPUT and NCOLS
!     must be positive. Other values are errors.
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
!    2.    For IND(J)=2, require X(J) .le. BU(J).
!    3.    For IND(J)=3, require X(J) .ge. BL(J) and
!                                X(J) .le. BU(J).
!    4.    For IND(J)=4, no bounds on X(J) are required.
!     The values of BL(*),BL(*) are modified by the subprogram. Values
!     other than 1,2,3 or 4 for IND(J) are errors. In the case IND(J)=3
!     (upper and lower bounds) the condition BL(J) .gt. BU(J) is an
!     error.
!
!    -------
!    IOPT(*)
!    -------
!     This is the array where the user can specify nonstandard options
!     for DBOLSM. Most of the time this feature can be ignored by
!     setting the input value IOPT(1)=99. Occasionally users may have
!     needs that require use of the following subprogram options. For
!     details about how to use the options see below: IOPT(*) CONTENTS.
!
!     Option Number   Brief Statement of Purpose
!     ----- ------   ----- --------- -- -------
!           1         Move the IOPT(*) processing pointer.
!           2         Change rank determination tolerance.
!           3         Change blow-up factor that determines the
!                     size of variables being dropped from active
!                     status.
!           4         Reset the maximum number of iterations to use
!                     in solving the problem.
!           5         The data matrix is triangularized before the
!                     problem is solved whenever (NCOLS/MINPUT) .lt.
!                     FAC. Change the value of FAC.
!           6         Redefine the weighting matrix used for
!                     linear independence checking.
!           7         Debug output is desired.
!          99         No more options to change.
!
!    ----
!    X(*)
!    ----
!     This array is used to pass data associated with options 1,2,3 and
!     5. Ignore this input parameter if none of these options are used.
!     Otherwise see below: IOPT(*) CONTENTS.
!
!    ----------------
!    IBASIS(*),IBB(*)
!    ----------------
!     These arrays must be initialized by the user. The values
!         IBASIS(J)=J, J=1,...,NCOLS
!         IBB(J)   =1, J=1,...,NCOLS
!     are appropriate except when using nonstandard features.
!
!    ------
!    SCL(*)
!    ------
!     This is the array of scaling factors to use on the columns of the
!     matrix E. These values must be defined by the user. To suppress
!     any column scaling set SCL(J)=1.0, J=1,...,NCOLS.
!
!    OUTPUT
!    ------
!
!    ----------
!    X(*),RNORM
!    ----------
!     The array X(*) contains a solution (if MODE .ge. 0 or .eq. -22)
!     for the constrained least squares problem. The value RNORM is the
!     minimum residual vector length.
!
!    ----
!    MODE
!    ----
!     The sign of mode determines whether the subprogram has completed
!     normally, or encountered an error condition or abnormal status.
!     A value of MODE .ge. 0 signifies that the subprogram has completed
!     normally. The value of MODE (.ge. 0) is the number of variables
!     in an active status: not at a bound nor at the value ZERO, for
!     the case of free variables. A negative value of MODE will be one
!     of the 18 cases -38,-37,...,-22, or -1. Values .lt. -1 correspond
!     to an abnormal completion of the subprogram. To understand the
!     abnormal completion codes see below: ERROR MESSAGES for DBOLSM
!     An approximate solution will be returned to the user only when
!     maximum iterations is reached, MODE=-22.
!
!    -----------
!    RW(*),WW(*)
!    -----------
!     These are working arrays each with NCOLS entries. The array RW(*)
!     contains the working (scaled, nonactive) solution values. The
!     array WW(*) contains the working (scaled, active) gradient vector
!     values.
!
!    ----------------
!    IBASIS(*),IBB(*)
!    ----------------
!     These arrays contain information about the status of the solution
!     when MODE .ge. 0. The indices IBASIS(K), K=1,...,MODE, show the
!     nonactive variables; indices IBASIS(K), K=MODE+1,..., NCOLS are
!     the active variables. The value (IBB(J)-1) is the number of times
!     variable J was reflected from its upper bound. (Normally the user
!     can ignore these parameters.)
!
!    IOPT(*) CONTENTS
!    ------- --------
!     The option array allows a user to modify internal variables in
!     the subprogram without recompiling the source code. A central
!     goal of the initial software design was to do a good job for most
!     people. Thus the use of options will be restricted to a select
!     group of users. The processing of the option array proceeds as
!     follows: a pointer, here called LP, is initially set to the value
!     1. The value is updated as the options are processed.  At the
!     pointer position the option number is extracted and used for
!     locating other information that allows for options to be changed.
!     The portion of the array IOPT(*) that is used for each option is
!     fixed; the user and the subprogram both know how many locations
!     are needed for each option. A great deal of error checking is
!     done by the subprogram on the contents of the option array.
!     Nevertheless it is still possible to give the subprogram optional
!     input that is meaningless. For example, some of the options use
!     the location X(NCOLS+IOFF) for passing data. The user must manage
!     the allocation of these locations when more than one piece of
!     option data is being passed to the subprogram.
!
!   1
!   -
!     Move the processing pointer (either forward or backward) to the
!     location IOPT(LP+1). The processing pointer is moved to location
!     LP+2 of IOPT(*) in case IOPT(LP)=-1.  For example to skip over
!     locations 3,...,NCOLS+2 of IOPT(*),
!
!       IOPT(1)=1
!       IOPT(2)=NCOLS+3
!       (IOPT(I), I=3,...,NCOLS+2 are not defined here.)
!       IOPT(NCOLS+3)=99
!       call DBOLSM
!
!     CAUTION: Misuse of this option can yield some very hard-to-find
!     bugs.  Use it with care.
!
!   2
!   -
!     The algorithm that solves the bounded least squares problem
!     iteratively drops columns from the active set. This has the
!     effect of joining a new column vector to the QR factorization of
!     the rectangular matrix consisting of the partially triangularized
!     nonactive columns. After triangularizing this matrix a test is
!     made on the size of the pivot element. The column vector is
!     rejected as dependent if the magnitude of the pivot element is
!     .le. TOL* magnitude of the column in components strictly above
!     the pivot element. Nominally the value of this (rank) tolerance
!     is TOL = SQRT(R1MACH(4)). To change only the value of TOL, for
!     example,
!
!       X(NCOLS+1)=TOL
!       IOPT(1)=2
!       IOPT(2)=1
!       IOPT(3)=99
!       call DBOLSM
!
!     Generally, if LP is the processing pointer for IOPT(*),
!
!       X(NCOLS+IOFF)=TOL
!       IOPT(LP)=2
!       IOPT(LP+1)=IOFF
!        .
!       call DBOLSM
!
!     The required length of IOPT(*) is increased by 2 if option 2 is
!     used; The required length of X(*) is increased by 1. A value of
!     IOFF .le. 0 is an error. A value of TOL .le. R1MACH(4) gives a
!     warning message; it is not considered an error.
!
!   3
!   -
!     A solution component is left active (not used) if, roughly
!     speaking, it seems too large. Mathematically the new component is
!     left active if the magnitude is .ge.((vector norm of F)/(matrix
!     norm of E))/BLOWUP. Nominally the factor BLOWUP = SQRT(R1MACH(4)).
!     To change only the value of BLOWUP, for example,
!
!       X(NCOLS+2)=BLOWUP
!       IOPT(1)=3
!       IOPT(2)=2
!       IOPT(3)=99
!       call DBOLSM
!
!     Generally, if LP is the processing pointer for IOPT(*),
!
!       X(NCOLS+IOFF)=BLOWUP
!       IOPT(LP)=3
!       IOPT(LP+1)=IOFF
!        .
!       call DBOLSM
!
!     The required length of IOPT(*) is increased by 2 if option 3 is
!     used; the required length of X(*) is increased by 1. A value of
!     IOFF .le. 0 is an error. A value of BLOWUP .le. 0.0 is an error.
!
!   4
!   -
!     Normally the algorithm for solving the bounded least squares
!     problem requires between NCOLS/3 and NCOLS drop-add steps to
!     converge. (this remark is based on examining a small number of
!     test cases.) The amount of arithmetic for such problems is
!     typically about twice that required for linear least squares if
!     there are no bounds and if plane rotations are used in the
!     solution method. Convergence of the algorithm, while
!     mathematically certain, can be much slower than indicated. To
!     avoid this potential but unlikely event ITMAX drop-add steps are
!     permitted. Nominally ITMAX=5*(MAX(MINPUT,NCOLS)). To change the
!     value of ITMAX, for example,
!
!       IOPT(1)=4
!       IOPT(2)=ITMAX
!       IOPT(3)=99
!       call DBOLSM
!
!     Generally, if LP is the processing pointer for IOPT(*),
!
!       IOPT(LP)=4
!       IOPT(LP+1)=ITMAX
!        .
!       call DBOLSM
!
!     The value of ITMAX must be .gt. 0. Other values are errors. Use
!     of this option increases the required length of IOPT(*) by 2.
!
!   5
!   -
!     For purposes of increased efficiency the MINPUT by NCOLS+1 data
!     matrix [E:F] is triangularized as a first step whenever MINPUT
!     satisfies FAC*MINPUT .gt. NCOLS. Nominally FAC=0.75. To change the
!     value of FAC,
!
!       X(NCOLS+3)=FAC
!       IOPT(1)=5
!       IOPT(2)=3
!       IOPT(3)=99
!       call DBOLSM
!
!     Generally, if LP is the processing pointer for IOPT(*),
!
!       X(NCOLS+IOFF)=FAC
!       IOPT(LP)=5
!       IOPT(LP+1)=IOFF
!        .
!       call DBOLSM
!
!     The value of FAC must be nonnegative. Other values are errors.
!     Resetting FAC=0.0 suppresses the initial triangularization step.
!     Use of this option increases the required length of IOPT(*) by 2;
!     The required length of of X(*) is increased by 1.
!
!   6
!   -
!     The norm used in testing the magnitudes of the pivot element
!     compared to the mass of the column above the pivot line can be
!     changed. The type of change that this option allows is to weight
!     the components with an index larger than MVAL by the parameter
!     WT. Normally MVAL=0 and WT=1. To change both the values MVAL and
!     WT, where LP is the processing pointer for IOPT(*),
!
!       X(NCOLS+IOFF)=WT
!       IOPT(LP)=6
!       IOPT(LP+1)=IOFF
!       IOPT(LP+2)=MVAL
!
!     Use of this option increases the required length of IOPT(*) by 3.
!     The length of X(*) is increased by 1. Values of MVAL must be
!     nonnegative and not greater than MINPUT. Other values are errors.
!     The value of WT must be positive. Any other value is an error. If
!     either error condition is present a message will be printed.
!
!   7
!   -
!     Debug output, showing the detailed add-drop steps for the
!     constrained least squares problem, is desired. This option is
!     intended to be used to locate suspected bugs.
!
!   99
!   --
!     There are no more options to change.
!
!     The values for options are 1,...,7,99, and are the only ones
!     permitted. Other values are errors. Options -99,-1,...,-7 mean
!     that the repective options 99,1,...,7 are left at their default
!     values. An example is the option to modify the (rank) tolerance:
!
!       X(NCOLS+1)=TOL
!       IOPT(1)=-2
!       IOPT(2)=1
!       IOPT(3)=99
!
!    Error Messages for DBOLSM
!    ----- -------- --- ---------
!    -22    MORE THAN ITMAX = ... ITERATIONS SOLVING BOUNDED LEAST
!           SQUARES PROBLEM.
!
!    -23    THE OPTION NUMBER = ... IS NOT DEFINED.
!
!    -24    THE OFFSET = ... BEYOND POSTION NCOLS = ... MUST BE POSITIVE
!           FOR OPTION NUMBER 2.
!
!    -25    THE TOLERANCE FOR RANK DETERMINATION = ... IS LESS THAN
!           MACHINE PRECISION = ....
!
!    -26    THE OFFSET = ... BEYOND POSITION NCOLS = ... MUST BE POSTIVE
!           FOR OPTION NUMBER 3.
!
!    -27    THE RECIPROCAL OF THE BLOW-UP FACTOR FOR REJECTING VARIABLES
!           MUST BE POSITIVE. NOW = ....
!
!    -28    THE MAXIMUM NUMBER OF ITERATIONS = ... MUST BE POSITIVE.
!
!    -29    THE OFFSET = ... BEYOND POSITION NCOLS = ... MUST BE POSTIVE
!           FOR OPTION NUMBER 5.
!
!    -30    THE FACTOR (NCOLS/MINPUT) WHERE PRETRIANGULARIZING IS
!           PERFORMED MUST BE NONNEGATIVE. NOW = ....
!
!    -31    THE NUMBER OF ROWS = ... MUST BE POSITIVE.
!
!    -32    THE NUMBER OF COLUMNS = ... MUST BE POSTIVE.
!
!    -33    THE ROW DIMENSION OF W(,) = ... MUST BE  >=  THE NUMBER OF
!           ROWS = ....
!
!    -34    FOR J = ... THE CONSTRAINT INDICATOR MUST BE 1-4.
!
!    -35    FOR J = ... THE LOWER BOUND = ... IS  >  THE UPPER BOUND =
!           ....
!
!    -36    THE INPUT ORDER OF COLUMNS = ... IS NOT BETWEEN 1 AND NCOLS
!           = ....
!
!    -37    THE BOUND POLARITY FLAG IN COMPONENT J = ... MUST BE
!           POSITIVE. NOW = ....
!
!    -38    THE ROW SEPARATOR TO APPLY WEIGHTING (...) MUST LIE BETWEEN
!           0 AND MINPUT = .... WEIGHT = ... MUST BE POSITIVE.
!
!***SEE ALSO  DBOCLS, DBOLS
!***ROUTINES CALLED  D1MACH, DAXPY, DCOPY, DDOT, DMOUT, DNRM2, DROT,
!                    DROTG, DSWAP, DVOUT, IVOUT, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   821220  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!   920422  Fixed usage of MINPUT.  (WRB)
!   901009  Editorial changes, code now reads from top to bottom.  (RWC)
!***END PROLOGUE  DBOLSM
!
!     PURPOSE
!     -------
!     THIS IS THE MAIN SUBPROGRAM THAT SOLVES THE BOUNDED
!     LEAST SQUARES PROBLEM.  THE PROBLEM SOLVED HERE IS:
!
!     SOLVE E*X =  F  (LEAST SQUARES SENSE)
!     WITH BOUNDS ON SELECTED X VALUES.
!
!     TO CHANGE THIS SUBPROGRAM FROM SINGLE TO DOUBLE PRECISION BEGIN
!     EDITING AT THE CARD 'C++'.
!     CHANGE THE SUBPROGRAM NAME TO DBOLSM AND THE STRINGS
!     /SAXPY/ TO /DAXPY/, /SCOPY/ TO /DCOPY/,
!     /SDOT/ TO /DDOT/, /SNRM2/ TO /DNRM2/,
!     /SROT/ TO /DROT/, /SROTG/ TO /DROTG/, /R1MACH/ TO /D1MACH/,
!     /SVOUT/ TO /DVOUT/, /SMOUT/ TO /DMOUT/,
!     /SSWAP/ TO /DSWAP/, /E0/ TO /D0/,
!     /REAL            / TO /DOUBLE PRECISION/.
!++
!
  DOUBLE PRECISION W(MDW,*),BL(*),BU(*)
  DOUBLE PRECISION X(*),RW(*),WW(*),SCL(*)
  DOUBLE PRECISION ALPHA,BETA,BOU,COLABV,COLBLO
  DOUBLE PRECISION CL1,CL2,CL3,ONE,BIG
  DOUBLE PRECISION FAC,RNORM,SC,SS,T,TOLIND,WT
  DOUBLE PRECISION TWO,T1,T2,WBIG,WLARGE,WMAG,XNEW
  DOUBLE PRECISION ZERO,DDOT,DNRM2
  DOUBLE PRECISION D1MACH,TOLSZE
  INTEGER IBASIS(*),IBB(*),IND(*),IOPT(*)
  LOGICAL FOUND,CONSTR
  CHARACTER*8 XERN1, XERN2
  CHARACTER*16 XERN3, XERN4
!
  PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
!
  INEXT(IDUM) = MIN(IDUM+1,MROWS)
!***FIRST EXECUTABLE STATEMENT  DBOLSM
!
!     Verify that the problem dimensions are defined properly.
!
  if (MINPUT <= 0) THEN
      WRITE (XERN1, '(I8)') MINPUT
      call XERMSG ('SLATEC', 'DBOLSM', 'THE NUMBER OF ROWS = ' // &
         XERN1 // ' MUST BE POSITIVE.', 31, 1)
      MODE = -31
      return
  end if
!
  if (NCOLS <= 0) THEN
      WRITE (XERN1, '(I8)') NCOLS
      call XERMSG ('SLATEC', 'DBOLSM', 'THE NUMBER OF COLUMNS = ' // &
         XERN1 // ' MUST BE POSITIVE.', 32, 1)
      MODE = -32
      return
  end if
!
  if (MDW < MINPUT) THEN
      WRITE (XERN1, '(I8)') MDW
      WRITE (XERN2, '(I8)') MINPUT
      call XERMSG ('SLATEC', 'DBOLSM', &
         'THE ROW DIMENSION OF W(,) = ' // XERN1 // &
         ' MUST BE  >=  THE NUMBER OF ROWS = ' // XERN2, 33, 1)
      MODE = -33
      return
  end if
!
!     Verify that bound information is correct.
!
  DO 10 J = 1,NCOLS
      if (IND(J) < 1 .OR. IND(J) > 4) THEN
          WRITE (XERN1, '(I8)') J
          WRITE (XERN2, '(I8)') IND(J)
          call XERMSG ('SLATEC', 'DBOLSM', 'FOR J = ' // XERN1 // &
             ' THE CONSTRAINT INDICATOR MUST BE 1-4', 34, 1)
          MODE = -34
          return
      ENDIF
   10 CONTINUE
!
  DO 20 J = 1,NCOLS
      if (IND(J) == 3) THEN
          if (BU(J) < BL(J)) THEN
              WRITE (XERN1, '(I8)') J
              WRITE (XERN3, '(1PD15.6)') BL(J)
              WRITE (XERN4, '(1PD15.6)') BU(J)
              call XERMSG ('SLATEC', 'DBOLSM', 'FOR J = ' // XERN1 &
                 // ' THE LOWER BOUND = ' // XERN3 // &
                 ' IS  >  THE UPPER BOUND = ' // XERN4, 35, 1)
              MODE = -35
              return
          ENDIF
      ENDIF
   20 CONTINUE
!
!     Check that permutation and polarity arrays have been set.
!
  DO 30 J = 1,NCOLS
      if (IBASIS(J) < 1 .OR. IBASIS(J) > NCOLS) THEN
          WRITE (XERN1, '(I8)') IBASIS(J)
          WRITE (XERN2, '(I8)') NCOLS
          call XERMSG ('SLATEC', 'DBOLSM', &
             'THE INPUT ORDER OF COLUMNS = ' // XERN1 // &
             ' IS NOT BETWEEN 1 AND NCOLS = ' // XERN2, 36, 1)
          MODE = -36
          return
      ENDIF
!
      if (IBB(J) <= 0) THEN
          WRITE (XERN1, '(I8)') J
          WRITE (XERN2, '(I8)') IBB(J)
          call XERMSG ('SLATEC', 'DBOLSM', &
             'THE BOUND POLARITY FLAG IN COMPONENT J = ' // XERN1 // &
             ' MUST BE POSITIVE.$$NOW = ' // XERN2, 37, 1)
          MODE = -37
          return
      ENDIF
   30 CONTINUE
!
!     Process the option array.
!
  FAC = 0.75D0
  TOLIND = SQRT(D1MACH(4))
  TOLSZE = SQRT(D1MACH(4))
  ITMAX = 5*MAX(MINPUT,NCOLS)
  WT = ONE
  MVAL = 0
  IPRINT = 0
!
!     Changes to some parameters can occur through the option array,
!     IOPT(*).  Process this array looking carefully for input data
!     errors.
!
  LP = 0
  LDS = 0
!
!     Test for no more options.
!
  590 LP = LP + LDS
  IP = IOPT(LP+1)
  JP = ABS(IP)
  if (IP == 99) THEN
      go to 470
  ELSE if (JP == 99) THEN
      LDS = 1
  ELSE if (JP == 1) THEN
!
!         Move the IOPT(*) processing pointer.
!
      if (IP > 0) THEN
          LP = IOPT(LP+2) - 1
          LDS = 0
      ELSE
          LDS = 2
      ENDIF
  ELSE if (JP == 2) THEN
!
!         Change tolerance for rank determination.
!
      if (IP > 0) THEN
          IOFF = IOPT(LP+2)
          if (IOFF <= 0) THEN
              WRITE (XERN1, '(I8)') IOFF
              WRITE (XERN2, '(I8)') NCOLS
              call XERMSG ('SLATEC', 'DBOLSM', 'THE OFFSET = ' // &
                 XERN1 // ' BEYOND POSITION NCOLS = ' // XERN2 // &
                 ' MUST BE POSITIVE FOR OPTION NUMBER 2.', 24, 1)
              MODE = -24
              return
          ENDIF
!
          TOLIND = X(NCOLS+IOFF)
          if (TOLIND < D1MACH(4)) THEN
              WRITE (XERN3, '(1PD15.6)') TOLIND
              WRITE (XERN4, '(1PD15.6)') D1MACH(4)
              call XERMSG ('SLATEC', 'DBOLSM', &
                 'THE TOLERANCE FOR RANK DETERMINATION = ' // XERN3 &
                 // ' IS LESS THAN MACHINE PRECISION = ' // XERN4, &
                 25, 0)
              MODE = -25
          ENDIF
      ENDIF
      LDS = 2
  ELSE if (JP == 3) THEN
!
!         Change blowup factor for allowing variables to become
!         inactive.
!
      if (IP > 0) THEN
          IOFF = IOPT(LP+2)
          if (IOFF <= 0) THEN
              WRITE (XERN1, '(I8)') IOFF
              WRITE (XERN2, '(I8)') NCOLS
              call XERMSG ('SLATEC', 'DBOLSM', 'THE OFFSET = ' // &
                 XERN1 // ' BEYOND POSITION NCOLS = ' // XERN2 // &
                 ' MUST BE POSITIVE FOR OPTION NUMBER 3.', 26, 1)
              MODE = -26
              return
          ENDIF
!
          TOLSZE = X(NCOLS+IOFF)
          if (TOLSZE <= ZERO) THEN
              WRITE (XERN3, '(1PD15.6)') TOLSZE
              call XERMSG ('SLATEC', 'DBOLSM', 'THE RECIPROCAL ' // &
                 'OF THE BLOW-UP FACTOR FOR REJECTING VARIABLES ' // &
                 'MUST BE POSITIVE.$$NOW = ' // XERN3, 27, 1)
              MODE = -27
              return
          ENDIF
      ENDIF
      LDS = 2
  ELSE if (JP == 4) THEN
!
!         Change the maximum number of iterations allowed.
!
      if (IP > 0) THEN
          ITMAX = IOPT(LP+2)
          if (ITMAX <= 0) THEN
              WRITE (XERN1, '(I8)') ITMAX
              call XERMSG ('SLATEC', 'DBOLSM', &
                 'THE MAXIMUM NUMBER OF ITERATIONS = ' // XERN1 // &
                 ' MUST BE POSITIVE.', 28, 1)
              MODE = -28
              return
          ENDIF
      ENDIF
      LDS = 2
  ELSE if (JP == 5) THEN
!
!         Change the factor for pretriangularizing the data matrix.
!
      if (IP > 0) THEN
          IOFF = IOPT(LP+2)
          if (IOFF <= 0) THEN
              WRITE (XERN1, '(I8)') IOFF
              WRITE (XERN2, '(I8)') NCOLS
              call XERMSG ('SLATEC', 'DBOLSM', 'THE OFFSET = ' // &
                 XERN1 // ' BEYOND POSITION NCOLS = ' // XERN2 // &
                 ' MUST BE POSITIVE FOR OPTION NUMBER 5.', 29, 1)
              MODE = -29
              return
          ENDIF
!
          FAC = X(NCOLS+IOFF)
          if (FAC < ZERO) THEN
              WRITE (XERN3, '(1PD15.6)') FAC
              call XERMSG ('SLATEC', 'DBOLSM', &
                 'THE FACTOR (NCOLS/MINPUT) WHERE PRE-' // &
                 'TRIANGULARIZING IS PERFORMED MUST BE NON-' // &
                 'NEGATIVE.$$NOW = ' // XERN3, 30, 0)
              MODE = -30
              return
          ENDIF
      ENDIF
      LDS = 2
  ELSE if (JP == 6) THEN
!
!         Change the weighting factor (from 1.0) to apply to components
!         numbered .gt. MVAL (initially set to 1.)  This trick is needed
!         for applications of this subprogram to the heavily weighted
!         least squares problem that come from equality constraints.
!
      if (IP > 0) THEN
          IOFF = IOPT(LP+2)
          MVAL = IOPT(LP+3)
          WT = X(NCOLS+IOFF)
      ENDIF
!
      if (MVAL < 0 .OR. MVAL > MINPUT .OR. WT <= ZERO) THEN
          WRITE (XERN1, '(I8)') MVAL
          WRITE (XERN2, '(I8)') MINPUT
          WRITE (XERN3, '(1PD15.6)') WT
          call XERMSG ('SLATEC', 'DBOLSM', &
             'THE ROW SEPARATOR TO APPLY WEIGHTING (' // XERN1 // &
             ') MUST LIE BETWEEN 0 AND MINPUT = ' // XERN2 // &
             '.$$WEIGHT = ' // XERN3 // ' MUST BE POSITIVE.', 38, 0)
          MODE = -38
          return
      ENDIF
      LDS = 3
  ELSE if (JP == 7) THEN
!
!         Turn on debug output.
!
      if (IP > 0) IPRINT = 1
      LDS = 2
  ELSE
      WRITE (XERN1, '(I8)') IP
      call XERMSG ('SLATEC', 'DBOLSM', 'THE OPTION NUMBER = ' // &
         XERN1 // ' IS NOT DEFINED.', 23, 1)
      MODE = -23
      return
  end if
  go to 590
!
!     Pretriangularize rectangular arrays of certain sizes for
!     increased efficiency.
!
  470 if (FAC*MINPUT > NCOLS) THEN
      DO 490 J = 1,NCOLS+1
          DO 480 I = MINPUT,J+MVAL+1,-1
              call DROTG(W(I-1,J),W(I,J),SC,SS)
              W(I,J) = ZERO
              call DROT(NCOLS-J+1,W(I-1,J+1),MDW,W(I,J+1),MDW,SC,SS)
  480         CONTINUE
  490     CONTINUE
      MROWS = NCOLS + MVAL + 1
  ELSE
      MROWS = MINPUT
  end if
!
!     Set the X(*) array to zero so all components are defined.
!
  call DCOPY(NCOLS,ZERO,0,X,1)
!
!     The arrays IBASIS(*) and IBB(*) are initialized by the calling
!     program and the column scaling is defined in the calling program.
!     'BIG' is plus infinity on this machine.
!
  BIG = D1MACH(2)
  DO 550 J = 1,NCOLS
      if (IND(J) == 1) THEN
          BU(J) = BIG
      ELSE if (IND(J) == 2) THEN
          BL(J) = -BIG
      ELSE if (IND(J) == 4) THEN
          BL(J) = -BIG
          BU(J) = BIG
      ENDIF
  550 CONTINUE
!
  DO 570 J = 1,NCOLS
      if ((BL(J) <= ZERO.AND.ZERO <= BU(J).AND.ABS(BU(J)) <  &
          ABS(BL(J))) .OR. BU(J) < ZERO) THEN
          T = BU(J)
          BU(J) = -BL(J)
          BL(J) = -T
          SCL(J) = -SCL(J)
          DO 560 I = 1,MROWS
              W(I,J) = -W(I,J)
  560         CONTINUE
      ENDIF
!
!         Indices in set T(=TIGHT) are denoted by negative values
!         of IBASIS(*).
!
      if (BL(J) >= ZERO) THEN
          IBASIS(J) = -IBASIS(J)
          T = -BL(J)
          BU(J) = BU(J) + T
          call DAXPY(MROWS,T,W(1,J),1,W(1,NCOLS+1),1)
      ENDIF
  570 CONTINUE
!
  NSETB = 0
  ITER = 0
!
  if (IPRINT > 0) THEN
      call DMOUT(MROWS,NCOLS+1,MDW,W,'('' PRETRI. INPUT MATRIX'')', &
                 -4)
      call DVOUT(NCOLS,BL,'('' LOWER BOUNDS'')',-4)
      call DVOUT(NCOLS,BU,'('' UPPER BOUNDS'')',-4)
  end if
!
  580 ITER = ITER + 1
  if (ITER > ITMAX) THEN
     WRITE (XERN1, '(I8)') ITMAX
     call XERMSG ('SLATEC', 'DBOLSM', 'MORE THAN ITMAX = ' // XERN1 &
        // ' ITERATIONS SOLVING BOUNDED LEAST SQUARES PROBLEM.', &
        22, 1)
     MODE = -22
!
!        Rescale and translate variables.
!
     IGOPR = 1
     go to 130
  end if
!
!     Find a variable to become non-active.
!                                                 T
!     Compute (negative) of gradient vector, W = E *(F-E*X).
!
  call DCOPY(NCOLS,ZERO,0,WW,1)
  DO 200 J = NSETB+1,NCOLS
      JCOL = ABS(IBASIS(J))
      WW(J) = DDOT(MROWS-NSETB,W(INEXT(NSETB),J),1, &
              W(INEXT(NSETB),NCOLS+1),1)*ABS(SCL(JCOL))
  200 CONTINUE
!
  if (IPRINT > 0) THEN
      call DVOUT(NCOLS,WW,'('' GRADIENT VALUES'')',-4)
      call IVOUT(NCOLS,IBASIS,'('' INTERNAL VARIABLE ORDER'')',-4)
      call IVOUT(NCOLS,IBB,'('' BOUND POLARITY'')',-4)
  end if
!
!     If active set = number of total rows, quit.
!
  210 if (NSETB == MROWS) THEN
      FOUND = .FALSE.
      go to 120
  end if
!
!     Choose an extremal component of gradient vector for a candidate
!     to become non-active.
!
  WLARGE = -BIG
  WMAG = -BIG
  DO 220 J = NSETB+1,NCOLS
      T = WW(J)
      if (T == BIG) go to 220
      ITEMP = IBASIS(J)
      JCOL = ABS(ITEMP)
      T1 = DNRM2(MVAL-NSETB,W(INEXT(NSETB),J),1)
      if (ITEMP < 0) THEN
          if (MOD(IBB(JCOL),2) == 0) T = -T
          if (T < ZERO) go to 220
          if (MVAL > NSETB) T = T1
          if (T > WLARGE) THEN
              WLARGE = T
              JLARGE = J
          ENDIF
      ELSE
          if (MVAL > NSETB) T = T1
          if (ABS(T) > WMAG) THEN
              WMAG = ABS(T)
              JMAG = J
          ENDIF
      ENDIF
  220 CONTINUE
!
!     Choose magnitude of largest component of gradient for candidate.
!
  JBIG = 0
  WBIG = ZERO
  if (WLARGE > ZERO) THEN
      JBIG = JLARGE
      WBIG = WLARGE
  end if
!
  if (WMAG >= WBIG) THEN
      JBIG = JMAG
      WBIG = WMAG
  end if
!
  if (JBIG == 0) THEN
      FOUND = .FALSE.
      if (IPRINT > 0) THEN
          call IVOUT(0,I,'('' FOUND NO VARIABLE TO ENTER'')',-4)
      ENDIF
      go to 120
  end if
!
!     See if the incoming column is sufficiently independent.  This
!     test is made before an elimination is performed.
!
  if (IPRINT > 0) &
      call IVOUT(1,JBIG,'('' TRY TO BRING IN THIS COL.'')',-4)
!
  if (MVAL <= NSETB) THEN
      CL1 = DNRM2(MVAL,W(1,JBIG),1)
      CL2 = ABS(WT)*DNRM2(NSETB-MVAL,W(INEXT(MVAL),JBIG),1)
      CL3 = ABS(WT)*DNRM2(MROWS-NSETB,W(INEXT(NSETB),JBIG),1)
      call DROTG(CL1,CL2,SC,SS)
      COLABV = ABS(CL1)
      COLBLO = CL3
  ELSE
      CL1 = DNRM2(NSETB,W(1,JBIG),1)
      CL2 = DNRM2(MVAL-NSETB,W(INEXT(NSETB),JBIG),1)
      CL3 = ABS(WT)*DNRM2(MROWS-MVAL,W(INEXT(MVAL),JBIG),1)
      COLABV = CL1
      call DROTG(CL2,CL3,SC,SS)
      COLBLO = ABS(CL2)
  end if
!
  if (COLBLO <= TOLIND*COLABV) THEN
      WW(JBIG) = BIG
      if (IPRINT > 0) &
          call IVOUT(0,I,'('' VARIABLE IS DEPENDENT, NOT USED.'')', &
             -4)
      go to 210
  end if
!
!     Swap matrix columns NSETB+1 and JBIG, plus pointer information,
!     and gradient values.
!
  NSETB = NSETB + 1
  if (NSETB /= JBIG) THEN
      call DSWAP(MROWS,W(1,NSETB),1,W(1,JBIG),1)
      call DSWAP(1,WW(NSETB),1,WW(JBIG),1)
      ITEMP = IBASIS(NSETB)
      IBASIS(NSETB) = IBASIS(JBIG)
      IBASIS(JBIG) = ITEMP
  end if
!
!     Eliminate entries below the pivot line in column NSETB.
!
  if (MROWS > NSETB) THEN
      DO 230 I = MROWS,NSETB+1,-1
          if (I == MVAL+1) go to 230
          call DROTG(W(I-1,NSETB),W(I,NSETB),SC,SS)
          W(I,NSETB) = ZERO
          call DROT(NCOLS-NSETB+1,W(I-1,NSETB+1),MDW,W(I,NSETB+1), &
                    MDW,SC,SS)
  230     CONTINUE
!
      if (MVAL >= NSETB .AND. MVAL < MROWS) THEN
          call DROTG(W(NSETB,NSETB),W(MVAL+1,NSETB),SC,SS)
          W(MVAL+1,NSETB) = ZERO
          call DROT(NCOLS-NSETB+1,W(NSETB,NSETB+1),MDW, &
                    W(MVAL+1,NSETB+1),MDW,SC,SS)
      ENDIF
  end if
!
  if (W(NSETB,NSETB) == ZERO) THEN
      WW(NSETB) = BIG
      NSETB = NSETB - 1
      if (IPRINT > 0) THEN
          call IVOUT(0,I,'('' PIVOT IS ZERO, NOT USED.'')',-4)
      ENDIF
      go to 210
  end if
!
!     Check that new variable is moving in the right direction.
!
  ITEMP = IBASIS(NSETB)
  JCOL = ABS(ITEMP)
  XNEW = (W(NSETB,NCOLS+1)/W(NSETB,NSETB))/ABS(SCL(JCOL))
  if (ITEMP < 0) THEN
!
!         if ( WW(NSETB) >= ZERO.AND.XNEW <= ZERO) exit(quit)
!         if ( WW(NSETB) <= ZERO.AND.XNEW >= ZERO) exit(quit)
!
      if ((WW(NSETB) >= ZERO.AND.XNEW <= ZERO) .OR. &
          (WW(NSETB) <= ZERO.AND.XNEW >= ZERO)) go to 240
  end if
  FOUND = .TRUE.
  go to 120
!
  240 WW(NSETB) = BIG
  NSETB = NSETB - 1
  if (IPRINT > 0) &
      call IVOUT(0,I,'('' VARIABLE HAS BAD DIRECTION, NOT USED.'')', &
         -4)
  go to 210
!
!     Solve the triangular system.
!
  270 call DCOPY(NSETB,W(1,NCOLS+1),1,RW,1)
  DO 280 J = NSETB,1,-1
      RW(J) = RW(J)/W(J,J)
      JCOL = ABS(IBASIS(J))
      T = RW(J)
      if (MOD(IBB(JCOL),2) == 0) RW(J) = -RW(J)
      call DAXPY(J-1,-T,W(1,J),1,RW,1)
      RW(J) = RW(J)/ABS(SCL(JCOL))
  280 CONTINUE
!
  if (IPRINT > 0) THEN
      call DVOUT(NSETB,RW,'('' SOLN. VALUES'')',-4)
      call IVOUT(NSETB,IBASIS,'('' COLS. USED'')',-4)
  end if
!
  if (LGOPR == 2) THEN
      call DCOPY(NSETB,RW,1,X,1)
      DO 450 J = 1,NSETB
          ITEMP = IBASIS(J)
          JCOL = ABS(ITEMP)
          if (ITEMP < 0) THEN
              BOU = ZERO
          ELSE
              BOU = BL(JCOL)
          ENDIF
!
          if ((-BOU) /= BIG) BOU = BOU/ABS(SCL(JCOL))
          if (X(J) <= BOU) THEN
              JDROP1 = J
              go to 340
          ENDIF
!
          BOU = BU(JCOL)
          if (BOU /= BIG) BOU = BOU/ABS(SCL(JCOL))
          if (X(J) >= BOU) THEN
              JDROP2 = J
              go to 340
          ENDIF
  450     CONTINUE
      go to 340
  end if
!
!     See if the unconstrained solution (obtained by solving the
!     triangular system) satisfies the problem bounds.
!
  ALPHA = TWO
  BETA = TWO
  X(NSETB) = ZERO
  DO 310 J = 1,NSETB
      ITEMP = IBASIS(J)
      JCOL = ABS(ITEMP)
      T1 = TWO
      T2 = TWO
      if (ITEMP < 0) THEN
          BOU = ZERO
      ELSE
          BOU = BL(JCOL)
      ENDIF
      if ((-BOU) /= BIG) BOU = BOU/ABS(SCL(JCOL))
      if (RW(J) <= BOU) T1 = (X(J)-BOU)/ (X(J)-RW(J))
      BOU = BU(JCOL)
      if (BOU /= BIG) BOU = BOU/ABS(SCL(JCOL))
      if (RW(J) >= BOU) T2 = (BOU-X(J))/ (RW(J)-X(J))
!
!     If not, then compute a step length so that the variables remain
!     feasible.
!
      if (T1 < ALPHA) THEN
          ALPHA = T1
          JDROP1 = J
      ENDIF
!
      if (T2 < BETA) THEN
          BETA = T2
          JDROP2 = J
      ENDIF
  310 CONTINUE
!
  CONSTR = ALPHA  <  TWO .OR. BETA  <  TWO
  if (.NOT.CONSTR) THEN
!
!         Accept the candidate because it satisfies the stated bounds
!         on the variables.
!
      call DCOPY(NSETB,RW,1,X,1)
      go to 580
  end if
!
!     Take a step that is as large as possible with all variables
!     remaining feasible.
!
  DO 330 J = 1,NSETB
      X(J) = X(J) + MIN(ALPHA,BETA)* (RW(J)-X(J))
  330 CONTINUE
!
  if (ALPHA <= BETA) THEN
      JDROP2 = 0
  ELSE
      JDROP1 = 0
  end if
!
  340 if (JDROP1+JDROP2 <= 0 .OR. NSETB <= 0) go to 580
  350 JDROP = JDROP1 + JDROP2
  ITEMP = IBASIS(JDROP)
  JCOL = ABS(ITEMP)
  if (JDROP2 > 0) THEN
!
!         Variable is at an upper bound.  Subtract multiple of this
!         column from right hand side.
!
      T = BU(JCOL)
      if (ITEMP > 0) THEN
          BU(JCOL) = T - BL(JCOL)
          BL(JCOL) = -T
          ITEMP = -ITEMP
          SCL(JCOL) = -SCL(JCOL)
          DO 360 I = 1,JDROP
              W(I,JDROP) = -W(I,JDROP)
  360         CONTINUE
      ELSE
          IBB(JCOL) = IBB(JCOL) + 1
          if (MOD(IBB(JCOL),2) == 0) T = -T
      ENDIF
!
!     Variable is at a lower bound.
!
  ELSE
      if (ITEMP < ZERO) THEN
          T = ZERO
      ELSE
          T = -BL(JCOL)
          BU(JCOL) = BU(JCOL) + T
          ITEMP = -ITEMP
      ENDIF
  end if
!
  call DAXPY(JDROP,T,W(1,JDROP),1,W(1,NCOLS+1),1)
!
!     Move certain columns left to achieve upper Hessenberg form.
!
  call DCOPY(JDROP,W(1,JDROP),1,RW,1)
  DO 370 J = JDROP+1,NSETB
      IBASIS(J-1) = IBASIS(J)
      X(J-1) = X(J)
      call DCOPY(J,W(1,J),1,W(1,J-1),1)
  370 CONTINUE
!
  IBASIS(NSETB) = ITEMP
  W(1,NSETB) = ZERO
  call DCOPY(MROWS-JDROP,W(1,NSETB),0,W(JDROP+1,NSETB),1)
  call DCOPY(JDROP,RW,1,W(1,NSETB),1)
!
!     Transform the matrix from upper Hessenberg form to upper
!     triangular form.
!
  NSETB = NSETB - 1
  DO 390 I = JDROP,NSETB
!
!         Look for small pivots and avoid mixing weighted and
!         nonweighted rows.
!
      if (I == MVAL) THEN
          T = ZERO
          DO 380 J = I,NSETB
              JCOL = ABS(IBASIS(J))
              T1 = ABS(W(I,J)*SCL(JCOL))
              if (T1 > T) THEN
                  JBIG = J
                  T = T1
              ENDIF
  380         CONTINUE
          go to 400
      ENDIF
      call DROTG(W(I,I),W(I+1,I),SC,SS)
      W(I+1,I) = ZERO
      call DROT(NCOLS-I+1,W(I,I+1),MDW,W(I+1,I+1),MDW,SC,SS)
  390 CONTINUE
  go to 430
!
!     The triangularization is completed by giving up the Hessenberg
!     form and triangularizing a rectangular matrix.
!
  400 call DSWAP(MROWS,W(1,I),1,W(1,JBIG),1)
  call DSWAP(1,WW(I),1,WW(JBIG),1)
  call DSWAP(1,X(I),1,X(JBIG),1)
  ITEMP = IBASIS(I)
  IBASIS(I) = IBASIS(JBIG)
  IBASIS(JBIG) = ITEMP
  JBIG = I
  DO 420 J = JBIG,NSETB
      DO 410 I = J+1,MROWS
          call DROTG(W(J,J),W(I,J),SC,SS)
          W(I,J) = ZERO
          call DROT(NCOLS-J+1,W(J,J+1),MDW,W(I,J+1),MDW,SC,SS)
  410     CONTINUE
  420 CONTINUE
!
!     See if the remaining coefficients are feasible.  They should be
!     because of the way MIN(ALPHA,BETA) was chosen.  Any that are not
!     feasible will be set to their bounds and appropriately translated.
!
  430 JDROP1 = 0
  JDROP2 = 0
  LGOPR = 2
  go to 270
!
!     Find a variable to become non-active.
!
  120 if (FOUND) THEN
      LGOPR = 1
      go to 270
  end if
!
!     Rescale and translate variables.
!
  IGOPR = 2
  130 call DCOPY(NSETB,X,1,RW,1)
  call DCOPY(NCOLS,ZERO,0,X,1)
  DO 140 J = 1,NSETB
      JCOL = ABS(IBASIS(J))
      X(JCOL) = RW(J)*ABS(SCL(JCOL))
  140 CONTINUE
!
  DO 150 J = 1,NCOLS
      if (MOD(IBB(J),2) == 0) X(J) = BU(J) - X(J)
  150 CONTINUE
!
  DO 160 J = 1,NCOLS
      JCOL = IBASIS(J)
      if (JCOL < 0) X(-JCOL) = BL(-JCOL) + X(-JCOL)
  160 CONTINUE
!
  DO 170 J = 1,NCOLS
      if (SCL(J) < ZERO) X(J) = -X(J)
  170 CONTINUE
!
  I = MAX(NSETB,MVAL)
  RNORM = DNRM2(MROWS-I,W(INEXT(I),NCOLS+1),1)
!
  if (IGOPR == 2) MODE = NSETB
  return
end
