subroutine DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL, &
     IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
!
!! DDASSL solves a system of differential/algebraic equations...
!  of the form G(T,Y,YPRIME) = 0.
!
!***LIBRARY   SLATEC (DASSL)
!***CATEGORY  I1A2
!***TYPE      DOUBLE PRECISION (SDASSL-S, DDASSL-D)
!***KEYWORDS  BACKWARD DIFFERENTIATION FORMULAS, DASSL,
!             DIFFERENTIAL/ALGEBRAIC, IMPLICIT DIFFERENTIAL SYSTEMS
!***AUTHOR  Petzold, Linda R., (LLNL)
!             Computing and Mathematics Research Division
!             Lawrence Livermore National Laboratory
!             L - 316, P.O. Box 808,
!             Livermore, CA.    94550
!***DESCRIPTION
!
! *Usage:
!
!      EXTERNAL RES, JAC
!      INTEGER NEQ, INFO(N), IDID, LRW, LIW, IWORK(LIW), IPAR
!      DOUBLE PRECISION T, Y(NEQ), YPRIME(NEQ), TOUT, RTOL, ATOL,
!     *   RWORK(LRW), RPAR
!
!      call DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
!     *   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
!
!
! *Arguments:
!  (In the following, all real arrays should be type DOUBLE PRECISION.)
!
!  RES:EXT     This is a subroutine which you provide to define the
!              differential/algebraic system.
!
!  NEQ:IN      This is the number of equations to be solved.
!
!  T:INOUT     This is the current value of the independent variable.
!
!  Y(*):INOUT  This array contains the solution components at T.
!
!  YPRIME(*):INOUT  This array contains the derivatives of the solution
!              components at T.
!
!  TOUT:IN     This is a point at which a solution is desired.
!
!  INFO(N):IN  The basic task of the code is to solve the system from T
!              to TOUT and return an answer at TOUT.  INFO is an integer
!              array which is used to communicate exactly how you want
!              this task to be carried out.  (See below for details.)
!              N must be greater than or equal to 15.
!
!  RTOL,ATOL:INOUT  These quantities represent relative and absolute
!              error tolerances which you provide to indicate how
!              accurately you wish the solution to be computed.  You
!              may choose them to be both scalars or else both vectors.
!              Caution:  In Fortran 77, a scalar is not the same as an
!                        array of length 1.  Some compilers may object
!                        to using scalars for RTOL,ATOL.
!
!  IDID:OUT    This scalar quantity is an indicator reporting what the
!              code did.  You must monitor this integer variable to
!              decide  what action to take next.
!
!  RWORK:WORK  A real work array of length LRW which provides the
!              code with needed storage space.
!
!  LRW:IN      The length of RWORK.  (See below for required length.)
!
!  IWORK:WORK  An integer work array of length LIW which provides the
!              code with needed storage space.
!
!  LIW:IN      The length of IWORK.  (See below for required length.)
!
!  RPAR,IPAR:IN  These are real and integer parameter arrays which
!              you can use for communication between your calling
!              program and the RES subroutine (and the JAC subroutine)
!
!  JAC:EXT     This is the name of a subroutine which you may choose
!              to provide for defining a matrix of partial derivatives
!              described below.
!
!  Quantities which may be altered by DDASSL are:
!     T, Y(*), YPRIME(*), INFO(1), RTOL, ATOL,
!     IDID, RWORK(*) AND IWORK(*)
!
! *Description
!
!  Subroutine DDASSL uses the backward differentiation formulas of
!  orders one through five to solve a system of the above form for Y and
!  YPRIME.  Values for Y and YPRIME at the initial time must be given as
!  input.  These values must be consistent, (that is, if T,Y,YPRIME are
!  the given initial values, they must satisfy G(T,Y,YPRIME) = 0.).  The
!  subroutine solves the system from T to TOUT.  It is easy to continue
!  the solution to get results at additional TOUT.  This is the interval
!  mode of operation.  Intermediate results can also be obtained easily
!  by using the intermediate-output capability.
!
!  The following detailed description is divided into subsections:
!    1. Input required for the first call to DDASSL.
!    2. Output after any return from DDASSL.
!    3. What to do to continue the integration.
!    4. Error messages.
!
!
!  -------- INPUT -- WHAT TO DO ON THE FIRST call TO DDASSL ------------
!
!  The first call of the code is defined to be the start of each new
!  problem. Read through the descriptions of all the following items,
!  provide sufficient storage space for designated arrays, set
!  appropriate variables for the initialization of the problem, and
!  give information about how you want the problem to be solved.
!
!
!  RES -- Provide a subroutine of the form
!             SUBROUTINE RES(T,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
!         to define the system of differential/algebraic
!         equations which is to be solved. For the given values
!         of T,Y and YPRIME, the subroutine should
!         return the residual of the differential/algebraic
!         system
!             DELTA = G(T,Y,YPRIME)
!         (DELTA(*) is a vector of length NEQ which is
!         output for RES.)
!
!         Subroutine RES must not alter T,Y or YPRIME.
!         You must declare the name RES in an external
!         statement in your program that calls DDASSL.
!         You must dimension Y,YPRIME and DELTA in RES.
!
!         IRES is an integer flag which is always equal to
!         zero on input. Subroutine RES should alter IRES
!         only if it encounters an illegal value of Y or
!         a stop condition. Set IRES = -1 if an input value
!         is illegal, and DDASSL will try to solve the problem
!         without getting IRES = -1. If IRES = -2, DDASSL
!         will return control to the calling program
!         with IDID = -11.
!
!         RPAR and IPAR are real and integer parameter arrays which
!         you can use for communication between your calling program
!         and subroutine RES. They are not altered by DDASSL. If you
!         do not need RPAR or IPAR, ignore these parameters by treat-
!         ing them as dummy arguments. If you do choose to use them,
!         dimension them in your calling program and in RES as arrays
!         of appropriate length.
!
!  NEQ -- Set it to the number of differential equations.
!         (NEQ  >=  1)
!
!  T -- Set it to the initial point of the integration.
!         T must be defined as a variable.
!
!  Y(*) -- Set this vector to the initial values of the NEQ solution
!         components at the initial point. You must dimension Y of
!         length at least NEQ in your calling program.
!
!  YPRIME(*) -- Set this vector to the initial values of the NEQ
!         first derivatives of the solution components at the initial
!         point.  You must dimension YPRIME at least NEQ in your
!         calling program. If you do not know initial values of some
!         of the solution components, see the explanation of INFO(11).
!
!  TOUT -- Set it to the first point at which a solution
!         is desired. You can not take TOUT = T.
!         integration either forward in T (TOUT  >  T) or
!         backward in T (TOUT  <  T) is permitted.
!
!         The code advances the solution from T to TOUT using
!         step sizes which are automatically selected so as to
!         achieve the desired accuracy. If you wish, the code will
!         return with the solution and its derivative at
!         intermediate steps (intermediate-output mode) so that
!         you can monitor them, but you still must provide TOUT in
!         accord with the basic aim of the code.
!
!         The first step taken by the code is a critical one
!         because it must reflect how fast the solution changes near
!         the initial point. The code automatically selects an
!         initial step size which is practically always suitable for
!         the problem. By using the fact that the code will not step
!         past TOUT in the first step, you could, if necessary,
!         restrict the length of the initial step size.
!
!         For some problems it may not be permissible to integrate
!         past a point TSTOP because a discontinuity occurs there
!         or the solution or its derivative is not defined beyond
!         TSTOP. When you have declared a TSTOP point (SEE INFO(4)
!         and RWORK(1)), you have told the code not to integrate
!         past TSTOP. In this case any TOUT beyond TSTOP is invalid
!         input.
!
!  INFO(*) -- Use the INFO array to give the code more details about
!         how you want your problem solved.  This array should be
!         dimensioned of length 15, though DDASSL uses only the first
!         eleven entries.  You must respond to all of the following
!         items, which are arranged as questions.  The simplest use
!         of the code corresponds to answering all questions as yes,
!         i.e. setting all entries of INFO to 0.
!
!       INFO(1) - This parameter enables the code to initialize
!              itself. You must set it to indicate the start of every
!              new problem.
!
!          **** Is this the first call for this problem ...
!                Yes - Set INFO(1) = 0
!                 No - Not applicable here.
!                      See below for continuation calls.  ****
!
!       INFO(2) - How much accuracy you want of your solution
!              is specified by the error tolerances RTOL and ATOL.
!              The simplest use is to take them both to be scalars.
!              To obtain more flexibility, they can both be vectors.
!              The code must be told your choice.
!
!          **** Are both error tolerances RTOL, ATOL scalars ...
!                Yes - Set INFO(2) = 0
!                      and input scalars for both RTOL and ATOL
!                 No - Set INFO(2) = 1
!                      and input arrays for both RTOL and ATOL ****
!
!       INFO(3) - The code integrates from T in the direction
!              of TOUT by steps. If you wish, it will return the
!              computed solution and derivative at the next
!              intermediate step (the intermediate-output mode) or
!              TOUT, whichever comes first. This is a good way to
!              proceed if you want to see the behavior of the solution.
!              If you must have solutions at a great many specific
!              TOUT points, this code will compute them efficiently.
!
!          **** Do you want the solution only at
!                TOUT (and not at the next intermediate step) ...
!                 Yes - Set INFO(3) = 0
!                  No - Set INFO(3) = 1 ****
!
!       INFO(4) - To handle solutions at a great many specific
!              values TOUT efficiently, this code may integrate past
!              TOUT and interpolate to obtain the result at TOUT.
!              Sometimes it is not possible to integrate beyond some
!              point TSTOP because the equation changes there or it is
!              not defined past TSTOP. Then you must tell the code
!              not to go past.
!
!           **** Can the integration be carried out without any
!                restrictions on the independent variable T ...
!                 Yes - Set INFO(4)=0
!                  No - Set INFO(4)=1
!                       and define the stopping point TSTOP by
!                       setting RWORK(1)=TSTOP ****
!
!       INFO(5) - To solve differential/algebraic problems it is
!              necessary to use a matrix of partial derivatives of the
!              system of differential equations. If you do not
!              provide a subroutine to evaluate it analytically (see
!              description of the item JAC in the call list), it will
!              be approximated by numerical differencing in this code.
!              although it is less trouble for you to have the code
!              compute partial derivatives by numerical differencing,
!              the solution will be more reliable if you provide the
!              derivatives via JAC. Sometimes numerical differencing
!              is cheaper than evaluating derivatives in JAC and
!              sometimes it is not - this depends on your problem.
!
!           **** Do you want the code to evaluate the partial
!                derivatives automatically by numerical differences ...
!                   Yes - Set INFO(5)=0
!                    No - Set INFO(5)=1
!                  and provide subroutine JAC for evaluating the
!                  matrix of partial derivatives ****
!
!       INFO(6) - DDASSL will perform much better if the matrix of
!              partial derivatives, DG/DY + CJ*DG/DYPRIME,
!              (here CJ is a scalar determined by DDASSL)
!              is banded and the code is told this. In this
!              case, the storage needed will be greatly reduced,
!              numerical differencing will be performed much cheaper,
!              and a number of important algorithms will execute much
!              faster. The differential equation is said to have
!              half-bandwidths ML (lower) and MU (upper) if equation i
!              involves only unknowns Y(J) with
!                             I-ML  <=  J  <=  I+MU
!              for all I=1,2,...,NEQ. Thus, ML and MU are the widths
!              of the lower and upper parts of the band, respectively,
!              with the main diagonal being excluded. If you do not
!              indicate that the equation has a banded matrix of partial
!              derivatives, the code works with a full matrix of NEQ**2
!              elements (stored in the conventional way). Computations
!              with banded matrices cost less time and storage than with
!              full matrices if 2*ML+MU  <  NEQ. If you tell the
!              code that the matrix of partial derivatives has a banded
!              structure and you want to provide subroutine JAC to
!              compute the partial derivatives, then you must be careful
!              to store the elements of the matrix in the special form
!              indicated in the description of JAC.
!
!          **** Do you want to solve the problem using a full
!               (dense) matrix (and not a special banded
!               structure) ...
!                Yes - Set INFO(6)=0
!                 No - Set INFO(6)=1
!                       and provide the lower (ML) and upper (MU)
!                       bandwidths by setting
!                       IWORK(1)=ML
!                       IWORK(2)=MU ****
!
!
!        INFO(7) -- You can specify a maximum (absolute value of)
!              stepsize, so that the code
!              will avoid passing over very
!              large regions.
!
!          ****  Do you want the code to decide
!                on its own maximum stepsize?
!                Yes - Set INFO(7)=0
!                 No - Set INFO(7)=1
!                      and define HMAX by setting
!                      RWORK(2)=HMAX ****
!
!        INFO(8) -- Differential/algebraic problems
!              may occasionally suffer from
!              severe scaling difficulties on the
!              first step. If you know a great deal
!              about the scaling of your problem, you can
!              help to alleviate this problem by
!              specifying an initial stepsize HO.
!
!          ****  Do you want the code to define
!                its own initial stepsize?
!                Yes - Set INFO(8)=0
!                 No - Set INFO(8)=1
!                      and define HO by setting
!                      RWORK(3)=HO ****
!
!        INFO(9) -- If storage is a severe problem,
!              you can save some locations by
!              restricting the maximum order MAXORD.
!              the default value is 5. for each
!              order decrease below 5, the code
!              requires NEQ fewer locations, however
!              it is likely to be slower. In any
!              case, you must have 1  <=  MAXORD  <=  5
!          ****  Do you want the maximum order to
!                default to 5?
!                Yes - Set INFO(9)=0
!                 No - Set INFO(9)=1
!                      and define MAXORD by setting
!                      IWORK(3)=MAXORD ****
!
!        INFO(10) --If you know that the solutions to your equations
!               will always be nonnegative, it may help to set this
!               parameter. However, it is probably best to
!               try the code without using this option first,
!               and only to use this option if that doesn't
!               work very well.
!           ****  Do you want the code to solve the problem without
!                 invoking any special nonnegativity constraints?
!                  Yes - Set INFO(10)=0
!                   No - Set INFO(10)=1
!
!        INFO(11) --DDASSL normally requires the initial T,
!               Y, and YPRIME to be consistent. That is,
!               you must have G(T,Y,YPRIME) = 0 at the initial
!               time. If you do not know the initial
!               derivative precisely, you can let DDASSL try
!               to compute it.
!          ****   Are the initial T, Y, YPRIME consistent?
!                 Yes - Set INFO(11) = 0
!                  No - Set INFO(11) = 1,
!                       and set YPRIME to an initial approximation
!                       to YPRIME.  (If you have no idea what
!                       YPRIME should be, set it to zero. Note
!                       that the initial Y should be such
!                       that there must exist a YPRIME so that
!                       G(T,Y,YPRIME) = 0.)
!
!  RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL
!         error tolerances to tell the code how accurately you
!         want the solution to be computed.  They must be defined
!         as variables because the code may change them.  You
!         have two choices --
!               Both RTOL and ATOL are scalars. (INFO(2)=0)
!               Both RTOL and ATOL are vectors. (INFO(2)=1)
!         in either case all components must be non-negative.
!
!         The tolerances are used by the code in a local error
!         test at each step which requires roughly that
!               ABS(LOCAL ERROR)  <=  RTOL*ABS(Y)+ATOL
!         for each vector component.
!         (More specifically, a root-mean-square norm is used to
!         measure the size of vectors, and the error test uses the
!         magnitude of the solution at the beginning of the step.)
!
!         The true (global) error is the difference between the
!         true solution of the initial value problem and the
!         computed approximation.  Practically all present day
!         codes, including this one, control the local error at
!         each step and do not even attempt to control the global
!         error directly.
!         Usually, but not always, the true accuracy of the
!         computed Y is comparable to the error tolerances. This
!         code will usually, but not always, deliver a more
!         accurate solution if you reduce the tolerances and
!         integrate again.  By comparing two such solutions you
!         can get a fairly reliable idea of the true error in the
!         solution at the bigger tolerances.
!
!         Setting ATOL=0. results in a pure relative error test on
!         that component.  Setting RTOL=0. results in a pure
!         absolute error test on that component.  A mixed test
!         with non-zero RTOL and ATOL corresponds roughly to a
!         relative error test when the solution component is much
!         bigger than ATOL and to an absolute error test when the
!         solution component is smaller than the threshhold ATOL.
!
!         The code will not attempt to compute a solution at an
!         accuracy unreasonable for the machine being used.  It will
!         advise you if you ask for too much accuracy and inform
!         you as to the maximum accuracy it believes possible.
!
!  RWORK(*) --  Dimension this real work array of length LRW in your
!         calling program.
!
!  LRW -- Set it to the declared length of the RWORK array.
!               You must have
!                    LRW  >=  40+(MAXORD+4)*NEQ+NEQ**2
!               for the full (dense) JACOBIAN case (when INFO(6)=0), or
!                    LRW  >=  40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
!               for the banded user-defined JACOBIAN case
!               (when INFO(5)=1 and INFO(6)=1), or
!                     LRW  >=  40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
!                           +2*(NEQ/(ML+MU+1)+1)
!               for the banded finite-difference-generated JACOBIAN case
!               (when INFO(5)=0 and INFO(6)=1)
!
!  IWORK(*) --  Dimension this integer work array of length LIW in
!         your calling program.
!
!  LIW -- Set it to the declared length of the IWORK array.
!               You must have LIW  >=  20+NEQ
!
!  RPAR, IPAR -- These are parameter arrays, of real and integer
!         type, respectively.  You can use them for communication
!         between your program that calls DDASSL and the
!         RES subroutine (and the JAC subroutine).  They are not
!         altered by DDASSL.  If you do not need RPAR or IPAR,
!         ignore these parameters by treating them as dummy
!         arguments.  If you do choose to use them, dimension
!         them in your calling program and in RES (and in JAC)
!         as arrays of appropriate length.
!
!  JAC -- If you have set INFO(5)=0, you can ignore this parameter
!         by treating it as a dummy argument.  Otherwise, you must
!         provide a subroutine of the form
!             SUBROUTINE JAC(T,Y,YPRIME,PD,CJ,RPAR,IPAR)
!         to define the matrix of partial derivatives
!             PD=DG/DY+CJ*DG/DYPRIME
!         CJ is a scalar which is input to JAC.
!         For the given values of T,Y,YPRIME, the
!         subroutine must evaluate the non-zero partial
!         derivatives for each equation and each solution
!         component, and store these values in the
!         matrix PD.  The elements of PD are set to zero
!         before each call to JAC so only non-zero elements
!         need to be defined.
!
!         Subroutine JAC must not alter T,Y,(*),YPRIME(*), or CJ.
!         You must declare the name JAC in an EXTERNAL statement in
!         your program that calls DDASSL.  You must dimension Y,
!         YPRIME and PD in JAC.
!
!         The way you must store the elements into the PD matrix
!         depends on the structure of the matrix which you
!         indicated by INFO(6).
!               *** INFO(6)=0 -- Full (dense) matrix ***
!                   Give PD a first dimension of NEQ.
!                   When you evaluate the (non-zero) partial derivative
!                   of equation I with respect to variable J, you must
!                   store it in PD according to
!                   PD(I,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)"
!               *** INFO(6)=1 -- Banded JACOBIAN with ML lower and MU
!                   upper diagonal bands (refer to INFO(6) description
!                   of ML and MU) ***
!                   Give PD a first dimension of 2*ML+MU+1.
!                   when you evaluate the (non-zero) partial derivative
!                   of equation I with respect to variable J, you must
!                   store it in PD according to
!                   IROW = I - J + ML + MU + 1
!                   PD(IROW,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)"
!
!         RPAR and IPAR are real and integer parameter arrays
!         which you can use for communication between your calling
!         program and your JACOBIAN subroutine JAC. They are not
!         altered by DDASSL. If you do not need RPAR or IPAR,
!         ignore these parameters by treating them as dummy
!         arguments. If you do choose to use them, dimension
!         them in your calling program and in JAC as arrays of
!         appropriate length.
!
!
!  OPTIONALLY REPLACEABLE NORM ROUTINE:
!
!     DDASSL uses a weighted norm DDANRM to measure the size
!     of vectors such as the estimated error in each step.
!     A FUNCTION subprogram
!       DOUBLE PRECISION FUNCTION DDANRM(NEQ,V,WT,RPAR,IPAR)
!       DIMENSION V(NEQ),WT(NEQ)
!     is used to define this norm. Here, V is the vector
!     whose norm is to be computed, and WT is a vector of
!     weights.  A DDANRM routine has been included with DDASSL
!     which computes the weighted root-mean-square norm
!     given by
!       DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
!     this norm is suitable for most problems. In some
!     special cases, it may be more convenient and/or
!     efficient to define your own norm by writing a function
!     subprogram to be called instead of DDANRM. This should,
!     however, be attempted only after careful thought and
!     consideration.
!
!
!  -------- OUTPUT -- AFTER ANY RETURN FROM DDASSL ---------------------
!
!  The principal aim of the code is to return a computed solution at
!  TOUT, although it is also possible to obtain intermediate results
!  along the way. To find out whether the code achieved its goal
!  or if the integration process was interrupted before the task was
!  completed, you must check the IDID parameter.
!
!
!  T -- The solution was successfully advanced to the
!               output value of T.
!
!  Y(*) -- Contains the computed solution approximation at T.
!
!  YPRIME(*) -- Contains the computed derivative
!               approximation at T.
!
!  IDID -- Reports what the code did.
!
!                     *** Task completed ***
!                Reported by positive values of IDID
!
!           IDID = 1 -- A step was successfully taken in the
!                   intermediate-output mode. The code has not
!                   yet reached TOUT.
!
!           IDID = 2 -- The integration to TSTOP was successfully
!                   completed (T=TSTOP) by stepping exactly to TSTOP.
!
!           IDID = 3 -- The integration to TOUT was successfully
!                   completed (T=TOUT) by stepping past TOUT.
!                   Y(*) is obtained by interpolation.
!                   YPRIME(*) is obtained by interpolation.
!
!                    *** Task interrupted ***
!                Reported by negative values of IDID
!
!           IDID = -1 -- A large amount of work has been expended.
!                   (About 500 steps)
!
!           IDID = -2 -- The error tolerances are too stringent.
!
!           IDID = -3 -- The local error test cannot be satisfied
!                   because you specified a zero component in ATOL
!                   and the corresponding computed solution
!                   component is zero. Thus, a pure relative error
!                   test is impossible for this component.
!
!           IDID = -6 -- DDASSL had repeated error test
!                   failures on the last attempted step.
!
!           IDID = -7 -- The corrector could not converge.
!
!           IDID = -8 -- The matrix of partial derivatives
!                   is singular.
!
!           IDID = -9 -- The corrector could not converge.
!                   there were repeated error test failures
!                   in this step.
!
!           IDID =-10 -- The corrector could not converge
!                   because IRES was equal to minus one.
!
!           IDID =-11 -- IRES equal to -2 was encountered
!                   and control is being returned to the
!                   calling program.
!
!           IDID =-12 -- DDASSL failed to compute the initial
!                   YPRIME.
!
!
!
!           IDID = -13,..,-32 -- Not applicable for this code
!
!                    *** Task terminated ***
!                Reported by the value of IDID=-33
!
!           IDID = -33 -- The code has encountered trouble from which
!                   it cannot recover. A message is printed
!                   explaining the trouble and control is returned
!                   to the calling program. For example, this occurs
!                   when invalid input is detected.
!
!  RTOL, ATOL -- These quantities remain unchanged except when
!               IDID = -2. In this case, the error tolerances have been
!               increased by the code to values which are estimated to
!               be appropriate for continuing the integration. However,
!               the reported solution at T was obtained using the input
!               values of RTOL and ATOL.
!
!  RWORK, IWORK -- Contain information which is usually of no
!               interest to the user but necessary for subsequent calls.
!               However, you may find use for
!
!               RWORK(3)--Which contains the step size H to be
!                       attempted on the next step.
!
!               RWORK(4)--Which contains the current value of the
!                       independent variable, i.e., the farthest point
!                       integration has reached. This will be different
!                       from T only when interpolation has been
!                       performed (IDID=3).
!
!               RWORK(7)--Which contains the stepsize used
!                       on the last successful step.
!
!               IWORK(7)--Which contains the order of the method to
!                       be attempted on the next step.
!
!               IWORK(8)--Which contains the order of the method used
!                       on the last step.
!
!               IWORK(11)--Which contains the number of steps taken so
!                        far.
!
!               IWORK(12)--Which contains the number of calls to RES
!                        so far.
!
!               IWORK(13)--Which contains the number of evaluations of
!                        the matrix of partial derivatives needed so
!                        far.
!
!               IWORK(14)--Which contains the total number
!                        of error test failures so far.
!
!               IWORK(15)--Which contains the total number
!                        of convergence test failures so far.
!                        (includes singular iteration matrix
!                        failures.)
!
!
!  -------- INPUT -- WHAT TO DO TO CONTINUE THE INTEGRATION ------------
!                    (CALLS AFTER THE FIRST)
!
!  This code is organized so that subsequent calls to continue the
!  integration involve little (if any) additional effort on your
!  part. You must monitor the IDID parameter in order to determine
!  what to do next.
!
!  Recalling that the principal task of the code is to integrate
!  from T to TOUT (the interval mode), usually all you will need
!  to do is specify a new TOUT upon reaching the current TOUT.
!
!  Do not alter any quantity not specifically permitted below,
!  in particular do not alter NEQ,T,Y(*),YPRIME(*),RWORK(*),IWORK(*)
!  or the differential equation in subroutine RES. Any such
!  alteration constitutes a new problem and must be treated as such,
!  i.e., you must start afresh.
!
!  You cannot change from vector to scalar error control or vice
!  versa (INFO(2)), but you can change the size of the entries of
!  RTOL, ATOL. Increasing a tolerance makes the equation easier
!  to integrate. Decreasing a tolerance will make the equation
!  harder to integrate and should generally be avoided.
!
!  You can switch from the intermediate-output mode to the
!  interval mode (INFO(3)) or vice versa at any time.
!
!  If it has been necessary to prevent the integration from going
!  past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
!  code will not integrate to any TOUT beyond the currently
!  specified TSTOP. Once TSTOP has been reached you must change
!  the value of TSTOP or set INFO(4)=0. You may change INFO(4)
!  or TSTOP at any time but you must supply the value of TSTOP in
!  RWORK(1) whenever you set INFO(4)=1.
!
!  Do not change INFO(5), INFO(6), IWORK(1), or IWORK(2)
!  unless you are going to restart the code.
!
!                 *** Following a completed task ***
!  If
!     IDID = 1, call the code again to continue the integration
!                  another step in the direction of TOUT.
!
!     IDID = 2 or 3, define a new TOUT and call the code again.
!                  TOUT must be different from T. You cannot change
!                  the direction of integration without restarting.
!
!                 *** Following an interrupted task ***
!               To show the code that you realize the task was
!               interrupted and that you want to continue, you
!               must take appropriate action and set INFO(1) = 1
!  If
!    IDID = -1, The code has taken about 500 steps.
!                  If you want to continue, set INFO(1) = 1 and
!                  call the code again. An additional 500 steps
!                  will be allowed.
!
!    IDID = -2, The error tolerances RTOL, ATOL have been
!                  increased to values the code estimates appropriate
!                  for continuing. You may want to change them
!                  yourself. If you are sure you want to continue
!                  with relaxed error tolerances, set INFO(1)=1 and
!                  call the code again.
!
!    IDID = -3, A solution component is zero and you set the
!                  corresponding component of ATOL to zero. If you
!                  are sure you want to continue, you must first
!                  alter the error criterion to use positive values
!                  for those components of ATOL corresponding to zero
!                  solution components, then set INFO(1)=1 and call
!                  the code again.
!
!    IDID = -4,-5  --- Cannot occur with this code.
!
!    IDID = -6, Repeated error test failures occurred on the
!                  last attempted step in DDASSL. A singularity in the
!                  solution may be present. If you are absolutely
!                  certain you want to continue, you should restart
!                  the integration. (Provide initial values of Y and
!                  YPRIME which are consistent)
!
!    IDID = -7, Repeated convergence test failures occurred
!                  on the last attempted step in DDASSL. An inaccurate
!                  or ill-conditioned JACOBIAN may be the problem. If
!                  you are absolutely certain you want to continue, you
!                  should restart the integration.
!
!    IDID = -8, The matrix of partial derivatives is singular.
!                  Some of your equations may be redundant.
!                  DDASSL cannot solve the problem as stated.
!                  It is possible that the redundant equations
!                  could be removed, and then DDASSL could
!                  solve the problem. It is also possible
!                  that a solution to your problem either
!                  does not exist or is not unique.
!
!    IDID = -9, DDASSL had multiple convergence test
!                  failures, preceded by multiple error
!                  test failures, on the last attempted step.
!                  It is possible that your problem
!                  is ill-posed, and cannot be solved
!                  using this code. Or, there may be a
!                  discontinuity or a singularity in the
!                  solution. If you are absolutely certain
!                  you want to continue, you should restart
!                  the integration.
!
!    IDID =-10, DDASSL had multiple convergence test failures
!                  because IRES was equal to minus one.
!                  If you are absolutely certain you want
!                  to continue, you should restart the
!                  integration.
!
!    IDID =-11, IRES=-2 was encountered, and control is being
!                  returned to the calling program.
!
!    IDID =-12, DDASSL failed to compute the initial YPRIME.
!                  This could happen because the initial
!                  approximation to YPRIME was not very good, or
!                  if a YPRIME consistent with the initial Y
!                  does not exist. The problem could also be caused
!                  by an inaccurate or singular iteration matrix.
!
!    IDID = -13,..,-32  --- Cannot occur with this code.
!
!
!                 *** Following a terminated task ***
!
!  If IDID= -33, you cannot continue the solution of this problem.
!                  An attempt to do so will result in your
!                  run being terminated.
!
!
!  -------- ERROR MESSAGES ---------------------------------------------
!
!      The SLATEC error print routine XERMSG is called in the event of
!   unsuccessful completion of a task.  Most of these are treated as
!   "recoverable errors", which means that (unless the user has directed
!   otherwise) control will be returned to the calling program for
!   possible action after the message has been printed.
!
!   In the event of a negative value of IDID other than -33, an appro-
!   priate message is printed and the "error number" printed by XERMSG
!   is the value of IDID.  There are quite a number of illegal input
!   errors that can lead to a returned value IDID=-33.  The conditions
!   and their printed "error numbers" are as follows:
!
!   Error number       Condition
!
!        1       Some element of INFO vector is not zero or one.
!        2       NEQ .le. 0
!        3       MAXORD not in range.
!        4       LRW is less than the required length for RWORK.
!        5       LIW is less than the required length for IWORK.
!        6       Some element of RTOL is .lt. 0
!        7       Some element of ATOL is .lt. 0
!        8       All elements of RTOL and ATOL are zero.
!        9       INFO(4)=1 and TSTOP is behind TOUT.
!       10       HMAX .lt. 0.0
!       11       TOUT is behind T.
!       12       INFO(8)=1 and H0=0.0
!       13       Some element of WT is .le. 0.0
!       14       TOUT is too close to T to start integration.
!       15       INFO(4)=1 and TSTOP is behind T.
!       16       --( Not used in this version )--
!       17       ML illegal.  Either .lt. 0 or .gt. NEQ
!       18       MU illegal.  Either .lt. 0 or .gt. NEQ
!       19       TOUT = T.
!
!   If DDASSL is called again without any action taken to remove the
!   cause of an unsuccessful return, XERMSG will be called with a fatal
!   error flag, which will cause unconditional termination of the
!   program.  There are two such fatal errors:
!
!   Error number -998:  The last step was terminated with a negative
!       value of IDID other than -33, and no appropriate action was
!       taken.
!
!   Error number -999:  The previous call was terminated because of
!       illegal input (IDID=-33) and there is illegal input in the
!       present call, as well.  (Suspect infinite loop.)
!
!  ---------------------------------------------------------------------
!
!***REFERENCES  A DESCRIPTION OF DASSL: A DIFFERENTIAL/ALGEBRAIC
!                 SYSTEM SOLVER, L. R. PETZOLD, SAND82-8637,
!                 SANDIA NATIONAL LABORATORIES, SEPTEMBER 1982.
!***ROUTINES CALLED  D1MACH, DDAINI, DDANRM, DDASTP, DDATRP, DDAWTS,
!                    XERMSG
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   880387  Code changes made.  All common statements have been
!           replaced by a DATA statement, which defines pointers into
!           RWORK, and PARAMETER statements which define pointers
!           into IWORK.  As well the documentation has gone through
!           grammatical changes.
!   881005  The prologue has been changed to mixed case.
!           The subordinate routines had revision dates changed to
!           this date, although the documentation for these routines
!           is all upper case.  No code changes.
!   890511  Code changes made.  The DATA statement in the declaration
!           section of DDASSL was replaced with a PARAMETER
!           statement.  Also the statement S = 100.D0 was removed
!           from the top of the Newton iteration in DDASTP.
!           The subordinate routines had revision dates changed to
!           this date.
!   890517  The revision date syntax was replaced with the revision
!           history syntax.  Also the "DECK" comment was added to
!           the top of all subroutines.  These changes are consistent
!           with new SLATEC guidelines.
!           The subordinate routines had revision dates changed to
!           this date.  No code changes.
!   891013  Code changes made.
!           Removed all occurrences of FLOAT or DBLE.  All operations
!           are now performed with "mixed-mode" arithmetic.
!           Also, specific function names were replaced with generic
!           function names to be consistent with new SLATEC guidelines.
!           In particular:
!              Replaced DSQRT with SQRT everywhere.
!              Replaced DABS with ABS everywhere.
!              Replaced DMIN1 with MIN everywhere.
!              Replaced MIN0 with MIN everywhere.
!              Replaced DMAX1 with MAX everywhere.
!              Replaced MAX0 with MAX everywhere.
!              Replaced DSIGN with SIGN everywhere.
!           Also replaced REVISION DATE with REVISION HISTORY in all
!           subordinate routines.
!   901004  Miscellaneous changes to prologue to complete conversion
!           to SLATEC 4.0 format.  No code changes.  (F.N.Fritsch)
!   901009  Corrected GAMS classification code and converted subsidiary
!           routines to 4.0 format.  No code changes.  (F.N.Fritsch)
!   901010  Converted XERRWV calls to XERMSG calls.  (R.Clemens, AFWL)
!   901019  Code changes made.
!           Merged SLATEC 4.0 changes with previous changes made
!           by C. Ulrich.  Below is a history of the changes made by
!           C. Ulrich. (Changes in subsidiary routines are implied
!           by this history)
!           891228  Bug was found and repaired inside the DDASSL
!                   and DDAINI routines.  DDAINI was incorrectly
!                   returning the initial T with Y and YPRIME
!                   computed at T+H.  The routine now returns T+H
!                   rather than the initial T.
!                   Cosmetic changes made to DDASTP.
!           900904  Three modifications were made to fix a bug (inside
!                   DDASSL) re interpolation for continuation calls and
!                   cases where TN is very close to TSTOP:
!
!                   1) In testing for whether H is too large, just
!                      compare H to (TSTOP - TN), rather than
!                      (TSTOP - TN) * (1-4*UROUND), and set H to
!                      TSTOP - TN.  This will force DDASTP to step
!                      exactly to TSTOP under certain situations
!                      (i.e. when H returned from DDASTP would otherwise
!                      take TN beyond TSTOP).
!
!                   2) Inside the DDASTP loop, interpolate exactly to
!                      TSTOP if TN is very close to TSTOP (rather than
!                      interpolating to within roundoff of TSTOP).
!
!                   3) Modified IDID description for IDID = 2 to say
!                      that the solution is returned by stepping exactly
!                      to TSTOP, rather than TOUT.  (In some cases the
!                      solution is actually obtained by extrapolating
!                      over a distance near unit roundoff to TSTOP,
!                      but this small distance is deemed acceptable in
!                      these circumstances.)
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue, removed unreferenced labels,
!           and improved XERMSG calls.  (FNF)
!   901030  Added ERROR MESSAGES section and reworked other sections to
!           be of more uniform format.  (FNF)
!   910624  Fixed minor bug related to HMAX (six lines after label
!           525).  (LRP)
!***END PROLOGUE  DDASSL
!
!**End
!
!     Declare arguments.
!
  INTEGER  NEQ, INFO(15), IDID, LRW, IWORK(*), LIW, IPAR(*)
  DOUBLE PRECISION &
     T, Y(*), YPRIME(*), TOUT, RTOL(*), ATOL(*), RWORK(*), &
     RPAR(*)
  EXTERNAL  RES, JAC
!
!     Declare externals.
!
  EXTERNAL  D1MACH, DDAINI, DDANRM, DDASTP, DDATRP, DDAWTS, XERMSG
  DOUBLE PRECISION  D1MACH, DDANRM
!
!     Declare local variables.
!
  INTEGER  I, ITEMP, LALPHA, LBETA, LCJ, LCJOLD, LCTF, LDELTA, &
     LENIW, LENPD, LENRW, LE, LETF, LGAMMA, LH, LHMAX, LHOLD, LIPVT, &
     LJCALC, LK, LKOLD, LIWM, LML, LMTYPE, LMU, LMXORD, LNJE, LNPD, &
     LNRE, LNS, LNST, LNSTL, LPD, LPHASE, LPHI, LPSI, LROUND, LS, &
     LSIGMA, LTN, LTSTOP, LWM, LWT, MBAND, MSAVE, MXORD, NPD, NTEMP, &
     NZFLG
  DOUBLE PRECISION &
     ATOLI, H, HMAX, HMIN, HO, R, RH, RTOLI, TDIST, TN, TNEXT, &
     TSTOP, UROUND, YPNORM
  LOGICAL  DONE
!       Auxiliary variables for conversion of values to be included in
!       error messages.
  CHARACTER*8  XERN1, XERN2
  CHARACTER*16 XERN3, XERN4
!
!     SET POINTERS INTO IWORK
  PARAMETER (LML=1, LMU=2, LMXORD=3, LMTYPE=4, LNST=11, &
    LNRE=12, LNJE=13, LETF=14, LCTF=15, LNPD=16, &
    LIPVT=21, LJCALC=5, LPHASE=6, LK=7, LKOLD=8, &
    LNS=9, LNSTL=10, LIWM=1)
!
!     SET RELATIVE OFFSET INTO RWORK
  PARAMETER (NPD=1)
!
!     SET POINTERS INTO RWORK
  PARAMETER (LTSTOP=1, LHMAX=2, LH=3, LTN=4, &
    LCJ=5, LCJOLD=6, LHOLD=7, LS=8, LROUND=9, &
    LALPHA=11, LBETA=17, LGAMMA=23, &
    LPSI=29, LSIGMA=35, LDELTA=41)
!
!***FIRST EXECUTABLE STATEMENT  DDASSL
  if ( INFO(1) /= 0)go to 100
!
!-----------------------------------------------------------------------
!     THIS BLOCK IS EXECUTED FOR THE INITIAL call ONLY.
!     IT CONTAINS CHECKING OF INPUTS AND INITIALIZATIONS.
!-----------------------------------------------------------------------
!
!     FIRST CHECK INFO ARRAY TO MAKE SURE ALL ELEMENTS OF INFO
!     ARE EITHER ZERO OR ONE.
  DO 10 I=2,11
     if ( INFO(I) /= 0.AND.INFO(I) /= 1)go to 701
10       CONTINUE
!
  if ( NEQ <= 0)go to 702
!
!     CHECK AND COMPUTE MAXIMUM ORDER
  MXORD=5
  if ( INFO(9) == 0)go to 20
     MXORD=IWORK(LMXORD)
     if ( MXORD < 1.OR.MXORD > 5)go to 703
20       IWORK(LMXORD)=MXORD
!
!     COMPUTE MTYPE,LENPD,LENRW.CHECK ML AND MU.
  if ( INFO(6) /= 0)go to 40
     LENPD=NEQ**2
     LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
     if ( INFO(5) /= 0)go to 30
        IWORK(LMTYPE)=2
        go to 60
30          IWORK(LMTYPE)=1
        go to 60
40    if ( IWORK(LML) < 0.OR.IWORK(LML) >= NEQ)go to 717
  if ( IWORK(LMU) < 0.OR.IWORK(LMU) >= NEQ)go to 718
  LENPD=(2*IWORK(LML)+IWORK(LMU)+1)*NEQ
  if ( INFO(5) /= 0)go to 50
     IWORK(LMTYPE)=5
     MBAND=IWORK(LML)+IWORK(LMU)+1
     MSAVE=(NEQ/MBAND)+1
     LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD+2*MSAVE
     go to 60
50       IWORK(LMTYPE)=4
     LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
!
!     CHECK LENGTHS OF RWORK AND IWORK
60    LENIW=20+NEQ
  IWORK(LNPD)=LENPD
  if ( LRW < LENRW)go to 704
  if ( LIW < LENIW)go to 705
!
!     CHECK TO SEE THAT TOUT IS DIFFERENT FROM T
  if ( TOUT  ==  T)go to 719
!
!     CHECK HMAX
  if ( INFO(7) == 0)go to 70
     HMAX=RWORK(LHMAX)
     if ( HMAX <= 0.0D0)go to 710
70    CONTINUE
!
!     INITIALIZE COUNTERS
  IWORK(LNST)=0
  IWORK(LNRE)=0
  IWORK(LNJE)=0
!
  IWORK(LNSTL)=0
  IDID=1
  go to 200
!
!-----------------------------------------------------------------------
!     THIS BLOCK IS FOR CONTINUATION CALLS
!     ONLY. HERE WE CHECK INFO(1), AND if THE
!     LAST STEP WAS INTERRUPTED WE CHECK WHETHER
!     APPROPRIATE ACTION WAS TAKEN.
!-----------------------------------------------------------------------
!
100   CONTINUE
  if ( INFO(1) == 1)go to 110
  if ( INFO(1) /= -1)go to 701
!
!     if WE ARE HERE, THE LAST STEP WAS INTERRUPTED
!     BY AN ERROR CONDITION FROM DDASTP, AND
!     APPROPRIATE ACTION WAS NOT TAKEN. THIS
!     IS A FATAL ERROR.
  WRITE (XERN1, '(I8)') IDID
  call XERMSG ('SLATEC', 'DDASSL', &
     'THE LAST STEP TERMINATED WITH A NEGATIVE VALUE OF IDID = ' // &
     XERN1 // ' AND NO APPROPRIATE ACTION WAS TAKEN.  ' // &
     'RUN TERMINATED', -998, 2)
  return
110   CONTINUE
  IWORK(LNSTL)=IWORK(LNST)
!
!-----------------------------------------------------------------------
!     THIS BLOCK IS EXECUTED ON ALL CALLS.
!     THE ERROR TOLERANCE PARAMETERS ARE
!     CHECKED, AND THE WORK ARRAY POINTERS
!     ARE SET.
!-----------------------------------------------------------------------
!
200   CONTINUE
!     CHECK RTOL,ATOL
  NZFLG=0
  RTOLI=RTOL(1)
  ATOLI=ATOL(1)
  DO 210 I=1,NEQ
     if ( INFO(2) == 1)RTOLI=RTOL(I)
     if ( INFO(2) == 1)ATOLI=ATOL(I)
     if ( RTOLI > 0.0D0.OR.ATOLI > 0.0D0)NZFLG=1
     if ( RTOLI < 0.0D0)go to 706
     if ( ATOLI < 0.0D0)go to 707
210      CONTINUE
  if ( NZFLG == 0)go to 708
!
!     SET UP RWORK STORAGE.IWORK STORAGE IS FIXED
!     IN DATA STATEMENT.
  LE=LDELTA+NEQ
  LWT=LE+NEQ
  LPHI=LWT+NEQ
  LPD=LPHI+(IWORK(LMXORD)+1)*NEQ
  LWM=LPD
  NTEMP=NPD+IWORK(LNPD)
  if ( INFO(1) == 1)go to 400
!
!-----------------------------------------------------------------------
!     THIS BLOCK IS EXECUTED ON THE INITIAL CALL
!     ONLY. SET THE INITIAL STEP SIZE, AND
!     THE ERROR WEIGHT VECTOR, AND PHI.
!     COMPUTE INITIAL YPRIME, if NECESSARY.
!-----------------------------------------------------------------------
!
  TN=T
  IDID=1
!
!     SET ERROR WEIGHT VECTOR WT
  call DDAWTS(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
  DO 305 I = 1,NEQ
     if ( RWORK(LWT+I-1) <= 0.0D0) go to 713
305      CONTINUE
!
!     COMPUTE UNIT ROUNDOFF AND HMIN
  UROUND = D1MACH(4)
  RWORK(LROUND) = UROUND
  HMIN = 4.0D0*UROUND*MAX(ABS(T),ABS(TOUT))
!
!     CHECK INITIAL INTERVAL TO SEE THAT IT IS LONG ENOUGH
  TDIST = ABS(TOUT - T)
  if ( TDIST  <  HMIN) go to 714
!
!     CHECK HO, if THIS WAS INPUT
  if (INFO(8)  ==  0) go to 310
     HO = RWORK(LH)
     if ((TOUT - T)*HO  <  0.0D0) go to 711
     if (HO  ==  0.0D0) go to 712
     go to 320
310    CONTINUE
!
!     COMPUTE INITIAL STEPSIZE, TO BE USED BY EITHER
!     DDASTP OR DDAINI, DEPENDING ON INFO(11)
  HO = 0.001D0*TDIST
  YPNORM = DDANRM(NEQ,YPRIME,RWORK(LWT),RPAR,IPAR)
  if (YPNORM  >  0.5D0/HO) HO = 0.5D0/YPNORM
  HO = SIGN(HO,TOUT-T)
!     ADJUST HO if NECESSARY TO MEET HMAX BOUND
320   if (INFO(7)  ==  0) go to 330
     RH = ABS(HO)/RWORK(LHMAX)
     if (RH  >  1.0D0) HO = HO/RH
!     COMPUTE TSTOP, if APPLICABLE
330   if (INFO(4)  ==  0) go to 340
     TSTOP = RWORK(LTSTOP)
     if ((TSTOP - T)*HO  <  0.0D0) go to 715
     if ((T + HO - TSTOP)*HO  >  0.0D0) HO = TSTOP - T
     if ((TSTOP - TOUT)*HO  <  0.0D0) go to 709
!
!     COMPUTE INITIAL DERIVATIVE, UPDATING TN AND Y, if APPLICABLE
340   if (INFO(11)  ==  0) go to 350
  call DDAINI(TN,Y,YPRIME,NEQ, &
    RES,JAC,HO,RWORK(LWT),IDID,RPAR,IPAR, &
    RWORK(LPHI),RWORK(LDELTA),RWORK(LE), &
    RWORK(LWM),IWORK(LIWM),HMIN,RWORK(LROUND), &
    INFO(10),NTEMP)
  if (IDID  <  0) go to 390
!
!     LOAD H WITH HO.  STORE H IN RWORK(LH)
350   H = HO
  RWORK(LH) = H
!
!     LOAD Y AND H*YPRIME INTO PHI(*,1) AND PHI(*,2)
  ITEMP = LPHI + NEQ
  DO 370 I = 1,NEQ
     RWORK(LPHI + I - 1) = Y(I)
370      RWORK(ITEMP + I - 1) = H*YPRIME(I)
!
390   go to 500
!
!-------------------------------------------------------
!     THIS BLOCK IS FOR CONTINUATION CALLS ONLY. ITS
!     PURPOSE IS TO CHECK STOP CONDITIONS BEFORE
!     TAKING A STEP.
!     ADJUST H if NECESSARY TO MEET HMAX BOUND
!-------------------------------------------------------
!
400   CONTINUE
  UROUND=RWORK(LROUND)
  DONE = .FALSE.
  TN=RWORK(LTN)
  H=RWORK(LH)
  if ( INFO(7)  ==  0) go to 410
     RH = ABS(H)/RWORK(LHMAX)
     if ( RH  >  1.0D0) H = H/RH
410   CONTINUE
  if ( T  ==  TOUT) go to 719
  if ( (T - TOUT)*H  >  0.0D0) go to 711
  if ( INFO(4)  ==  1) go to 430
  if ( INFO(3)  ==  1) go to 420
  if ( (TN-TOUT)*H < 0.0D0)go to 490
  call DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD), &
    RWORK(LPHI),RWORK(LPSI))
  T=TOUT
  IDID = 3
  DONE = .TRUE.
  go to 490
420   if ( (TN-T)*H  <=  0.0D0) go to 490
  if ( (TN - TOUT)*H  >  0.0D0) go to 425
  call DDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD), &
    RWORK(LPHI),RWORK(LPSI))
  T = TN
  IDID = 1
  DONE = .TRUE.
  go to 490
425   CONTINUE
  call DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD), &
    RWORK(LPHI),RWORK(LPSI))
  T = TOUT
  IDID = 3
  DONE = .TRUE.
  go to 490
430   if ( INFO(3)  ==  1) go to 440
  TSTOP=RWORK(LTSTOP)
  if ( (TN-TSTOP)*H > 0.0D0) go to 715
  if ( (TSTOP-TOUT)*H < 0.0D0)go to 709
  if ( (TN-TOUT)*H < 0.0D0)go to 450
  call DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD), &
     RWORK(LPHI),RWORK(LPSI))
  T=TOUT
  IDID = 3
  DONE = .TRUE.
  go to 490
440   TSTOP = RWORK(LTSTOP)
  if ( (TN-TSTOP)*H  >  0.0D0) go to 715
  if ( (TSTOP-TOUT)*H  <  0.0D0) go to 709
  if ( (TN-T)*H  <=  0.0D0) go to 450
  if ( (TN - TOUT)*H  >  0.0D0) go to 445
  call DDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD), &
    RWORK(LPHI),RWORK(LPSI))
  T = TN
  IDID = 1
  DONE = .TRUE.
  go to 490
445   CONTINUE
  call DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD), &
    RWORK(LPHI),RWORK(LPSI))
  T = TOUT
  IDID = 3
  DONE = .TRUE.
  go to 490
450   CONTINUE
!     CHECK WHETHER WE ARE WITHIN ROUNDOFF OF TSTOP
  if ( ABS(TN-TSTOP) > 100.0D0*UROUND* &
     (ABS(TN)+ABS(H)))go to 460
  call DDATRP(TN,TSTOP,Y,YPRIME,NEQ,IWORK(LKOLD), &
    RWORK(LPHI),RWORK(LPSI))
  IDID=2
  T=TSTOP
  DONE = .TRUE.
  go to 490
460   TNEXT=TN+H
  if ( (TNEXT-TSTOP)*H <= 0.0D0)go to 490
  H=TSTOP-TN
  RWORK(LH)=H
!
490   if (DONE) go to 580
!
!-------------------------------------------------------
!     THE NEXT BLOCK CONTAINS THE call TO THE
!     ONE-STEP INTEGRATOR DDASTP.
!     THIS IS A LOOPING POINT FOR THE INTEGRATION STEPS.
!     CHECK FOR TOO MANY STEPS.
!     UPDATE WT.
!     CHECK FOR TOO MUCH ACCURACY REQUESTED.
!     COMPUTE MINIMUM STEPSIZE.
!-------------------------------------------------------
!
500   CONTINUE
!     CHECK FOR FAILURE TO COMPUTE INITIAL YPRIME
  if (IDID  ==  -12) go to 527
!
!     CHECK FOR TOO MANY STEPS
  if ( (IWORK(LNST)-IWORK(LNSTL)) < 500) &
     go to 510
       IDID=-1
       go to 527
!
!     UPDATE WT
510   call DDAWTS(NEQ,INFO(2),RTOL,ATOL,RWORK(LPHI), &
    RWORK(LWT),RPAR,IPAR)
  DO 520 I=1,NEQ
     if ( RWORK(I+LWT-1) > 0.0D0)go to 520
       IDID=-3
       go to 527
520   CONTINUE
!
!     TEST FOR TOO MUCH ACCURACY REQUESTED.
  R=DDANRM(NEQ,RWORK(LPHI),RWORK(LWT),RPAR,IPAR)* &
     100.0D0*UROUND
  if ( R <= 1.0D0)go to 525
!     MULTIPLY RTOL AND ATOL BY R AND RETURN
  if ( INFO(2) == 1)go to 523
       RTOL(1)=R*RTOL(1)
       ATOL(1)=R*ATOL(1)
       IDID=-2
       go to 527
523   DO 524 I=1,NEQ
       RTOL(I)=R*RTOL(I)
524        ATOL(I)=R*ATOL(I)
  IDID=-2
  go to 527
525   CONTINUE
!
!     COMPUTE MINIMUM STEPSIZE
  HMIN=4.0D0*UROUND*MAX(ABS(TN),ABS(TOUT))
!
!     TEST H VS. HMAX
  if (INFO(7)  /=  0) THEN
     RH = ABS(H)/RWORK(LHMAX)
     if (RH  >  1.0D0) H = H/RH
  end if
!
  call DDASTP(TN,Y,YPRIME,NEQ, &
     RES,JAC,H,RWORK(LWT),INFO(1),IDID,RPAR,IPAR, &
     RWORK(LPHI),RWORK(LDELTA),RWORK(LE), &
     RWORK(LWM),IWORK(LIWM), &
     RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA), &
     RWORK(LPSI),RWORK(LSIGMA), &
     RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD), &
     RWORK(LS),HMIN,RWORK(LROUND), &
     IWORK(LPHASE),IWORK(LJCALC),IWORK(LK), &
     IWORK(LKOLD),IWORK(LNS),INFO(10),NTEMP)
527   if ( IDID < 0)go to 600
!
!--------------------------------------------------------
!     THIS BLOCK HANDLES THE CASE OF A SUCCESSFUL RETURN
!     FROM DDASTP (IDID=1).  TEST FOR STOP CONDITIONS.
!--------------------------------------------------------
!
  if ( INFO(4) /= 0)go to 540
       if ( INFO(3) /= 0)go to 530
         if ( (TN-TOUT)*H < 0.0D0)go to 500
         call DDATRP(TN,TOUT,Y,YPRIME,NEQ, &
           IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
         IDID=3
         T=TOUT
         go to 580
530          if ( (TN-TOUT)*H >= 0.0D0)go to 535
         T=TN
         IDID=1
         go to 580
535          call DDATRP(TN,TOUT,Y,YPRIME,NEQ, &
           IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
         IDID=3
         T=TOUT
         go to 580
540   if ( INFO(3) /= 0)go to 550
  if ( (TN-TOUT)*H < 0.0D0)go to 542
     call DDATRP(TN,TOUT,Y,YPRIME,NEQ, &
       IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
     T=TOUT
     IDID=3
     go to 580
542   if ( ABS(TN-TSTOP) <= 100.0D0*UROUND* &
     (ABS(TN)+ABS(H)))go to 545
  TNEXT=TN+H
  if ( (TNEXT-TSTOP)*H <= 0.0D0)go to 500
  H=TSTOP-TN
  go to 500
545   call DDATRP(TN,TSTOP,Y,YPRIME,NEQ, &
    IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
  IDID=2
  T=TSTOP
  go to 580
550   if ( (TN-TOUT)*H >= 0.0D0)go to 555
  if ( ABS(TN-TSTOP) <= 100.0D0*UROUND*(ABS(TN)+ABS(H)))go to 552
  T=TN
  IDID=1
  go to 580
552   call DDATRP(TN,TSTOP,Y,YPRIME,NEQ, &
    IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
  IDID=2
  T=TSTOP
  go to 580
555   call DDATRP(TN,TOUT,Y,YPRIME,NEQ, &
     IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
  T=TOUT
  IDID=3
  go to 580
!
!--------------------------------------------------------
!     ALL SUCCESSFUL RETURNS FROM DDASSL ARE MADE FROM
!     THIS BLOCK.
!--------------------------------------------------------
!
580   CONTINUE
  RWORK(LTN)=TN
  RWORK(LH)=H
  return
!
!-----------------------------------------------------------------------
!     THIS BLOCK HANDLES ALL UNSUCCESSFUL
!     returnS OTHER THAN FOR ILLEGAL INPUT.
!-----------------------------------------------------------------------
!
600   CONTINUE
  ITEMP=-IDID
  go to (610,620,630,690,690,640,650,660,670,675, &
    680,685), ITEMP
!
!     THE MAXIMUM NUMBER OF STEPS WAS TAKEN BEFORE
!     REACHING TOUT
610   WRITE (XERN3, '(1P,D15.6)') TN
  call XERMSG ('SLATEC', 'DDASSL', &
     'AT CURRENT T = ' // XERN3 // ' 500 STEPS TAKEN ON THIS ' // &
     'CALL BEFORE REACHING TOUT', IDID, 1)
  go to 690
!
!     TOO MUCH ACCURACY FOR MACHINE PRECISION
620   WRITE (XERN3, '(1P,D15.6)') TN
  call XERMSG ('SLATEC', 'DDASSL', &
     'AT T = ' // XERN3 // ' TOO MUCH ACCURACY REQUESTED FOR ' // &
     'PRECISION OF MACHINE. RTOL AND ATOL WERE INCREASED TO ' // &
     'APPROPRIATE VALUES', IDID, 1)
  go to 690
!
!     WT(I)  <=  0.0 FOR SOME I (NOT AT START OF PROBLEM)
630   WRITE (XERN3, '(1P,D15.6)') TN
  call XERMSG ('SLATEC', 'DDASSL', &
     'AT T = ' // XERN3 // ' SOME ELEMENT OF WT HAS BECOME  <=  ' // &
     '0.0', IDID, 1)
  go to 690
!
!     ERROR TEST FAILED REPEATEDLY OR WITH H=HMIN
640   WRITE (XERN3, '(1P,D15.6)') TN
  WRITE (XERN4, '(1P,D15.6)') H
  call XERMSG ('SLATEC', 'DDASSL', &
     'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 // &
     ' THE ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN', &
     IDID, 1)
  go to 690
!
!     CORRECTOR CONVERGENCE FAILED REPEATEDLY OR WITH H=HMIN
650   WRITE (XERN3, '(1P,D15.6)') TN
  WRITE (XERN4, '(1P,D15.6)') H
  call XERMSG ('SLATEC', 'DDASSL', &
     'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 // &
     ' THE CORRECTOR FAILED TO CONVERGE REPEATEDLY OR WITH ' // &
     'ABS(H)=HMIN', IDID, 1)
  go to 690
!
!     THE ITERATION MATRIX IS SINGULAR
660   WRITE (XERN3, '(1P,D15.6)') TN
  WRITE (XERN4, '(1P,D15.6)') H
  call XERMSG ('SLATEC', 'DDASSL', &
     'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 // &
     ' THE ITERATION MATRIX IS SINGULAR', IDID, 1)
  go to 690
!
!     CORRECTOR FAILURE PRECEDED BY ERROR TEST FAILURES.
670   WRITE (XERN3, '(1P,D15.6)') TN
  WRITE (XERN4, '(1P,D15.6)') H
  call XERMSG ('SLATEC', 'DDASSL', &
     'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 // &
     ' THE CORRECTOR COULD NOT CONVERGE.  ALSO, THE ERROR TEST ' // &
     'FAILED REPEATEDLY.', IDID, 1)
  go to 690
!
!     CORRECTOR FAILURE BECAUSE IRES = -1
675   WRITE (XERN3, '(1P,D15.6)') TN
  WRITE (XERN4, '(1P,D15.6)') H
  call XERMSG ('SLATEC', 'DDASSL', &
     'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 // &
     ' THE CORRECTOR COULD NOT CONVERGE BECAUSE IRES WAS EQUAL ' // &
     'TO MINUS ONE', IDID, 1)
  go to 690
!
!     FAILURE BECAUSE IRES = -2
680   WRITE (XERN3, '(1P,D15.6)') TN
  WRITE (XERN4, '(1P,D15.6)') H
  call XERMSG ('SLATEC', 'DDASSL', &
     'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 // &
     ' IRES WAS EQUAL TO MINUS TWO', IDID, 1)
  go to 690
!
!     FAILED TO COMPUTE INITIAL YPRIME
685   WRITE (XERN3, '(1P,D15.6)') TN
  WRITE (XERN4, '(1P,D15.6)') HO
  call XERMSG ('SLATEC', 'DDASSL', &
     'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 // &
     ' THE INITIAL YPRIME COULD NOT BE COMPUTED', IDID, 1)
  go to 690
!
690   CONTINUE
  INFO(1)=-1
  T=TN
  RWORK(LTN)=TN
  RWORK(LH)=H
  return
!
!-----------------------------------------------------------------------
!     THIS BLOCK HANDLES ALL ERROR RETURNS DUE
!     TO ILLEGAL INPUT, AS DETECTED BEFORE CALLING
!     DDASTP. FIRST THE ERROR MESSAGE ROUTINE IS
!     CALLED. if THIS HAPPENS TWICE IN
!     SUCCESSION, EXECUTION IS TERMINATED
!
!-----------------------------------------------------------------------
701   call XERMSG ('SLATEC', 'DDASSL', &
     'SOME ELEMENT OF INFO VECTOR IS NOT ZERO OR ONE', 1, 1)
  go to 750
!
702   WRITE (XERN1, '(I8)') NEQ
  call XERMSG ('SLATEC', 'DDASSL', &
     'NEQ = ' // XERN1 // '  <=  0', 2, 1)
  go to 750
!
703   WRITE (XERN1, '(I8)') MXORD
  call XERMSG ('SLATEC', 'DDASSL', &
     'MAXORD = ' // XERN1 // ' NOT IN RANGE', 3, 1)
  go to 750
!
704   WRITE (XERN1, '(I8)') LENRW
  WRITE (XERN2, '(I8)') LRW
  call XERMSG ('SLATEC', 'DDASSL', &
     'RWORK LENGTH NEEDED, LENRW = ' // XERN1 // &
     ', EXCEEDS LRW = ' // XERN2, 4, 1)
  go to 750
!
705   WRITE (XERN1, '(I8)') LENIW
  WRITE (XERN2, '(I8)') LIW
  call XERMSG ('SLATEC', 'DDASSL', &
     'IWORK LENGTH NEEDED, LENIW = ' // XERN1 // &
     ', EXCEEDS LIW = ' // XERN2, 5, 1)
  go to 750
!
706   call XERMSG ('SLATEC', 'DDASSL', &
     'SOME ELEMENT OF RTOL IS  <  0', 6, 1)
  go to 750
!
707   call XERMSG ('SLATEC', 'DDASSL', &
     'SOME ELEMENT OF ATOL IS  <  0', 7, 1)
  go to 750
!
708   call XERMSG ('SLATEC', 'DDASSL', &
     'ALL ELEMENTS OF RTOL AND ATOL ARE ZERO', 8, 1)
  go to 750
!
709   WRITE (XERN3, '(1P,D15.6)') TSTOP
  WRITE (XERN4, '(1P,D15.6)') TOUT
  call XERMSG ('SLATEC', 'DDASSL', &
     'INFO(4) = 1 AND TSTOP = ' // XERN3 // ' BEHIND TOUT = ' // &
     XERN4, 9, 1)
  go to 750
!
710   WRITE (XERN3, '(1P,D15.6)') HMAX
  call XERMSG ('SLATEC', 'DDASSL', &
     'HMAX = ' // XERN3 // '  <  0.0', 10, 1)
  go to 750
!
711   WRITE (XERN3, '(1P,D15.6)') TOUT
  WRITE (XERN4, '(1P,D15.6)') T
  call XERMSG ('SLATEC', 'DDASSL', &
     'TOUT = ' // XERN3 // ' BEHIND T = ' // XERN4, 11, 1)
  go to 750
!
712   call XERMSG ('SLATEC', 'DDASSL', &
     'INFO(8)=1 AND H0=0.0', 12, 1)
  go to 750
!
713   call XERMSG ('SLATEC', 'DDASSL', &
     'SOME ELEMENT OF WT IS  <=  0.0', 13, 1)
  go to 750
!
714   WRITE (XERN3, '(1P,D15.6)') TOUT
  WRITE (XERN4, '(1P,D15.6)') T
  call XERMSG ('SLATEC', 'DDASSL', &
     'TOUT = ' // XERN3 // ' TOO CLOSE TO T = ' // XERN4 // &
     ' TO START INTEGRATION', 14, 1)
  go to 750
!
715   WRITE (XERN3, '(1P,D15.6)') TSTOP
  WRITE (XERN4, '(1P,D15.6)') T
  call XERMSG ('SLATEC', 'DDASSL', &
     'INFO(4)=1 AND TSTOP = ' // XERN3 // ' BEHIND T = ' // XERN4, &
     15, 1)
  go to 750
!
717   WRITE (XERN1, '(I8)') IWORK(LML)
  call XERMSG ('SLATEC', 'DDASSL', &
     'ML = ' // XERN1 // ' ILLEGAL.  EITHER  <  0 OR  >  NEQ', &
     17, 1)
  go to 750
!
718   WRITE (XERN1, '(I8)') IWORK(LMU)
  call XERMSG ('SLATEC', 'DDASSL', &
     'MU = ' // XERN1 // ' ILLEGAL.  EITHER  <  0 OR  >  NEQ', &
     18, 1)
  go to 750
!
719   WRITE (XERN3, '(1P,D15.6)') TOUT
  call XERMSG ('SLATEC', 'DDASSL', &
    'TOUT = T = ' // XERN3, 19, 1)
  go to 750
!
750   IDID=-33
  if ( INFO(1) == -1) THEN
     call XERMSG ('SLATEC', 'DDASSL', &
        'REPEATED OCCURRENCES OF ILLEGAL INPUT$$' // &
        'RUN TERMINATED. APPARENT INFINITE LOOP', -999, 2)
  end if
!
  INFO(1)=-1
  return
!-----------END OF SUBROUTINE DDASSL------------------------------------
end
