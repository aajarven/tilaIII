subroutine CDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, &
     EWT, IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK, &
     LENW, IWORK, LENIW, JACOBN, FA, NDE, MXSTEP, G, USERS, IERFLG)
!
!! CDRIV3 solves N ordinary differential equations of the form ...
!  dY(I)/dT = F(Y(I),T), given the initial conditions Y(I) = YI.  The ...
!  program has options to allow the solution of both stiff and non-stiff ...
!  differential equations.  Other important options are available.  CDRIV3
!  allows complex-valued differential equations.
!
!***LIBRARY   SLATEC (SDRIVE)
!***CATEGORY  I1A2, I1A1B
!***TYPE      COMPLEX (SDRIV3-S, DDRIV3-D, CDRIV3-C)
!***KEYWORDS  COMPLEX VALUED, GEAR'S METHOD, INITIAL VALUE PROBLEMS,
!             ODE, ORDINARY DIFFERENTIAL EQUATIONS, SDRIVE, STIFF
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***DESCRIPTION
!
!  I.  ABSTRACT  .......................................................
!
!    The primary function of CDRIV3 is to solve N ordinary differential
!    equations of the form dY(I)/dT = F(Y(I),T), given the initial
!    conditions Y(I) = YI.  The program has options to allow the
!    solution of both stiff and non-stiff differential equations.  In
!    addition, CDRIV3 may be used to solve:
!      1. The initial value problem, A*dY(I)/dT = F(Y(I),T), where A is
!         a non-singular matrix depending on Y and T.
!      2. The hybrid differential/algebraic initial value problem,
!         A*dY(I)/dT = F(Y(I),T), where A is a vector (whose values may
!         depend upon Y and T) some of whose components will be zero
!         corresponding to those equations which are algebraic rather
!         than differential.
!    CDRIV3 is to be called once for each output point of T.
!
!  II.  PARAMETERS  ....................................................
!
!    The user should use parameter names in the call sequence of CDRIV3
!    for those quantities whose value may be altered by CDRIV3.  The
!    parameters in the call sequence are:
!
!    N      = (Input) The number of dependent functions whose solution
!             is desired.  N must not be altered during a problem.
!
!    T      = (Real) The independent variable.  On input for the first
!             call, T is the initial point.  On output, T is the point
!             at which the solution is given.
!
!    Y      = (Complex) The vector of dependent variables.  Y is used as
!             input on the first call, to set the initial values.  On
!             output, Y is the computed solution vector.  This array Y
!             is passed in the call sequence of the user-provided
!             routines F, JACOBN, FA, USERS, and G.  Thus parameters
!             required by those routines can be stored in this array in
!             components N+1 and above.  (Note: Changes by the user to
!             the first N components of this array will take effect only
!             after a restart, i.e., after setting NSTATE to 1 .)
!
!    F      = A subroutine supplied by the user.  The name must be
!             declared EXTERNAL in the user's calling program.  This
!             subroutine is of the form:
!                   SUBROUTINE F (N, T, Y, YDOT)
!                   COMPLEX Y(*), YDOT(*)
!                     .
!                     .
!                   YDOT(1) = ...
!                     .
!                     .
!                   YDOT(N) = ...
!                   END (Sample)
!             This computes YDOT = F(Y,T), the right hand side of the
!             differential equations.  Here Y is a vector of length at
!             least N.  The actual length of Y is determined by the
!             user's declaration in the program which calls CDRIV3.
!             Thus the dimensioning of Y in F, while required by FORTRAN
!             convention, does not actually allocate any storage.  When
!             this subroutine is called, the first N components of Y are
!             intermediate approximations to the solution components.
!             The user should not alter these values.  Here YDOT is a
!             vector of length N.  The user should only compute YDOT(I)
!             for I from 1 to N.  Normally a return from F passes
!             control back to  CDRIV3.  However, if the user would like
!             to abort the calculation, i.e., return control to the
!             program which calls CDRIV3, he should set N to zero.
!             CDRIV3 will signal this by returning a value of NSTATE
!             equal to 6 .  Altering the value of N in F has no effect
!             on the value of N in the call sequence of CDRIV3.
!
!    NSTATE = An integer describing the status of integration.  The
!             meaning of NSTATE is as follows:
!               1  (Input) Means the first call to the routine.  This
!                  value must be set by the user.  On all subsequent
!                  calls the value of NSTATE should be tested by the
!                  user, but must not be altered.  (As a convenience to
!                  the user who may wish to put out the initial
!                  conditions, CDRIV3 can be called with NSTATE=1, and
!                  TOUT=T.  In this case the program will return with
!                  NSTATE unchanged, i.e., NSTATE=1.)
!               2  (Output) Means a successful integration.  If a normal
!                  continuation is desired (i.e., a further integration
!                  in the same direction), simply advance TOUT and call
!                  again.  All other parameters are automatically set.
!               3  (Output)(Unsuccessful) Means the integrator has taken
!                  MXSTEP steps without reaching TOUT.  The user can
!                  continue the integration by simply calling CDRIV3
!                  again.
!               4  (Output)(Unsuccessful) Means too much accuracy has
!                  been requested.  EPS has been increased to a value
!                  the program estimates is appropriate.  The user can
!                  continue the integration by simply calling CDRIV3
!                  again.
!               5  (Output) A root was found at a point less than TOUT.
!                  The user can continue the integration toward TOUT by
!                  simply calling CDRIV3 again.
!               6  (Output)(Unsuccessful) N has been set to zero in
!                  SUBROUTINE F.
!               7  (Output)(Unsuccessful) N has been set to zero in
!                  FUNCTION G.  See description of G below.
!               8  (Output)(Unsuccessful) N has been set to zero in
!                  SUBROUTINE JACOBN.  See description of JACOBN below.
!               9  (Output)(Unsuccessful) N has been set to zero in
!                  SUBROUTINE FA.  See description of FA below.
!              10  (Output)(Unsuccessful) N has been set to zero in
!                  SUBROUTINE USERS.  See description of USERS below.
!              11  (Output)(Successful) For NTASK = 2 or 3, T is beyond
!                  TOUT.  The solution was obtained by interpolation.
!                  The user can continue the integration by simply
!                  advancing TOUT and calling CDRIV3 again.
!              12  (Output)(Unsuccessful) The solution could not be
!                  obtained.  The value of IERFLG (see description
!                  below) for a "Recoverable" situation indicates the
!                  type of difficulty encountered: either an illegal
!                  value for a parameter or an inability to continue the
!                  solution.  For this condition the user should take
!                  corrective action and reset NSTATE to 1 before
!                  calling CDRIV3 again.  Otherwise the program will
!                  terminate the run.
!
!    TOUT   = (Input, Real) The point at which the solution is desired.
!             The position of TOUT relative to T on the first call
!             determines the direction of integration.
!
!    NTASK  = (Input) An index specifying the manner of returning the
!             solution, according to the following:
!               NTASK = 1  Means CDRIV3 will integrate past TOUT and
!                          interpolate the solution.  This is the most
!                          efficient mode.
!               NTASK = 2  Means CDRIV3 will return the solution after
!                          each internal integration step, or at TOUT,
!                          whichever comes first.  In the latter case,
!                          the program integrates exactly to TOUT.
!               NTASK = 3  Means CDRIV3 will adjust its internal step to
!                          reach TOUT exactly (useful if a singularity
!                          exists beyond TOUT.)
!
!    NROOT  = (Input) The number of equations whose roots are desired.
!             If NROOT is zero, the root search is not active.  This
!             option is useful for obtaining output at points which are
!             not known in advance, but depend upon the solution, e.g.,
!             when some solution component takes on a specified value.
!             The root search is carried out using the user-written
!             function G (see description of G below.)  CDRIV3 attempts
!             to find the value of T at which one of the equations
!             changes sign.  CDRIV3 can find at most one root per
!             equation per internal integration step, and will then
!             return the solution either at TOUT or at a root, whichever
!             occurs first in the direction of integration.  The initial
!             point is never reported as a root.  The index of the
!             equation whose root is being reported is stored in the
!             sixth element of IWORK.
!             NOTE: NROOT is never altered by this program.
!
!    EPS    = (Real) On input, the requested relative accuracy in all
!             solution components.  EPS = 0 is allowed.  On output, the
!             adjusted relative accuracy if the input value was too
!             small.  The value of EPS should be set as large as is
!             reasonable, because the amount of work done by CDRIV3
!             increases as EPS decreases.
!
!    EWT    = (Input, Real) Problem zero, i.e., the smallest, nonzero,
!             physically meaningful value for the solution.  (Array,
!             possibly of length one.  See following description of
!             IERROR.)  Setting EWT smaller than necessary can adversely
!             affect the running time.
!
!    IERROR = (Input) Error control indicator.  A value of 3 is
!             suggested for most problems.  Other choices and detailed
!             explanations of EWT and IERROR are given below for those
!             who may need extra flexibility.
!
!             These last three input quantities EPS, EWT and IERROR
!             control the accuracy of the computed solution.  EWT and
!             IERROR are used internally to compute an array YWT.  One
!             step error estimates divided by YWT(I) are kept less than
!             EPS in root mean square norm.
!                 IERROR (Set by the user) =
!                 1  Means YWT(I) = 1. (Absolute error control)
!                                   EWT is ignored.
!                 2  Means YWT(I) = ABS(Y(I)),  (Relative error control)
!                                   EWT is ignored.
!                 3  Means YWT(I) = MAX(ABS(Y(I)), EWT(1)).
!                 4  Means YWT(I) = MAX(ABS(Y(I)), EWT(I)).
!                    This choice is useful when the solution components
!                    have differing scales.
!                 5  Means YWT(I) = EWT(I).
!             If IERROR is 3, EWT need only be dimensioned one.
!             If IERROR is 4 or 5, the user must dimension EWT at least
!             N, and set its values.
!
!    MINT   = (Input) The integration method indicator.
!               MINT = 1  Means the Adams methods, and is used for
!                         non-stiff problems.
!               MINT = 2  Means the stiff methods of Gear (i.e., the
!                         backward differentiation formulas), and is
!                         used for stiff problems.
!               MINT = 3  Means the program dynamically selects the
!                         Adams methods when the problem is non-stiff
!                         and the Gear methods when the problem is
!                         stiff.  When using the Adams methods, the
!                         program uses a value of MITER=0; when using
!                         the Gear methods, the program uses the value
!                         of MITER provided by the user.  Only a value
!                         of IMPL = 0 and a value of MITER = 1, 2, 4, or
!                         5 is allowed for this option.  The user may
!                         not alter the value of MINT or MITER without
!                         restarting, i.e., setting NSTATE to 1.
!
!    MITER  = (Input) The iteration method indicator.
!               MITER = 0  Means functional iteration.  This value is
!                          suggested for non-stiff problems.
!               MITER = 1  Means chord method with analytic Jacobian.
!                          In this case, the user supplies subroutine
!                          JACOBN (see description below).
!               MITER = 2  Means chord method with Jacobian calculated
!                          internally by finite differences.
!               MITER = 3  Means chord method with corrections computed
!                          by the user-written routine USERS (see
!                          description of USERS below.)  This option
!                          allows all matrix algebra and storage
!                          decisions to be made by the user.  When using
!                          a value of MITER = 3, the subroutine FA is
!                          not required, even if IMPL is not 0.  For
!                          further information on using this option, see
!                          Section IV-E below.
!               MITER = 4  Means the same as MITER = 1 but the A and
!                          Jacobian matrices are assumed to be banded.
!               MITER = 5  Means the same as MITER = 2 but the A and
!                          Jacobian matrices are assumed to be banded.
!
!    IMPL   = (Input) The implicit method indicator.
!               IMPL = 0    Means solving dY(I)/dT = F(Y(I),T).
!               IMPL = 1    Means solving A*dY(I)/dT = F(Y(I),T), non-
!                           singular A (see description of FA below.)
!                           Only MINT = 1 or 2, and MITER = 1, 2, 3, 4,
!                           or 5 are allowed for this option.
!               IMPL = 2,3  Means solving certain systems of hybrid
!                           differential/algebraic equations (see
!                           description of FA below.)  Only MINT = 2 and
!                           MITER = 1, 2, 3, 4, or 5, are allowed for
!                           this option.
!               The value of IMPL must not be changed during a problem.
!
!    ML     = (Input) The lower half-bandwidth in the case of a banded
!             A or Jacobian matrix.  (I.e., maximum(R-C) for nonzero
!             A(R,C).)
!
!    MU     = (Input) The upper half-bandwidth in the case of a banded
!             A or Jacobian matrix.  (I.e., maximum(C-R).)
!
!    MXORD  = (Input) The maximum order desired. This is  <=  12 for
!             the Adams methods and  <=  5 for the Gear methods.  Normal
!             value is 12 and 5, respectively.  If MINT is 3, the
!             maximum order used will be MIN(MXORD, 12) when using the
!             Adams methods, and MIN(MXORD, 5) when using the Gear
!             methods.  MXORD must not be altered during a problem.
!
!    HMAX   = (Input, Real) The maximum magnitude of the step size that
!             will be used for the problem.  This is useful for ensuring
!             that important details are not missed.  If this is not the
!             case, a large value, such as the interval length, is
!             suggested.
!
!    WORK
!    LENW   = (Input)
!             WORK is an array of LENW complex words used
!             internally for temporary storage.  The user must allocate
!             space for this array in the calling program by a statement
!             such as
!                       COMPLEX WORK(...)
!             The following table gives the required minimum value for
!             the length of WORK, depending on the value of IMPL and
!             MITER.  LENW should be set to the value used.  The
!             contents of WORK should not be disturbed between calls to
!             CDRIV3.
!
!      IMPL =   0            1               2             3
!              ---------------------------------------------------------
! MITER =  0   (MXORD+4)*N   Not allowed     Not allowed   Not allowed
!              + 2*NROOT
!              + 250
!
!         1,2  N*N +         2*N*N +         N*N +         N*(N + NDE)
!              (MXORD+5)*N   (MXORD+5)*N     (MXORD+6)*N   + (MXORD+5)*N
!              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT
!              + 250         + 250           + 250         + 250
!
!          3   (MXORD+4)*N   (MXORD+4)*N     (MXORD+4)*N   (MXORD+4)*N
!              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT
!              + 250         + 250           + 250         + 250
!
!         4,5  (2*ML+MU+1)   2*(2*ML+MU+1)   (2*ML+MU+1)   (2*ML+MU+1)*
!              *N +          *N +            *N +          (N+NDE) +
!              (MXORD+5)*N   (MXORD+5)*N     (MXORD+6)*N   + (MXORD+5)*N
!              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT
!              + 250         + 250           + 250         + 250
!              ---------------------------------------------------------
!
!    IWORK
!    LENIW  = (Input)
!             IWORK is an integer array of length LENIW used internally
!             for temporary storage.  The user must allocate space for
!             this array in the calling program by a statement such as
!                       INTEGER IWORK(...)
!             The length of IWORK should be at least
!               50      if MITER is 0 or 3, or
!               N+50    if MITER is 1, 2, 4, or 5, or MINT is 3,
!             and LENIW should be set to the value used.  The contents
!             of IWORK should not be disturbed between calls to CDRIV3.
!
!    JACOBN = A subroutine supplied by the user, if MITER is 1 or 4.
!             If this is the case, the name must be declared EXTERNAL in
!             the user's calling program.  Given a system of N
!             differential equations, it is meaningful to speak about
!             the partial derivative of the I-th right hand side with
!             respect to the J-th dependent variable.  In general there
!             are N*N such quantities.  Often however the equations can
!             be ordered so that the I-th differential equation only
!             involves dependent variables with index near I, e.g., I+1,
!             I-2.  Such a system is called banded.  If, for all I, the
!             I-th equation depends on at most the variables
!               Y(I-ML), Y(I-ML+1), ... , Y(I), Y(I+1), ... , Y(I+MU)
!             then we call ML+MU+1 the bandwidth of the system.  In a
!             banded system many of the partial derivatives above are
!             automatically zero.  For the cases MITER = 1, 2, 4, and 5,
!             some of these partials are needed.  For the cases
!             MITER = 2 and 5 the necessary derivatives are
!             approximated numerically by CDRIV3, and we only ask the
!             user to tell CDRIV3 the value of ML and MU if the system
!             is banded.  For the cases MITER = 1 and 4 the user must
!             derive these partials algebraically and encode them in
!             subroutine JACOBN.  By computing these derivatives the
!             user can often save 20-30 per cent of the computing time.
!             Usually, however, the accuracy is not much affected and
!             most users will probably forego this option.  The optional
!             user-written subroutine JACOBN has the form:
!                   SUBROUTINE JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
!                   COMPLEX Y(*), DFDY(MATDIM,*)
!                     .
!                     .
!                     Calculate values of DFDY
!                     .
!                     .
!                   END (Sample)
!             Here Y is a vector of length at least N.  The actual
!             length of Y is determined by the user's declaration in the
!             program which calls CDRIV3.  Thus the dimensioning of Y in
!             JACOBN, while required by FORTRAN convention, does not
!             actually allocate any storage.  When this subroutine is
!             called, the first N components of Y are intermediate
!             approximations to the solution components.  The user
!             should not alter these values.  If the system is not
!             banded (MITER=1), the partials of the I-th equation with
!             respect to the J-th dependent function are to be stored in
!             DFDY(I,J).  Thus partials of the I-th equation are stored
!             in the I-th row of DFDY.  If the system is banded
!             (MITER=4), then the partials of the I-th equation with
!             respect to Y(J) are to be stored in DFDY(K,J), where
!             K=I-J+MU+1 .  Normally a return from JACOBN passes control
!             back to CDRIV3.  However, if the user would like to abort
!             the calculation, i.e., return control to the program which
!             calls CDRIV3, he should set N to zero.  CDRIV3 will signal
!             this by returning a value of NSTATE equal to +8(-8).
!             Altering the value of N in JACOBN has no effect on the
!             value of N in the call sequence of CDRIV3.
!
!    FA     = A subroutine supplied by the user if IMPL is not zero, and
!             MITER is not 3.  If so, the name must be declared EXTERNAL
!             in the user's calling program.  This subroutine computes
!             the array A, where A*dY(I)/dT = F(Y(I),T).
!             There are three cases:
!
!               IMPL=1.
!               Subroutine FA is of the form:
!                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
!                   COMPLEX Y(*), A(MATDIM,*)
!                     .
!                     .
!                     Calculate ALL values of A
!                     .
!                     .
!                   END (Sample)
!               In this case A is assumed to be a nonsingular matrix,
!               with the same structure as DFDY (see JACOBN description
!               above).  Programming considerations prevent complete
!               generality.  If MITER is 1 or 2, A is assumed to be full
!               and the user must compute and store all values of
!               A(I,J), I,J=1, ... ,N.  If MITER is 4 or 5, A is assumed
!               to be banded with lower and upper half bandwidth ML and
!               MU.  The left hand side of the I-th equation is a linear
!               combination of dY(I-ML)/dT, dY(I-ML+1)/dT, ... ,
!               dY(I)/dT, ... , dY(I+MU-1)/dT, dY(I+MU)/dT.  Thus in the
!               I-th equation, the coefficient of dY(J)/dT is to be
!               stored in A(K,J), where K=I-J+MU+1.
!               NOTE: The array A will be altered between calls to FA.
!
!               IMPL=2.
!               Subroutine FA is of the form:
!                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
!                   COMPLEX Y(*), A(*)
!                     .
!                     .
!                     Calculate non-zero values of A(1),...,A(NDE)
!                     .
!                     .
!                   END (Sample)
!               In this case it is assumed that the system is ordered by
!               the user so that the differential equations appear
!               first, and the algebraic equations appear last.  The
!               algebraic equations must be written in the form:
!               0 = F(Y(I),T).  When using this option it is up to the
!               user to provide initial values for the Y(I) that satisfy
!               the algebraic equations as well as possible.  It is
!               further assumed that A is a vector of length NDE.  All
!               of the components of A, which may depend on T, Y(I),
!               etc., must be set by the user to non-zero values.
!
!               IMPL=3.
!               Subroutine FA is of the form:
!                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
!                   COMPLEX Y(*), A(MATDIM,*)
!                     .
!                     .
!                     Calculate ALL values of A
!                     .
!                     .
!                   END (Sample)
!               In this case A is assumed to be a nonsingular NDE by NDE
!               matrix with the same structure as DFDY (see JACOBN
!               description above).  Programming considerations prevent
!               complete generality.  If MITER is 1 or 2, A is assumed
!               to be full and the user must compute and store all
!               values of A(I,J), I,J=1, ... ,NDE.  If MITER is 4 or 5,
!               A is assumed to be banded with lower and upper half
!               bandwidths ML and MU.  The left hand side of the I-th
!               equation is a linear combination of dY(I-ML)/dT,
!               dY(I-ML+1)/dT, ... , dY(I)/dT, ... , dY(I+MU-1)/dT,
!               dY(I+MU)/dT.  Thus in the I-th equation, the coefficient
!               of dY(J)/dT is to be stored in A(K,J), where K=I-J+MU+1.
!               It is assumed that the system is ordered by the user so
!               that the differential equations appear first, and the
!               algebraic equations appear last.  The algebraic
!               equations must be written in the form 0 = F(Y(I),T).
!               When using this option it is up to the user to provide
!               initial values for the Y(I) that satisfy the algebraic
!               equations as well as possible.
!               NOTE: For IMPL = 3, the array A will be altered between
!               calls to FA.
!             Here Y is a vector of length at least N.  The actual
!             length of Y is determined by the user's declaration in the
!             program which calls CDRIV3.  Thus the dimensioning of Y in
!             FA, while required by FORTRAN convention, does not
!             actually allocate any storage.  When this subroutine is
!             called, the first N components of Y are intermediate
!             approximations to the solution components.  The user
!             should not alter these values.  FA is always called
!             immediately after calling F, with the same values of T
!             and Y.  Normally a return from FA passes control back to
!             CDRIV3.  However, if the user would like to abort the
!             calculation, i.e., return control to the program which
!             calls CDRIV3, he should set N to zero.  CDRIV3 will signal
!             this by returning a value of NSTATE equal to +9(-9).
!             Altering the value of N in FA has no effect on the value
!             of N in the call sequence of CDRIV3.
!
!    NDE    = (Input) The number of differential equations.  This is
!             required only for IMPL = 2 or 3, with NDE  <  N.
!
!    MXSTEP = (Input) The maximum number of internal steps allowed on
!             one call to CDRIV3.
!
!    G      = A real FORTRAN function supplied by the user
!             if NROOT is not 0.  In this case, the name must be
!             declared EXTERNAL in the user's calling program.  G is
!             repeatedly called with different values of IROOT to obtain
!             the value of each of the NROOT equations for which a root
!             is desired.  G is of the form:
!                   REAL FUNCTION G (N, T, Y, IROOT)
!                   COMPLEX Y(*)
!                   go to (10, ...), IROOT
!              10   G = ...
!                     .
!                     .
!                   END (Sample)
!             Here, Y is a vector of length at least N, whose first N
!             components are the solution components at the point T.
!             The user should not alter these values.  The actual length
!             of Y is determined by the user's declaration in the
!             program which calls CDRIV3.  Thus the dimensioning of Y in
!             G, while required by FORTRAN convention, does not actually
!             allocate any storage.  Normally a return from G passes
!             control back to  CDRIV3.  However, if the user would like
!             to abort the calculation, i.e., return control to the
!             program which calls CDRIV3, he should set N to zero.
!             CDRIV3 will signal this by returning a value of NSTATE
!             equal to +7(-7).  In this case, the index of the equation
!             being evaluated is stored in the sixth element of IWORK.
!             Altering the value of N in G has no effect on the value of
!             N in the call sequence of CDRIV3.
!
!    USERS  = A subroutine supplied by the user, if MITER is 3.
!             If this is the case, the name must be declared EXTERNAL in
!             the user's calling program.  The routine USERS is called
!             by CDRIV3 when certain linear systems must be solved.  The
!             user may choose any method to form, store and solve these
!             systems in order to obtain the solution result that is
!             returned to CDRIV3.  In particular, this allows sparse
!             matrix methods to be used.  The call sequence for this
!             routine is:
!
!                SUBROUTINE USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL,
!               8                  IMPL, N, NDE, IFLAG)
!                COMPLEX Y(*), YH(*), YWT(*), SAVE1(*), SAVE2(*)
!                REAL T, H, EL
!
!             The input variable IFLAG indicates what action is to be
!             taken.  Subroutine USERS should perform the following
!             operations, depending on the value of IFLAG and IMPL.
!
!               IFLAG = 0
!                 IMPL = 0.  USERS is not called.
!                 IMPL = 1, 2 or 3.  Solve the system A*X = SAVE2,
!                   returning the result in SAVE2.  The array SAVE1 can
!                   be used as a work array.  For IMPL = 1, there are N
!                   components to the system, and for IMPL = 2 or 3,
!                   there are NDE components to the system.
!
!               IFLAG = 1
!                 IMPL = 0.  Compute, decompose and store the matrix
!                   (I - H*EL*J), where I is the identity matrix and J
!                   is the Jacobian matrix of the right hand side.  The
!                   array SAVE1 can be used as a work array.
!                 IMPL = 1, 2 or 3. Compute, decompose and store the
!                   matrix (A - H*EL*J).  The array SAVE1 can be used as
!                   a work array.
!
!               IFLAG = 2
!                 IMPL = 0.   Solve the system
!                     (I - H*EL*J)*X = H*SAVE2 - YH - SAVE1,
!                   returning the result in SAVE2.
!                 IMPL = 1, 2 or 3.  Solve the system
!                   (A - H*EL*J)*X = H*SAVE2 - A*(YH + SAVE1)
!                   returning the result in SAVE2.
!                 The array SAVE1 should not be altered.
!             If IFLAG is 0 and IMPL is 1 or 2 and the matrix A is
!             singular, or if IFLAG is 1 and one of the matrices
!             (I - H*EL*J), (A - H*EL*J) is singular, the INTEGER
!             variable IFLAG is to be set to -1 before RETURNing.
!             Normally a return from USERS passes control back to
!             CDRIV3.  However, if the user would like to abort the
!             calculation, i.e., return control to the program which
!             calls CDRIV3, he should set N to zero.  CDRIV3 will signal
!             this by returning a value of NSTATE equal to +10(-10).
!             Altering the value of N in USERS has no effect on the
!             value of N in the call sequence of CDRIV3.
!
!    IERFLG = An error flag.  The error number associated with a
!             diagnostic message (see Section III-A below) is the same
!             as the corresponding value of IERFLG.  The meaning of
!             IERFLG:
!               0  The routine completed successfully. (No message is
!                  issued.)
!               3  (Warning) The number of steps required to reach TOUT
!                  exceeds MXSTEP.
!               4  (Warning) The value of EPS is too small.
!              11  (Warning) For NTASK = 2 or 3, T is beyond TOUT.
!                  The solution was obtained by interpolation.
!              15  (Warning) The integration step size is below the
!                  roundoff level of T.  (The program issues this
!                  message as a warning but does not return control to
!                  the user.)
!              22  (Recoverable) N is not positive.
!              23  (Recoverable) MINT is less than 1 or greater than 3 .
!              24  (Recoverable) MITER is less than 0 or greater than
!                  5 .
!              25  (Recoverable) IMPL is less than 0 or greater than 3 .
!              26  (Recoverable) The value of NSTATE is less than 1 or
!                  greater than 12 .
!              27  (Recoverable) EPS is less than zero.
!              28  (Recoverable) MXORD is not positive.
!              29  (Recoverable) For MINT = 3, either MITER = 0 or 3, or
!                  IMPL = 0 .
!              30  (Recoverable) For MITER = 0, IMPL is not 0 .
!              31  (Recoverable) For MINT = 1, IMPL is 2 or 3 .
!              32  (Recoverable) Insufficient storage has been allocated
!                  for the WORK array.
!              33  (Recoverable) Insufficient storage has been allocated
!                  for the IWORK array.
!              41  (Recoverable) The integration step size has gone
!                  to zero.
!              42  (Recoverable) The integration step size has been
!                  reduced about 50 times without advancing the
!                  solution.  The problem setup may not be correct.
!              43  (Recoverable)  For IMPL greater than 0, the matrix A
!                  is singular.
!             999  (Fatal) The value of NSTATE is 12 .
!
!  III.  OTHER COMMUNICATION TO THE USER  ..............................
!
!    A. The solver communicates to the user through the parameters
!       above.  In addition it writes diagnostic messages through the
!       standard error handling program XERMSG.  A complete description
!       of XERMSG is given in "Guide to the SLATEC Common Mathematical
!       Library" by Kirby W. Fong et al..  At installations which do not
!       have this error handling package the short but serviceable
!       routine, XERMSG, available with this package, can be used.  That
!       program uses the file named OUTPUT to transmit messages.
!
!    B. The first three elements of WORK and the first five elements of
!       IWORK will contain the following statistical data:
!         AVGH     The average step size used.
!         HUSED    The step size last used (successfully).
!         AVGORD   The average order used.
!         IMXERR   The index of the element of the solution vector that
!                  contributed most to the last error test.
!         NQUSED   The order last used (successfully).
!         NSTEP    The number of steps taken since last initialization.
!         NFE      The number of evaluations of the right hand side.
!         NJE      The number of evaluations of the Jacobian matrix.
!
!  IV.  REMARKS  .......................................................
!
!    A. Other routines used:
!         CDNTP, CDZRO, CDSTP, CDNTL, CDPST, CDCOR, CDCST,
!         CDPSC, and CDSCL;
!         CGEFA, CGESL, CGBFA, CGBSL, and SCNRM2 (from LINPACK)
!         R1MACH (from the Bell Laboratories Machine Constants Package)
!         XERMSG (from the SLATEC Common Math Library)
!       The last seven routines above, not having been written by the
!       present authors, are not explicitly part of this package.
!
!    B. On any return from CDRIV3 all information necessary to continue
!       the calculation is contained in the call sequence parameters,
!       including the work arrays.  Thus it is possible to suspend one
!       problem, integrate another, and then return to the first.
!
!    C. If this package is to be used in an overlay situation, the user
!       must declare in the primary overlay the variables in the call
!       sequence to CDRIV3.
!
!    D. Changing parameters during an integration.
!       The value of NROOT, EPS, EWT, IERROR, MINT, MITER, or HMAX may
!       be altered by the user between calls to CDRIV3.  For example, if
!       too much accuracy has been requested (the program returns with
!       NSTATE = 4 and an increased value of EPS) the user may wish to
!       increase EPS further.  In general, prudence is necessary when
!       making changes in parameters since such changes are not
!       implemented until the next integration step, which is not
!       necessarily the next call to CDRIV3.  This can happen if the
!       program has already integrated to a point which is beyond the
!       new point TOUT.
!
!    E. As the price for complete control of matrix algebra, the CDRIV3
!       USERS option puts all responsibility for Jacobian matrix
!       evaluation on the user.  It is often useful to approximate
!       numerically all or part of the Jacobian matrix.  However this
!       must be done carefully.  The FORTRAN sequence below illustrates
!       the method we recommend.  It can be inserted directly into
!       subroutine USERS to approximate Jacobian elements in rows I1
!       to I2 and columns J1 to J2.
!             COMPLEX DFDY(N,N), R, SAVE1(N), SAVE2(N), Y(N), YJ, YWT(N)
!             REAL EPSJ, H, R1MACH, T, UROUND
!             UROUND = R1MACH(4)
!             EPSJ = SQRT(UROUND)
!             DO 30 J = J1,J2
!               if (ABS(Y(J))  >  ABS(YWT(J))) THEN
!                 R = EPSJ*Y(J)
!               ELSE
!                 R = EPSJ*YWT(J)
!               end if
!               if (R  ==  0.E0) R = YWT(J)
!               YJ = Y(J)
!               Y(J) = Y(J) + R
!               call F (N, T, Y, SAVE1)
!               if (N  ==  0) RETURN
!               Y(J) = YJ
!               DO 20 I = I1,I2
!        20       DFDY(I,J) = (SAVE1(I) - SAVE2(I))/R
!        30     CONTINUE
!       Many problems give rise to structured sparse Jacobians, e.g.,
!       block banded.  It is possible to approximate them with fewer
!       function evaluations than the above procedure uses; see Curtis,
!       Powell and Reid, J. Inst. Maths Applics, (1974), Vol. 13,
!       pp. 117-119.
!
!    F. When any of the routines JACOBN, FA, G, or USERS, is not
!       required, difficulties associated with unsatisfied externals can
!       be avoided by using the name of the routine which calculates the
!       right hand side of the differential equations in place of the
!       corresponding name in the call sequence of CDRIV3.
!
!***REFERENCES  C. W. Gear, Numerical Initial Value Problems in
!                 Ordinary Differential Equations, Prentice-Hall, 1971.
!***ROUTINES CALLED  CDNTP, CDSTP, CDZRO, CGBFA, CGBSL, CGEFA, CGESL,
!                    R1MACH, SCNRM2, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   900329  Initial submission to SLATEC.
!***END PROLOGUE  CDRIV3
  EXTERNAL F, JACOBN, FA, G, USERS
  COMPLEX WORK(*), Y(*)
  REAL AE, AVGH, AVGORD, BIG, EL(13,12), EPS, EWT(*), &
       G, GLAST, GNOW, H, HMAX, HOLD, HSIGN, HUSED, NROUND, RC, RE, &
       RMAX, R1MACH, SCNRM2, SIZE, SUM, T, TLAST, TOUT, TQ(3,12), &
       TREND, TROOT, UROUND
  INTEGER I, IA, IAVGH, IAVGRD, ICNVRG, IDFDY, IEL, IERFLG, IERROR, &
          IFAC, IFLAG, IGNOW, IH, IHMAX, IHOLD, IHSIGN, IHUSED, &
          IJROOT, IJSTPL, IJTASK, IMNT, IMNTLD, IMPL, IMTR, IMTRLD, &
          IMTRSV, IMXERR, IMXORD, IMXRDS, INDMXR, INDPRT, INDPVT, &
          INDTRT, INFE, INFO, INJE, INQ, INQUSE, INROOT, INRTLD, &
          INSTEP, INWAIT, IRC, IRMAX, IROOT, IMACH1, IMACH4, ISAVE1, &
          ISAVE2, IT, ITOUT, ITQ, ITREND, ITROOT, IWORK(*), IYH, &
          IYWT, J, JSTATE, JTROOT, LENCHK, LENIW, LENW, LIWCHK, &
          MATDIM, MAXORD, MINT, MITER, ML, MU, MXORD, MXSTEP, N, &
          NDE, NDECOM, NPAR, NROOT, NSTATE, NSTEPL, NTASK
  LOGICAL CONVRG
  CHARACTER INTGR1*8, INTGR2*8, RL1*16, RL2*16
  PARAMETER(NROUND = 20.E0)
  PARAMETER(IAVGH = 1, IHUSED = 2, IAVGRD = 3, &
            IEL = 4, IH = 160, IHMAX = 161, IHOLD = 162, &
            IHSIGN = 163, IRC = 164, IRMAX = 165, IT = 166, &
            ITOUT = 167, ITQ = 168, ITREND = 204, IMACH1 = 205, &
            IMACH4 = 206, IYH = 251, &
            INDMXR = 1, INQUSE = 2, INSTEP = 3, INFE = 4, INJE = 5, &
            INROOT = 6, ICNVRG = 7, IJROOT = 8, IJTASK = 9, &
            IMNTLD = 10, IMTRLD = 11, INQ = 12, INRTLD = 13, &
            INDTRT = 14, INWAIT = 15, IMNT = 16, IMTRSV = 17, &
            IMTR = 18, IMXRDS = 19, IMXORD = 20, INDPRT = 21, &
            IJSTPL = 22, INDPVT = 51)
!***FIRST EXECUTABLE STATEMENT  CDRIV3
  if (NSTATE  ==  12) THEN
    IERFLG = 999
    call XERMSG('SLATEC', 'CDRIV3', &
    'Illegal input.  The value of NSTATE is 12 .', IERFLG, 2)
    return
  ELSE if (NSTATE  <  1 .OR. NSTATE  >  12) THEN
    WRITE(INTGR1, '(I8)') NSTATE
    IERFLG = 26
    call XERMSG('SLATEC', 'CDRIV3', &
    'Illegal input.  Improper value for NSTATE(= '//INTGR1//').', &
    IERFLG, 1)
    NSTATE = 12
    return
  end if
  NPAR = N
  if (EPS  <  0.E0) THEN
    WRITE(RL1, '(E16.8)') EPS
    IERFLG = 27
    call XERMSG('SLATEC', 'CDRIV3', &
    'Illegal input.  EPS, '//RL1//', is negative.', IERFLG, 1)
    NSTATE = 12
    return
  end if
  if (N  <=  0) THEN
    WRITE(INTGR1, '(I8)') N
    IERFLG = 22
    call XERMSG('SLATEC', 'CDRIV3', &
    'Illegal input.  Number of equations, '//INTGR1// &
    ', is not positive.', IERFLG, 1)
    NSTATE = 12
    return
  end if
  if (MXORD  <=  0) THEN
    WRITE(INTGR1, '(I8)') MXORD
    IERFLG = 28
    call XERMSG('SLATEC', 'CDRIV3', &
    'Illegal input.  Maximum order, '//INTGR1// &
    ', is not positive.', IERFLG, 1)
    NSTATE = 12
    return
  end if
  if (MINT  <  1 .OR. MINT  >  3) THEN
    WRITE(INTGR1, '(I8)') MINT
    IERFLG = 23
    call XERMSG('SLATEC', 'CDRIV3', &
    'Illegal input.  Improper value for the integration method '// &
    'flag, '//INTGR1//' .', IERFLG, 1)
    NSTATE = 12
    return
  ELSE if (MITER  <  0 .OR. MITER  >  5) THEN
    WRITE(INTGR1, '(I8)') MITER
    IERFLG = 24
    call XERMSG('SLATEC', 'CDRIV3', &
    'Illegal input.  Improper value for MITER(= '//INTGR1//').', &
    IERFLG, 1)
    NSTATE = 12
    return
  ELSE if (IMPL  <  0 .OR. IMPL  >  3) THEN
    WRITE(INTGR1, '(I8)') IMPL
    IERFLG = 25
    call XERMSG('SLATEC', 'CDRIV3', &
    'Illegal input.  Improper value for IMPL(= '//INTGR1//').', &
    IERFLG, 1)
    NSTATE = 12
    return
  ELSE if (MINT  ==  3 .AND. &
    (MITER  ==  0 .OR. MITER  ==  3 .OR. IMPL  /=  0)) THEN
    WRITE(INTGR1, '(I8)') MITER
    WRITE(INTGR2, '(I8)') IMPL
    IERFLG = 29
    call XERMSG('SLATEC', 'CDRIV3', &
    'Illegal input.  For MINT = 3, the value of MITER, '//INTGR1// &
    ', and/or IMPL, '//INTGR2//', is not allowed.', IERFLG, 1)
    NSTATE = 12
    return
  ELSE if ((IMPL  >=  1 .AND. IMPL  <=  3) .AND. MITER  ==  0) THEN
    WRITE(INTGR1, '(I8)') IMPL
    IERFLG = 30
    call XERMSG('SLATEC', 'CDRIV3', &
    'Illegal input.  For MITER = 0, the value of IMPL, '//INTGR1// &
    ', is not allowed.', IERFLG, 1)
    NSTATE = 12
    return
  ELSE if ((IMPL  ==  2 .OR. IMPL  ==  3) .AND. MINT  ==  1) THEN
    WRITE(INTGR1, '(I8)') IMPL
    IERFLG = 31
    call XERMSG('SLATEC', 'CDRIV3', &
    'Illegal input.  For MINT = 1, the value of IMPL, '//INTGR1// &
    ', is not allowed.', IERFLG, 1)
    NSTATE = 12
    return
  end if
  if (MITER  ==  0 .OR. MITER  ==  3) THEN
    LIWCHK = INDPVT - 1
  ELSE if (MITER  ==  1 .OR. MITER  ==  2 .OR. MITER  ==  4 .OR. &
    MITER  ==  5) THEN
    LIWCHK = INDPVT + N - 1
  end if
  if (LENIW  <  LIWCHK) THEN
    WRITE(INTGR1, '(I8)') LIWCHK
    IERFLG = 33
    call XERMSG('SLATEC', 'CDRIV3', &
    'Illegal input.  Insufficient storage allocated for the '// &
    'IWORK array.  Based on the value of the input parameters '// &
    'involved, the required storage is '//INTGR1//' .', IERFLG, 1)
    NSTATE = 12
    return
  end if
!                                                Allocate the WORK array
!                                         IYH is the index of YH in WORK
  if (MINT  ==  1 .OR. MINT  ==  3) THEN
    MAXORD = MIN(MXORD, 12)
  ELSE if (MINT  ==  2) THEN
    MAXORD = MIN(MXORD, 5)
  end if
  IDFDY = IYH + (MAXORD + 1)*N
!                                             IDFDY is the index of DFDY
!
  if (MITER  ==  0 .OR. MITER  ==  3) THEN
    IYWT = IDFDY
  ELSE if (MITER  ==  1 .OR. MITER  ==  2) THEN
    IYWT = IDFDY + N*N
  ELSE if (MITER  ==  4 .OR. MITER  ==  5) THEN
    IYWT = IDFDY + (2*ML + MU + 1)*N
  end if
!                                               IYWT is the index of YWT
  ISAVE1 = IYWT + N
!                                           ISAVE1 is the index of SAVE1
  ISAVE2 = ISAVE1 + N
!                                           ISAVE2 is the index of SAVE2
  IGNOW = ISAVE2 + N
!                                             IGNOW is the index of GNOW
  ITROOT = IGNOW + NROOT
!                                           ITROOT is the index of TROOT
  IFAC = ITROOT + NROOT
!                                               IFAC is the index of FAC
  if (MITER  ==  2 .OR. MITER  ==  5 .OR. MINT  ==  3) THEN
    IA = IFAC + N
  ELSE
    IA = IFAC
  end if
!                                                   IA is the index of A
  if (IMPL  ==  0 .OR. MITER  ==  3) THEN
    LENCHK = IA - 1
  ELSE if (IMPL  ==  1 .AND. (MITER  ==  1 .OR. MITER  ==  2)) THEN
    LENCHK = IA - 1 + N*N
  ELSE if (IMPL  ==  1 .AND. (MITER  ==  4 .OR. MITER  ==  5)) THEN
    LENCHK = IA - 1 + (2*ML + MU + 1)*N
  ELSE if (IMPL  ==  2 .AND. MITER  /=  3) THEN
    LENCHK = IA - 1 + N
  ELSE if (IMPL  ==  3 .AND. (MITER  ==  1 .OR. MITER  ==  2)) THEN
    LENCHK = IA - 1 + N*NDE
  ELSE if (IMPL  ==  3 .AND. (MITER  ==  4 .OR. MITER  ==  5)) THEN
    LENCHK = IA - 1 + (2*ML + MU + 1)*NDE
  end if
  if (LENW  <  LENCHK) THEN
    WRITE(INTGR1, '(I8)') LENCHK
    IERFLG = 32
    call XERMSG('SLATEC', 'CDRIV3', &
    'Illegal input.  Insufficient storage allocated for the '// &
    'WORK array.  Based on the value of the input parameters '// &
    'involved, the required storage is '//INTGR1//' .', IERFLG, 1)
    NSTATE = 12
    return
  end if
  if (MITER  ==  0 .OR. MITER  ==  3) THEN
    MATDIM = 1
  ELSE if (MITER  ==  1 .OR. MITER  ==  2) THEN
    MATDIM = N
  ELSE if (MITER  ==  4 .OR. MITER  ==  5) THEN
    MATDIM = 2*ML + MU + 1
  end if
  if (IMPL  ==  0 .OR. IMPL  ==  1) THEN
    NDECOM = N
  ELSE if (IMPL  ==  2 .OR. IMPL  ==  3) THEN
    NDECOM = NDE
  end if
  if (NSTATE  ==  1) THEN
!                                                  Initialize parameters
    if (MINT  ==  1 .OR. MINT  ==  3) THEN
      IWORK(IMXORD) = MIN(MXORD, 12)
    ELSE if (MINT  ==  2) THEN
      IWORK(IMXORD) = MIN(MXORD, 5)
    end if
    IWORK(IMXRDS) = MXORD
    if (MINT  ==  1 .OR. MINT  ==  2) THEN
      IWORK(IMNT) = MINT
      IWORK(IMTR) = MITER
      IWORK(IMNTLD) = MINT
      IWORK(IMTRLD) = MITER
    ELSE if (MINT  ==  3) THEN
      IWORK(IMNT) = 1
      IWORK(IMTR) = 0
      IWORK(IMNTLD) = IWORK(IMNT)
      IWORK(IMTRLD) = IWORK(IMTR)
      IWORK(IMTRSV) = MITER
    end if
    WORK(IHMAX) = HMAX
    UROUND = R1MACH (4)
    WORK(IMACH4) = UROUND
    WORK(IMACH1) = R1MACH (1)
    if (NROOT  /=  0) THEN
      RE = UROUND
      AE = WORK(IMACH1)
    end if
    H = (TOUT - T)*(1.E0 - 4.E0*UROUND)
    H = SIGN(MIN(ABS(H), HMAX), H)
    WORK(IH) = H
    HSIGN = SIGN(1.E0, H)
    WORK(IHSIGN) = HSIGN
    IWORK(IJTASK) = 0
    AVGH = 0.E0
    AVGORD = 0.E0
    WORK(IAVGH) = 0.E0
    WORK(IHUSED) = 0.E0
    WORK(IAVGRD) = 0.E0
    IWORK(INDMXR) = 0
    IWORK(INQUSE) = 0
    IWORK(INSTEP) = 0
    IWORK(IJSTPL) = 0
    IWORK(INFE) = 0
    IWORK(INJE) = 0
    IWORK(INROOT) = 0
    WORK(IT) = T
    IWORK(ICNVRG) = 0
    IWORK(INDPRT) = 0
!                                                 Set initial conditions
    DO 30 I = 1,N
 30       WORK(I+IYH-1) = Y(I)
    if (T  ==  TOUT) RETURN
    go to 180
  ELSE
    UROUND = WORK(IMACH4)
    if (NROOT  /=  0) THEN
      RE = UROUND
      AE = WORK(IMACH1)
    end if
  end if
!                                             On a continuation, check
!                                             that output points have
!                                             been or will be overtaken.
  if (IWORK(ICNVRG)  ==  1) THEN
    CONVRG = .TRUE.
  ELSE
    CONVRG = .FALSE.
  end if
  AVGH = WORK(IAVGH)
  AVGORD = WORK(IAVGRD)
  HOLD = WORK(IHOLD)
  RC = WORK(IRC)
  RMAX = WORK(IRMAX)
  TREND = WORK(ITREND)
  DO 35 J = 1,12
    DO 35 I = 1,13
 35       EL(I,J) = WORK(I+IEL+(J-1)*13-1)
  DO 40 J = 1,12
    DO 40 I = 1,3
 40       TQ(I,J) = WORK(I+ITQ+(J-1)*3-1)
  T = WORK(IT)
  H = WORK(IH)
  HSIGN = WORK(IHSIGN)
  if (IWORK(IJTASK)  ==  0) go to 180
!
!                                   IWORK(IJROOT) flags unreported
!                                   roots, and is set to the value of
!                                   NTASK when a root was last selected.
!                                   It is set to zero when all roots
!                                   have been reported.  IWORK(INROOT)
!                                   contains the index and WORK(ITOUT)
!                                   contains the value of the root last
!                                   selected to be reported.
!                                   IWORK(INRTLD) contains the value of
!                                   NROOT and IWORK(INDTRT) contains
!                                   the value of ITROOT when the array
!                                   of roots was last calculated.
  if (NROOT  /=  0) THEN
    if (IWORK(IJROOT)  >  0) THEN
!                                      TOUT has just been reported.
!                                      If TROOT  <=  TOUT, report TROOT.
      if (NSTATE  /=  5) THEN
        if (TOUT*HSIGN  >=  REAL(WORK(ITOUT))*HSIGN) THEN
          TROOT = WORK(ITOUT)
          call CDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
          T = TROOT
          NSTATE = 5
          IERFLG = 0
          go to 580
        end if
!                                         A root has just been reported.
!                                         Select the next root.
      ELSE
        TROOT = T
        IROOT = 0
        DO 50 I = 1,IWORK(INRTLD)
          JTROOT = I + IWORK(INDTRT) - 1
          if (REAL(WORK(JTROOT))*HSIGN  <=  TROOT*HSIGN) THEN
!
!                                              Check for multiple roots.
!
            if (WORK(JTROOT)  ==  WORK(ITOUT) .AND. &
            I  >  IWORK(INROOT)) THEN
              IROOT = I
              TROOT = WORK(JTROOT)
              go to 60
            end if
            if (REAL(WORK(JTROOT))*HSIGN  >  &
            REAL(WORK(ITOUT))*HSIGN) THEN
              IROOT = I
              TROOT = WORK(JTROOT)
            end if
          end if
 50           CONTINUE
 60         IWORK(INROOT) = IROOT
        WORK(ITOUT) = TROOT
        IWORK(IJROOT) = NTASK
        if (NTASK  ==  1) THEN
          if (IROOT  ==  0) THEN
            IWORK(IJROOT) = 0
          ELSE
            if (TOUT*HSIGN  >=  TROOT*HSIGN) THEN
              call CDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH), &
                          Y)
              NSTATE = 5
              T = TROOT
              IERFLG = 0
              go to 580
            end if
          end if
        ELSE if (NTASK  ==  2 .OR. NTASK  ==  3) THEN
!
!                                     If there are no more roots, or the
!                                     user has altered TOUT to be less
!                                     than a root, set IJROOT to zero.
!
          if (IROOT  ==  0 .OR. (TOUT*HSIGN  <  TROOT*HSIGN)) THEN
            IWORK(IJROOT) = 0
          ELSE
            call CDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH), &
                        Y)
            NSTATE = 5
            T = TROOT
            IERFLG = 0
            go to 580
          end if
        end if
      end if
    end if
  end if
!
  if (NTASK  ==  1) THEN
    NSTATE = 2
    if (T*HSIGN  >=  TOUT*HSIGN) THEN
      call CDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
      T = TOUT
      IERFLG = 0
      go to 580
    end if
  ELSE if (NTASK  ==  2) THEN
!                                                      Check if TOUT has
!                                                      been reset  <  T
    if (T*HSIGN  >  TOUT*HSIGN) THEN
      WRITE(RL1, '(E16.8)') T
      WRITE(RL2, '(E16.8)') TOUT
      IERFLG = 11
      call XERMSG('SLATEC', 'CDRIV3', &
      'While integrating exactly to TOUT, T, '//RL1// &
      ', was beyond TOUT, '//RL2//' .  Solution obtained by '// &
      'interpolation.', IERFLG, 0)
      NSTATE = 11
      call CDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
      T = TOUT
      go to 580
    end if
!                                   Determine if TOUT has been overtaken
!
    if (ABS(TOUT - T) <= NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
      T = TOUT
      NSTATE = 2
      IERFLG = 0
      go to 560
    end if
!                                             If there are no more roots
!                                             to report, report T.
    if (NSTATE  ==  5) THEN
      NSTATE = 2
      IERFLG = 0
      go to 560
    end if
    NSTATE = 2
!                                                       See if TOUT will
!                                                       be overtaken.
    if ((T + H)*HSIGN  >  TOUT*HSIGN) THEN
      H = TOUT - T
      if ((T + H)*HSIGN  >  TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
      WORK(IH) = H
      if (H  ==  0.E0) go to 670
      IWORK(IJTASK) = -1
    end if
  ELSE if (NTASK  ==  3) THEN
    NSTATE = 2
    if (T*HSIGN  >  TOUT*HSIGN) THEN
      WRITE(RL1, '(E16.8)') T
      WRITE(RL2, '(E16.8)') TOUT
      IERFLG = 11
      call XERMSG('SLATEC', 'CDRIV3', &
      'While integrating exactly to TOUT, T, '//RL1// &
      ', was beyond TOUT, '//RL2//' .  Solution obtained by '// &
      'interpolation.', IERFLG, 0)
      NSTATE = 11
      call CDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
      T = TOUT
      go to 580
    end if
    if (ABS(TOUT - T) <= NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
      T = TOUT
      IERFLG = 0
      go to 560
    end if
    if ((T + H)*HSIGN  >  TOUT*HSIGN) THEN
      H = TOUT - T
      if ((T + H)*HSIGN  >  TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
      WORK(IH) = H
      if (H  ==  0.E0) go to 670
      IWORK(IJTASK) = -1
    end if
  end if
!                         Implement changes in MINT, MITER, and/or HMAX.
!
  if ((MINT  /=  IWORK(IMNTLD) .OR. MITER  /=  IWORK(IMTRLD)) .AND. &
    MINT  /=  3 .AND. IWORK(IMNTLD)  /=  3) IWORK(IJTASK) = -1
  if (HMAX  /=  WORK(IHMAX)) THEN
    H = SIGN(MIN(ABS(H), HMAX), H)
    if (H  /=  WORK(IH)) THEN
      IWORK(IJTASK) = -1
      WORK(IH) = H
    end if
    WORK(IHMAX) = HMAX
  end if
!
 180  NSTEPL = IWORK(INSTEP)
  DO 190 I = 1,N
 190    Y(I) = WORK(I+IYH-1)
  if (NROOT  /=  0) THEN
    DO 200 I = 1,NROOT
      WORK(I+IGNOW-1) = G (NPAR, T, Y, I)
      if (NPAR  ==  0) THEN
        IWORK(INROOT) = I
        NSTATE = 7
        return
      end if
 200     CONTINUE
  end if
  if (IERROR  ==  1) THEN
    DO 230 I = 1,N
 230      WORK(I+IYWT-1) = 1.E0
    go to 410
  ELSE if (IERROR  ==  5) THEN
    DO 250 I = 1,N
 250      WORK(I+IYWT-1) = EWT(I)
    go to 410
  end if
!                                       Reset YWT array.  Looping point.
 260  if (IERROR  ==  2) THEN
    DO 280 I = 1,N
      if (Y(I)  ==  0.E0) go to 290
 280      WORK(I+IYWT-1) = Y(I)
    go to 410
 290    if (IWORK(IJTASK)  ==  0) THEN
      call F (NPAR, T, Y, WORK(ISAVE2))
      if (NPAR  ==  0) THEN
        NSTATE = 6
        return
      end if
      IWORK(INFE) = IWORK(INFE) + 1
      if (MITER  ==  3 .AND. IMPL  /=  0) THEN
        IFLAG = 0
        call USERS (Y, WORK(IYH), WORK(IYWT), WORK(ISAVE1), &
                    WORK(ISAVE2), T, H, REAL(WORK(IEL)), IMPL, NPAR, &
                    NDECOM, IFLAG)
        if (IFLAG  ==  -1) go to 690
        if (NPAR  ==  0) THEN
          NSTATE = 10
          return
        end if
      ELSE if (IMPL  ==  1) THEN
        if (MITER  ==  1 .OR. MITER  ==  2) THEN
          call FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
          if (NPAR  ==  0) THEN
            NSTATE = 9
            return
          end if
          call CGEFA (WORK(IA), MATDIM, N, IWORK(INDPVT), INFO)
          if (INFO  /=  0) go to 690
          call CGESL (WORK(IA), MATDIM, N, IWORK(INDPVT), &
                      WORK(ISAVE2), 0)
        ELSE if (MITER  ==  4 .OR. MITER  ==  5) THEN
          call FA (NPAR, T, Y, WORK(IA+ML), MATDIM, ML, MU, NDECOM)
          if (NPAR  ==  0) THEN
            NSTATE = 9
            return
          end if
          call CGBFA (WORK(IA), MATDIM, N, ML, MU, IWORK(INDPVT), &
                      INFO)
          if (INFO  /=  0) go to 690
          call CGBSL (WORK(IA), MATDIM, N, ML, MU, IWORK(INDPVT), &
                      WORK(ISAVE2), 0)
        end if
      ELSE if (IMPL  ==  2) THEN
        call FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
        if (NPAR  ==  0) THEN
          NSTATE = 9
          return
        end if
        DO 340 I = 1,NDECOM
          if (WORK(I+IA-1)  ==  0.E0) go to 690
 340          WORK(I+ISAVE2-1) = WORK(I+ISAVE2-1)/WORK(I+IA-1)
      ELSE if (IMPL  ==  3) THEN
        if (MITER  ==  1 .OR. MITER  ==  2) THEN
          call FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
          if (NPAR  ==  0) THEN
            NSTATE = 9
            return
          end if
          call CGEFA (WORK(IA), MATDIM, NDE, IWORK(INDPVT), INFO)
          if (INFO  /=  0) go to 690
          call CGESL (WORK(IA), MATDIM, NDE, IWORK(INDPVT), &
                      WORK(ISAVE2), 0)
        ELSE if (MITER  ==  4 .OR. MITER  ==  5) THEN
          call FA (NPAR, T, Y, WORK(IA+ML), MATDIM, ML, MU, NDECOM)
          if (NPAR  ==  0) THEN
            NSTATE = 9
            return
          end if
          call CGBFA (WORK(IA), MATDIM, NDE, ML, MU, IWORK(INDPVT), &
                      INFO)
          if (INFO  /=  0) go to 690
          call CGBSL (WORK(IA), MATDIM, NDE, ML, MU, IWORK(INDPVT), &
                      WORK(ISAVE2), 0)
        end if
      end if
    end if
    DO 360 J = I,N
      if (Y(J)  /=  0.E0) THEN
        WORK(J+IYWT-1) = Y(J)
      ELSE
        if (IWORK(IJTASK)  ==  0) THEN
          WORK(J+IYWT-1) = H*WORK(J+ISAVE2-1)
        ELSE
          WORK(J+IYWT-1) = WORK(J+IYH+N-1)
        end if
      end if
      if (WORK(J+IYWT-1)  ==  0.E0) WORK(J+IYWT-1) = UROUND
 360      CONTINUE
  ELSE if (IERROR  ==  3) THEN
    DO 380 I = 1,N
 380      WORK(I+IYWT-1) = MAX(EWT(1), ABS(Y(I)))
  ELSE if (IERROR  ==  4) THEN
    DO 400 I = 1,N
 400      WORK(I+IYWT-1) = MAX(EWT(I), ABS(Y(I)))
  end if
!
 410  DO 420 I = 1,N
 420    WORK(I+ISAVE2-1) = Y(I)/WORK(I+IYWT-1)
  SUM = SCNRM2(N, WORK(ISAVE2), 1)/SQRT(REAL(N))
  SUM = MAX(1.E0, SUM)
  if (EPS  <  SUM*UROUND) THEN
    EPS = SUM*UROUND*(1.E0 + 10.E0*UROUND)
    WRITE(RL1, '(E16.8)') T
    WRITE(RL2, '(E16.8)') EPS
    IERFLG = 4
    call XERMSG('SLATEC', 'CDRIV3', &
    'At T, '//RL1//', the requested accuracy, EPS, was not '// &
    'obtainable with the machine precision.  EPS has been '// &
    'increased to '//RL2//' .', IERFLG, 0)
    NSTATE = 4
    go to 560
  end if
  if (ABS(H)  >=  UROUND*ABS(T)) THEN
    IWORK(INDPRT) = 0
  ELSE if (IWORK(INDPRT)  ==  0) THEN
    WRITE(RL1, '(E16.8)') T
    WRITE(RL2, '(E16.8)') H
    IERFLG = 15
    call XERMSG('SLATEC', 'CDRIV3', &
    'At T, '//RL1//', the step size, '//RL2//', is smaller '// &
    'than the roundoff level of T.  This may occur if there is '// &
    'an abrupt change in the right hand side of the '// &
    'differential equations.', IERFLG, 0)
    IWORK(INDPRT) = 1
  end if
  if (NTASK /= 2) THEN
    if ((IWORK(INSTEP)-NSTEPL)  ==  MXSTEP) THEN
      WRITE(RL1, '(E16.8)') T
      WRITE(INTGR1, '(I8)') MXSTEP
      WRITE(RL2, '(E16.8)') TOUT
      IERFLG = 3
      call XERMSG('SLATEC', 'CDRIV3', &
      'At T, '//RL1//', '//INTGR1//' steps have been taken '// &
      'without reaching TOUT, '//RL2//' .', IERFLG, 0)
      NSTATE = 3
      go to 560
    end if
  end if
!
!     call CDSTP (EPS, F, FA, HMAX, IMPL, IERROR, JACOBN, MATDIM,
!    8            MAXORD, MINT, MITER, ML, MU, N, NDE, YWT, UROUND,
!    8            USERS,  AVGH, AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD,
!    8            NFE, NJE, NQUSED, NSTEP, T, Y, YH,  A, CONVRG,
!    8            DFDY, EL, FAC, HOLD, IPVT, JSTATE, JSTEPL, NQ, NWAIT,
!    8            RC, RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG, MTRSV,
!    8            MXRDSV)
!
  call CDSTP (EPS, F, FA, HMAX, IMPL, IERROR, JACOBN, MATDIM, &
              IWORK(IMXORD), IWORK(IMNT), IWORK(IMTR), ML, MU, NPAR, &
             NDECOM, WORK(IYWT), UROUND, USERS,  AVGH, AVGORD, H, &
             HUSED, IWORK(IJTASK), IWORK(IMNTLD), IWORK(IMTRLD), &
             IWORK(INFE), IWORK(INJE), IWORK(INQUSE), IWORK(INSTEP), &
             T, Y, WORK(IYH),  WORK(IA), CONVRG, WORK(IDFDY), EL, &
             WORK(IFAC), HOLD, IWORK(INDPVT), JSTATE, IWORK(IJSTPL), &
             IWORK(INQ), IWORK(INWAIT), RC, RMAX, WORK(ISAVE1), &
             WORK(ISAVE2), TQ, TREND, MINT, IWORK(IMTRSV), &
             IWORK(IMXRDS))
!
  WORK(IH) = H
  WORK(IT) = T
  go to (470, 670, 680, 690, 690, 660, 660, 660, 660, 660), JSTATE
 470  IWORK(IJTASK) = 1
!                                 Determine if a root has been overtaken
  if (NROOT  /=  0) THEN
    IROOT = 0
    DO 500 I = 1,NROOT
      GLAST = WORK(I+IGNOW-1)
      GNOW = G (NPAR, T, Y, I)
      if (NPAR  ==  0) THEN
        IWORK(INROOT) = I
        NSTATE = 7
        return
      end if
      WORK(I+IGNOW-1) = GNOW
      if (GLAST*GNOW  >  0.E0) THEN
        WORK(I+ITROOT-1) = T + H
      ELSE
        if (GNOW  ==  0.E0) THEN
          WORK(I+ITROOT-1) = T
          IROOT = I
        ELSE
          if (GLAST  ==  0.E0) THEN
            WORK(I+ITROOT-1) = T + H
          ELSE
            if (ABS(HUSED)  >=  UROUND*ABS(T)) THEN
              TLAST = T - HUSED
              IROOT = I
              TROOT = T
              call CDZRO (AE, G, H, NPAR, IWORK(INQ), IROOT, RE, T, &
                          WORK(IYH), UROUND,  TROOT, TLAST, &
                          GNOW, GLAST,  Y)
              DO 480 J = 1,N
 480                Y(J) = WORK(IYH+J-1)
              if (NPAR  ==  0) THEN
                IWORK(INROOT) = I
                NSTATE = 7
                return
              end if
              WORK(I+ITROOT-1) = TROOT
            ELSE
              WORK(I+ITROOT-1) = T
              IROOT = I
            end if
          end if
        end if
      end if
 500      CONTINUE
    if (IROOT  ==  0) THEN
      IWORK(IJROOT) = 0
!                                                  Select the first root
    ELSE
      IWORK(IJROOT) = NTASK
      IWORK(INRTLD) = NROOT
      IWORK(INDTRT) = ITROOT
      TROOT = T + H
      DO 510 I = 1,NROOT
        if (REAL(WORK(I+ITROOT-1))*HSIGN  <  TROOT*HSIGN) THEN
          TROOT = WORK(I+ITROOT-1)
          IROOT = I
        end if
 510        CONTINUE
      IWORK(INROOT) = IROOT
      WORK(ITOUT) = TROOT
      if (TROOT*HSIGN  <=  TOUT*HSIGN) THEN
        call CDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
        NSTATE = 5
        T = TROOT
        IERFLG = 0
        go to 580
      end if
    end if
  end if
!                               Test for NTASK condition to be satisfied
  NSTATE = 2
  if (NTASK  ==  1) THEN
    if (T*HSIGN  <  TOUT*HSIGN) go to 260
    call CDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
    T = TOUT
    IERFLG = 0
    go to 580
!                               TOUT is assumed to have been attained
!                               exactly if T is within twenty roundoff
!                               units of TOUT, relative to MAX(TOUT, T).
!
  ELSE if (NTASK  ==  2) THEN
    if (ABS(TOUT - T) <= NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
      T = TOUT
    ELSE
      if ((T + H)*HSIGN  >  TOUT*HSIGN) THEN
        H = TOUT - T
        if ((T + H)*HSIGN > TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
        WORK(IH) = H
        if (H  ==  0.E0) go to 670
        IWORK(IJTASK) = -1
      end if
    end if
  ELSE if (NTASK  ==  3) THEN
    if (ABS(TOUT - T) <= NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
      T = TOUT
    ELSE
      if ((T + H)*HSIGN  >  TOUT*HSIGN) THEN
        H = TOUT - T
        if ((T + H)*HSIGN > TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
        WORK(IH) = H
        if (H  ==  0.E0) go to 670
        IWORK(IJTASK) = -1
      end if
      go to 260
    end if
  end if
  IERFLG = 0
!                                      All returns are made through this
!                                      section.  IMXERR is determined.
 560  DO 570 I = 1,N
 570    Y(I) = WORK(I+IYH-1)
 580  if (CONVRG) THEN
    IWORK(ICNVRG) = 1
  ELSE
    IWORK(ICNVRG) = 0
  end if
  WORK(IAVGH) = AVGH
  WORK(IAVGRD) = AVGORD
  WORK(IHUSED) = HUSED
  WORK(IHOLD) = HOLD
  WORK(IRC) = RC
  WORK(IRMAX) = RMAX
  WORK(ITREND) = TREND
  DO 582 J = 1,12
    DO 582 I = 1,13
 582      WORK(I+IEL+(J-1)*13-1) = EL(I,J)
  DO 584 J = 1,12
    DO 584 I = 1,3
 584      WORK(I+ITQ+(J-1)*3-1) = TQ(I,J)
  if (IWORK(IJTASK)  ==  0) RETURN
  BIG = 0.E0
  IMXERR = 1
  DO  590 I = 1,N
!                                            SIZE = ABS(ERROR(I)/YWT(I))
    SIZE = ABS(WORK(I+ISAVE1-1)/WORK(I+IYWT-1))
    if (BIG  <  SIZE) THEN
      BIG = SIZE
      IMXERR = I
    end if
 590    CONTINUE
  IWORK(INDMXR) = IMXERR
  return
!
 660  NSTATE = JSTATE
  DO 662 I = 1,N
 662    Y(I) = WORK(I + IYH - 1)
  if (CONVRG) THEN
    IWORK(ICNVRG) = 1
  ELSE
    IWORK(ICNVRG) = 0
  end if
  WORK(IAVGH) = AVGH
  WORK(IAVGRD) = AVGORD
  WORK(IHUSED) = HUSED
  WORK(IHOLD) = HOLD
  WORK(IRC) = RC
  WORK(IRMAX) = RMAX
  WORK(ITREND) = TREND
  DO 664 J = 1,12
    DO 664 I = 1,13
 664      WORK(I+IEL+(J-1)*13-1) = EL(I,J)
  DO 666 J = 1,12
    DO 666 I = 1,3
 666      WORK(I+ITQ+(J-1)*3-1) = TQ(I,J)
  return
!                                        Fatal errors are processed here
!
 670  WRITE(RL1, '(E16.8)') T
  IERFLG = 41
  call XERMSG('SLATEC', 'CDRIV3', &
    'At T, '//RL1//', the attempted step size has gone to '// &
    'zero.  Often this occurs if the problem setup is incorrect.', &
    IERFLG, 1)
  NSTATE = 12
  return
!
 680  WRITE(RL1, '(E16.8)') T
  IERFLG = 42
  call XERMSG('SLATEC', 'CDRIV3', &
    'At T, '//RL1//', the step size has been reduced about 50 '// &
    'times without advancing the solution.  Often this occurs '// &
    'if the problem setup is incorrect.', IERFLG, 1)
  NSTATE = 12
  return
!
 690  WRITE(RL1, '(E16.8)') T
  IERFLG = 43
  call XERMSG('SLATEC', 'CDRIV3', &
    'At T, '//RL1//', while solving A*YDOT = F, A is singular.', &
    IERFLG, 1)
  NSTATE = 12
  return
end
