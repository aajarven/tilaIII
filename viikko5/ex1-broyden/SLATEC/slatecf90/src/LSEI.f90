subroutine LSEI (W, MDW, ME, MA, MG, N, PRGOPT, X, RNORME, RNORML, &
     MODE, WS, IP)
!
!! LSEI solves a linearly constrained least squares problem with ...
!            equality and inequality constraints, and optionally compute
!            a covariance matrix.
!
!***LIBRARY   SLATEC
!***CATEGORY  K1A2A, D9
!***TYPE      SINGLE PRECISION (LSEI-S, DLSEI-D)
!***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING,
!             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS,
!             QUADRATIC PROGRAMMING
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     Abstract
!
!     This subprogram solves a linearly constrained least squares
!     problem with both equality and inequality constraints, and, if the
!     user requests, obtains a covariance matrix of the solution
!     parameters.
!
!     Suppose there are given matrices E, A and G of respective
!     dimensions ME by N, MA by N and MG by N, and vectors F, B and H of
!     respective lengths ME, MA and MG.  This subroutine solves the
!     linearly constrained least squares problem
!
!                   EX = F, (E ME by N) (equations to be exactly
!                                       satisfied)
!                   AX = B, (A MA by N) (equations to be
!                                       approximately satisfied,
!                                       least squares sense)
!                   GX  >=  H,(G MG by N) (inequality constraints)
!
!     The inequalities GX  >=  H mean that every component of the
!     product GX must be  >=  the corresponding component of H.
!
!     In case the equality constraints cannot be satisfied, a
!     generalized inverse solution residual vector length is obtained
!     for F-EX.  This is the minimal length possible for F-EX.
!
!     Any values ME  >=  0, MA  >=  0, or MG  >=  0 are permitted.  The
!     rank of the matrix E is estimated during the computation.  We call
!     this value KRANKE.  It is an output parameter in IP(1) defined
!     below.  Using a generalized inverse solution of EX=F, a reduced
!     least squares problem with inequality constraints is obtained.
!     The tolerances used in these tests for determining the rank
!     of E and the rank of the reduced least squares problem are
!     given in Sandia Tech. Rept. SAND-78-1290.  They can be
!     modified by the user if new values are provided in
!     the option list of the array PRGOPT(*).
!
!     The user must dimension all arrays appearing in the call list..
!     W(MDW,N+1),PRGOPT(*),X(N),WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2)
!     where K=MAX(MA+MG,N).  This allows for a solution of a range of
!     problems in the given working space.  The dimension of WS(*)
!     given is a necessary overestimate.  Once a particular problem
!     has been run, the output parameter IP(3) gives the actual
!     dimension required for that problem.
!
!     The parameters for LSEI( ) are
!
!     Input..
!
!     W(*,*),MDW,   The array W(*,*) is doubly subscripted with
!     ME,MA,MG,N    first dimensioning parameter equal to MDW.
!                   For this discussion let us call M = ME+MA+MG.  Then
!                   MDW must satisfy MDW  >=  M.  The condition
!                   MDW  <  M is an error.
!
!                   The array W(*,*) contains the matrices and vectors
!
!                                  (E  F)
!                                  (A  B)
!                                  (G  H)
!
!                   in rows and columns 1,...,M and 1,...,N+1
!                   respectively.
!
!                   The integers ME, MA, and MG are the
!                   respective matrix row dimensions
!                   of E, A and G.  Each matrix has N columns.
!
!     PRGOPT(*)    This real-valued array is the option vector.
!                  If the user is satisfied with the nominal
!                  subprogram features set
!
!                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
!
!                  Otherwise PRGOPT(*) is a linked list consisting of
!                  groups of data of the following form
!
!                  LINK
!                  KEY
!                  DATA SET
!
!                  The parameters LINK and KEY are each one word.
!                  The DATA SET can be comprised of several words.
!                  The number of items depends on the value of KEY.
!                  The value of LINK points to the first
!                  entry of the next group of data within
!                  PRGOPT(*).  The exception is when there are
!                  no more options to change.  In that
!                  case, LINK=1 and the values KEY and DATA SET
!                  are not referenced.  The general layout of
!                  PRGOPT(*) is as follows.
!
!               ...PRGOPT(1) = LINK1 (link to first entry of next group)
!               .  PRGOPT(2) = KEY1 (key to the option change)
!               .  PRGOPT(3) = data value (data value for this change)
!               .       .
!               .       .
!               .       .
!               ...PRGOPT(LINK1)   = LINK2 (link to the first entry of
!               .                       next group)
!               .  PRGOPT(LINK1+1) = KEY2 (key to the option change)
!               .  PRGOPT(LINK1+2) = data value
!               ...     .
!               .       .
!               .       .
!               ...PRGOPT(LINK) = 1 (no more options to change)
!
!                  Values of LINK that are nonpositive are errors.
!                  A value of LINK  >  NLINK=100000 is also an error.
!                  This helps prevent using invalid but positive
!                  values of LINK that will probably extend
!                  beyond the program limits of PRGOPT(*).
!                  Unrecognized values of KEY are ignored.  The
!                  order of the options is arbitrary and any number
!                  of options can be changed with the following
!                  restriction.  To prevent cycling in the
!                  processing of the option array, a count of the
!                  number of options changed is maintained.
!                  Whenever this count exceeds NOPT=1000, an error
!                  message is printed and the subprogram returns.
!
!                  Options..
!
!                  KEY=1
!                         Compute in W(*,*) the N by N
!                  covariance matrix of the solution variables
!                  as an output parameter.  Nominally the
!                  covariance matrix will not be computed.
!                  (This requires no user input.)
!                  The data set for this option is a single value.
!                  It must be nonzero when the covariance matrix
!                  is desired.  If it is zero, the covariance
!                  matrix is not computed.  When the covariance matrix
!                  is computed, the first dimensioning parameter
!                  of the array W(*,*) must satisfy MDW  >=  MAX(M,N).
!
!                  KEY=10
!                         Suppress scaling of the inverse of the
!                  normal matrix by the scale factor RNORM**2/
!                  MAX(1, no. of degrees of freedom).  This option
!                  only applies when the option for computing the
!                  covariance matrix (KEY=1) is used.  With KEY=1 and
!                  KEY=10 used as options the unscaled inverse of the
!                  normal matrix is returned in W(*,*).
!                  The data set for this option is a single value.
!                  When it is nonzero no scaling is done.  When it is
!                  zero scaling is done.  The nominal case is to do
!                  scaling so if option (KEY=1) is used alone, the
!                  matrix will be scaled on output.
!
!                  KEY=2
!                         Scale the nonzero columns of the
!                         entire data matrix.
!                  (E)
!                  (A)
!                  (G)
!
!                  to have length one.  The data set for this
!                  option is a single value.  It must be
!                  nonzero if unit length column scaling
!                  is desired.
!
!                  KEY=3
!                         Scale columns of the entire data matrix
!                  (E)
!                  (A)
!                  (G)
!
!                  with a user-provided diagonal matrix.
!                  The data set for this option consists
!                  of the N diagonal scaling factors, one for
!                  each matrix column.
!
!                  KEY=4
!                         Change the rank determination tolerance for
!                  the equality constraint equations from
!                  the nominal value of SQRT(SRELPR).  This quantity can
!                  be no smaller than SRELPR, the arithmetic-
!                  storage precision.  The quantity SRELPR is the
!                  largest positive number such that T=1.+SRELPR
!                  satisfies T  ==  1.  The quantity used
!                  here is internally restricted to be at
!                  least SRELPR.  The data set for this option
!                  is the new tolerance.
!
!                  KEY=5
!                         Change the rank determination tolerance for
!                  the reduced least squares equations from
!                  the nominal value of SQRT(SRELPR).  This quantity can
!                  be no smaller than SRELPR, the arithmetic-
!                  storage precision.  The quantity used
!                  here is internally restricted to be at
!                  least SRELPR.  The data set for this option
!                  is the new tolerance.
!
!                  For example, suppose we want to change
!                  the tolerance for the reduced least squares
!                  problem, compute the covariance matrix of
!                  the solution parameters, and provide
!                  column scaling for the data matrix.  For
!                  these options the dimension of PRGOPT(*)
!                  must be at least N+9.  The Fortran statements
!                  defining these options would be as follows:
!
!                  PRGOPT(1)=4 (link to entry 4 in PRGOPT(*))
!                  PRGOPT(2)=1 (covariance matrix key)
!                  PRGOPT(3)=1 (covariance matrix wanted)
!
!                  PRGOPT(4)=7 (link to entry 7 in PRGOPT(*))
!                  PRGOPT(5)=5 (least squares equas.  tolerance key)
!                  PRGOPT(6)=... (new value of the tolerance)
!
!                  PRGOPT(7)=N+9 (link to entry N+9 in PRGOPT(*))
!                  PRGOPT(8)=3 (user-provided column scaling key)
!
!                  call SCOPY (N, D, 1, PRGOPT(9), 1)  (Copy the N
!                    scaling factors from the user array D(*)
!                    to PRGOPT(9)-PRGOPT(N+8))
!
!                  PRGOPT(N+9)=1 (no more options to change)
!
!                  The contents of PRGOPT(*) are not modified
!                  by the subprogram.
!                  The options for WNNLS( ) can also be included
!                  in this array.  The values of KEY recognized
!                  by WNNLS( ) are 6, 7 and 8.  Their functions
!                  are documented in the usage instructions for
!                  subroutine WNNLS( ).  Normally these options
!                  do not need to be modified when using LSEI( ).
!
!     IP(1),       The amounts of working storage actually
!     IP(2)        allocated for the working arrays WS(*) and
!                  IP(*), respectively.  These quantities are
!                  compared with the actual amounts of storage
!                  needed by LSEI( ).  Insufficient storage
!                  allocated for either WS(*) or IP(*) is an
!                  error.  This feature was included in LSEI( )
!                  because miscalculating the storage formulas
!                  for WS(*) and IP(*) might very well lead to
!                  subtle and hard-to-find execution errors.
!
!                  The length of WS(*) must be at least
!
!                  LW = 2*(ME+N)+K+(MG+2)*(N+7)
!
!                  where K = max(MA+MG,N)
!                  This test will not be made if IP(1) <= 0.
!
!                  The length of IP(*) must be at least
!
!                  LIP = MG+2*N+2
!                  This test will not be made if IP(2) <= 0.
!
!     Output..
!
!     X(*),RNORME,  The array X(*) contains the solution parameters
!     RNORML        if the integer output flag MODE = 0 or 1.
!                   The definition of MODE is given directly below.
!                   When MODE = 0 or 1, RNORME and RNORML
!                   respectively contain the residual vector
!                   Euclidean lengths of F - EX and B - AX.  When
!                   MODE=1 the equality constraint equations EX=F
!                   are contradictory, so RNORME  /=  0.  The residual
!                   vector F-EX has minimal Euclidean length.  For
!                   MODE  >=  2, none of these parameters is defined.
!
!     MODE          Integer flag that indicates the subprogram
!                   status after completion.  If MODE  >=  2, no
!                   solution has been computed.
!
!                   MODE =
!
!                   0  Both equality and inequality constraints
!                      are compatible and have been satisfied.
!
!                   1  Equality constraints are contradictory.
!                      A generalized inverse solution of EX=F was used
!                      to minimize the residual vector length F-EX.
!                      In this sense, the solution is still meaningful.
!
!                   2  Inequality constraints are contradictory.
!
!                   3  Both equality and inequality constraints
!                      are contradictory.
!
!                   The following interpretation of
!                   MODE=1,2 or 3 must be made.  The
!                   sets consisting of all solutions
!                   of the equality constraints EX=F
!                   and all vectors satisfying GX  >=  H
!                   have no points in common.  (In
!                   particular this does not say that
!                   each individual set has no points
!                   at all, although this could be the
!                   case.)
!
!                   4  Usage error occurred.  The value
!                      of MDW is  <  ME+MA+MG, MDW is
!                       <  N and a covariance matrix is
!                      requested, or the option vector
!                      PRGOPT(*) is not properly defined,
!                      or the lengths of the working arrays
!                      WS(*) and IP(*), when specified in
!                      IP(1) and IP(2) respectively, are not
!                      long enough.
!
!     W(*,*)        The array W(*,*) contains the N by N symmetric
!                   covariance matrix of the solution parameters,
!                   provided this was requested on input with
!                   the option vector PRGOPT(*) and the output
!                   flag is returned with MODE = 0 or 1.
!
!     IP(*)         The integer working array has three entries
!                   that provide rank and working array length
!                   information after completion.
!
!                      IP(1) = rank of equality constraint
!                              matrix.  Define this quantity
!                              as KRANKE.
!
!                      IP(2) = rank of reduced least squares
!                              problem.
!
!                      IP(3) = the amount of storage in the
!                              working array WS(*) that was
!                              actually used by the subprogram.
!                              The formula given above for the length
!                              of WS(*) is a necessary overestimate.
!                              If exactly the same problem matrices
!                              are used in subsequent executions,
!                              the declared dimension of WS(*) can
!                              be reduced to this output value.
!     User Designated
!     Working Arrays..
!
!     WS(*),IP(*)              These are respectively type real
!                              and type integer working arrays.
!                              Their required minimal lengths are
!                              given above.
!
!***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for
!                 linear least squares problems with equality and
!                 nonnegativity constraints, Report SAND77-0552, Sandia
!                 Laboratories, June 1978.
!               K. H. Haskell and R. J. Hanson, Selected algorithms for
!                 the linearly constrained least squares problem - a
!                 users guide, Report SAND78-1290, Sandia Laboratories,
!                 August 1979.
!               K. H. Haskell and R. J. Hanson, An algorithm for
!                 linear least squares problems with equality and
!                 nonnegativity constraints, Mathematical Programming
!                 21 (1981), pp. 98-118.
!               R. J. Hanson and K. H. Haskell, Two algorithms for the
!                 linearly constrained least squares problem, ACM
!                 Transactions on Mathematical Software, September 1982.
!***ROUTINES CALLED  H12, LSI, R1MACH, SASUM, SAXPY, SCOPY, SDOT, SNRM2,
!                    SSCAL, SSWAP, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890618  Completely restructured and extensively revised (WRB & RWC)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  LSEI
  INTEGER IP(3), MA, MDW, ME, MG, MODE, N
  REAL             PRGOPT(*), RNORME, RNORML, W(MDW,*), WS(*), X(*)
!
  EXTERNAL H12, LSI, R1MACH, SASUM, SAXPY, SCOPY, SDOT, SNRM2, &
     SSCAL, SSWAP, XERMSG
  REAL             R1MACH, SASUM, SDOT, SNRM2
!
  REAL             ENORM, FNORM, GAM, RB, RN, RNMAX, SIZE, SN, &
     SNMAX, SRELPR, T, TAU, UJ, UP, VJ, XNORM, XNRME
  INTEGER I, IMAX, J, JP1, K, KEY, KRANKE, LAST, LCHK, LINK, M, &
     MAPKE1, MDEQC, MEND, MEP1, N1, N2, NEXT, NLINK, NOPT, NP1, &
     NTIMES
  LOGICAL COV, FIRST
  CHARACTER*8 XERN1, XERN2, XERN3, XERN4
  SAVE FIRST, SRELPR
!
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  LSEI
!
!     Set the nominal tolerance used in the code for the equality
!     constraint equations.
!
  if (FIRST) SRELPR = R1MACH(4)
  FIRST = .FALSE.
  TAU = SQRT(SRELPR)
!
!     Check that enough storage was allocated in WS(*) and IP(*).
!
  MODE = 4
  if (MIN(N,ME,MA,MG)  <  0) THEN
     WRITE (XERN1, '(I8)') N
     WRITE (XERN2, '(I8)') ME
     WRITE (XERN3, '(I8)') MA
     WRITE (XERN4, '(I8)') MG
     call XERMSG ('SLATEC', 'LSEI', 'ALL OF THE VARIABLES N, ME,' // &
        ' MA, MG MUST BE  >=  0$$ENTERED ROUTINE WITH' // &
        '$$N  = ' // XERN1 // &
        '$$ME = ' // XERN2 // &
        '$$MA = ' // XERN3 // &
        '$$MG = ' // XERN4, 2, 1)
     return
  end if
!
  if (IP(1) > 0) THEN
     LCHK = 2*(ME+N) + MAX(MA+MG,N) + (MG+2)*(N+7)
     if (IP(1) < LCHK) THEN
        WRITE (XERN1, '(I8)') LCHK
        call XERMSG ('SLATEC', 'LSEI', 'INSUFFICIENT STORAGE ' // &
           'ALLOCATED FOR WS(*), NEED LW = ' // XERN1, 2, 1)
        return
     ENDIF
  end if
!
  if (IP(2) > 0) THEN
     LCHK = MG + 2*N + 2
     if (IP(2) < LCHK) THEN
        WRITE (XERN1, '(I8)') LCHK
        call XERMSG ('SLATEC', 'LSEI', 'INSUFFICIENT STORAGE ' // &
           'ALLOCATED FOR IP(*), NEED LIP = ' // XERN1, 2, 1)
        return
     ENDIF
  end if
!
!     Compute number of possible right multiplying Householder
!     transformations.
!
  M = ME + MA + MG
  if (N <= 0 .OR. M <= 0) THEN
     MODE = 0
     RNORME = 0
     RNORML = 0
     return
  end if
!
  if (MDW < M) THEN
     call XERMSG ('SLATEC', 'LSEI', 'MDW < ME+MA+MG IS AN ERROR', &
        2, 1)
     return
  end if
!
  NP1 = N + 1
  KRANKE = MIN(ME,N)
  N1 = 2*KRANKE + 1
  N2 = N1 + N
!
!     Set nominal values.
!
!     The nominal column scaling used in the code is
!     the identity scaling.
!
  call SCOPY (N, 1.E0, 0, WS(N1), 1)
!
!     No covariance matrix is nominally computed.
!
  COV = .FALSE.
!
!     Process option vector.
!     Define bound for number of options to change.
!
  NOPT = 1000
  NTIMES = 0
!
!     Define bound for positive values of LINK.
!
  NLINK = 100000
  LAST = 1
  LINK = PRGOPT(1)
  if (LINK == 0 .OR. LINK > NLINK) THEN
     call XERMSG ('SLATEC', 'LSEI', &
        'THE OPTION VECTOR IS UNDEFINED', 2, 1)
     return
  end if
!
  100 if (LINK > 1) THEN
     NTIMES = NTIMES + 1
     if (NTIMES > NOPT) THEN
        call XERMSG ('SLATEC', 'LSEI', &
           'THE LINKS IN THE OPTION VECTOR ARE CYCLING.', 2, 1)
        return
     ENDIF
!
     KEY = PRGOPT(LAST+1)
     if (KEY == 1) THEN
        COV = PRGOPT(LAST+2)  /=  0.E0
     ELSEIF (KEY == 2 .AND. PRGOPT(LAST+2) /= 0.E0) THEN
        DO 110 J = 1,N
           T = SNRM2(M,W(1,J),1)
           if (T /= 0.E0) T = 1.E0/T
           WS(J+N1-1) = T
  110       CONTINUE
     ELSEIF (KEY == 3) THEN
        call SCOPY (N, PRGOPT(LAST+2), 1, WS(N1), 1)
     ELSEIF (KEY == 4) THEN
        TAU = MAX(SRELPR,PRGOPT(LAST+2))
     ENDIF
!
     NEXT = PRGOPT(LINK)
     if (NEXT <= 0 .OR. NEXT > NLINK) THEN
     call XERMSG ('SLATEC', 'LSEI', &
        'THE OPTION VECTOR IS UNDEFINED', 2, 1)
        return
     ENDIF
!
     LAST = LINK
     LINK = NEXT
     go to 100
  end if
!
  DO 120 J = 1,N
     call SSCAL (M, WS(N1+J-1), W(1,J), 1)
  120 CONTINUE
!
  if (COV .AND. MDW < N) THEN
     call XERMSG ('SLATEC', 'LSEI', &
        'MDW  <  N WHEN COV MATRIX NEEDED, IS AN ERROR', 2, 1)
     return
  end if
!
!     Problem definition and option vector OK.
!
  MODE = 0
!
!     Compute norm of equality constraint matrix and right side.
!
  ENORM = 0.E0
  DO 130 J = 1,N
     ENORM = MAX(ENORM,SASUM(ME,W(1,J),1))
  130 CONTINUE
!
  FNORM = SASUM(ME,W(1,NP1),1)
  SNMAX = 0.E0
  RNMAX = 0.E0
  DO 150 I = 1,KRANKE
!
!        Compute maximum ratio of vector lengths. Partition is at
!        column I.
!
     DO 140 K = I,ME
        SN = SDOT(N-I+1,W(K,I),MDW,W(K,I),MDW)
        RN = SDOT(I-1,W(K,1),MDW,W(K,1),MDW)
        if (RN == 0.E0 .AND. SN > SNMAX) THEN
           SNMAX = SN
           IMAX = K
        ELSEIF (K == I .OR. SN*RNMAX > RN*SNMAX) THEN
           SNMAX = SN
           RNMAX = RN
           IMAX = K
        ENDIF
  140    CONTINUE
!
!        Interchange rows if necessary.
!
     if (I /= IMAX) call SSWAP (NP1, W(I,1), MDW, W(IMAX,1), MDW)
     if (SNMAX > RNMAX*TAU**2) THEN
!
!        Eliminate elements I+1,...,N in row I.
!
        call H12 (1, I, I+1, N, W(I,1), MDW, WS(I), W(I+1,1), MDW, &
                  1, M-I)
     ELSE
        KRANKE = I - 1
        go to 160
     ENDIF
  150 CONTINUE
!
!     Save diagonal terms of lower trapezoidal matrix.
!
  160 call SCOPY (KRANKE, W, MDW+1, WS(KRANKE+1), 1)
!
!     Use Householder transformation from left to achieve
!     KRANKE by KRANKE upper triangular form.
!
  if (KRANKE < ME) THEN
     DO 170 K = KRANKE,1,-1
!
!           Apply transformation to matrix cols. 1,...,K-1.
!
        call H12 (1, K, KRANKE+1, ME, W(1,K), 1, UP, W, 1, MDW, K-1)
!
!           Apply to rt side vector.
!
        call H12 (2, K, KRANKE+1, ME, W(1,K), 1, UP, W(1,NP1), 1, 1, &
                  1)
  170    CONTINUE
  end if
!
!     Solve for variables 1,...,KRANKE in new coordinates.
!
  call SCOPY (KRANKE, W(1, NP1), 1, X, 1)
  DO 180 I = 1,KRANKE
     X(I) = (X(I)-SDOT(I-1,W(I,1),MDW,X,1))/W(I,I)
  180 CONTINUE
!
!     Compute residuals for reduced problem.
!
  MEP1 = ME + 1
  RNORML = 0.E0
  DO 190 I = MEP1,M
     W(I,NP1) = W(I,NP1) - SDOT(KRANKE,W(I,1),MDW,X,1)
     SN = SDOT(KRANKE,W(I,1),MDW,W(I,1),MDW)
     RN = SDOT(N-KRANKE,W(I,KRANKE+1),MDW,W(I,KRANKE+1),MDW)
     if (RN <= SN*TAU**2 .AND. KRANKE < N) &
        call SCOPY (N-KRANKE, 0.E0, 0, W(I,KRANKE+1), MDW)
  190 CONTINUE
!
!     Compute equality constraint equations residual length.
!
  RNORME = SNRM2(ME-KRANKE,W(KRANKE+1,NP1),1)
!
!     Move reduced problem data upward if KRANKE < ME.
!
  if (KRANKE < ME) THEN
     DO 200 J = 1,NP1
        call SCOPY (M-ME, W(ME+1,J), 1, W(KRANKE+1,J), 1)
  200    CONTINUE
  end if
!
!     Compute solution of reduced problem.
!
  call LSI(W(KRANKE+1, KRANKE+1), MDW, MA, MG, N-KRANKE, PRGOPT, &
           X(KRANKE+1), RNORML, MODE, WS(N2), IP(2))
!
!     Test for consistency of equality constraints.
!
  if (ME > 0) THEN
     MDEQC = 0
     XNRME = SASUM(KRANKE,W(1,NP1),1)
     if (RNORME > TAU*(ENORM*XNRME+FNORM)) MDEQC = 1
     MODE = MODE + MDEQC
!
!        Check if solution to equality constraints satisfies inequality
!        constraints when there are no degrees of freedom left.
!
     if (KRANKE == N .AND. MG > 0) THEN
        XNORM = SASUM(N,X,1)
        MAPKE1 = MA + KRANKE + 1
        MEND = MA + KRANKE + MG
        DO 210 I = MAPKE1,MEND
           SIZE = SASUM(N,W(I,1),MDW)*XNORM + ABS(W(I,NP1))
           if (W(I,NP1) > TAU*SIZE) THEN
              MODE = MODE + 2
              go to 290
           ENDIF
  210       CONTINUE
     ENDIF
  end if
!
!     Replace diagonal terms of lower trapezoidal matrix.
!
  if (KRANKE > 0) THEN
     call SCOPY (KRANKE, WS(KRANKE+1), 1, W, MDW+1)
!
!        Reapply transformation to put solution in original coordinates.
!
     DO 220 I = KRANKE,1,-1
        call H12 (2, I, I+1, N, W(I,1), MDW, WS(I), X, 1, 1, 1)
  220    CONTINUE
!
!        Compute covariance matrix of equality constrained problem.
!
     if (COV) THEN
        DO 270 J = MIN(KRANKE,N-1),1,-1
           RB = WS(J)*W(J,J)
           if (RB /= 0.E0) RB = 1.E0/RB
           JP1 = J + 1
           DO 230 I = JP1,N
              W(I,J) = RB*SDOT(N-J,W(I,JP1),MDW,W(J,JP1),MDW)
  230          CONTINUE
!
           GAM = 0.5E0*RB*SDOT(N-J,W(JP1,J),1,W(J,JP1),MDW)
           call SAXPY (N-J, GAM, W(J,JP1), MDW, W(JP1,J), 1)
           DO 250 I = JP1,N
              DO 240 K = I,N
                 W(I,K) = W(I,K) + W(J,I)*W(K,J) + W(I,J)*W(J,K)
                 W(K,I) = W(I,K)
  240             CONTINUE
  250          CONTINUE
           UJ = WS(J)
           VJ = GAM*UJ
           W(J,J) = UJ*VJ + UJ*VJ
           DO 260 I = JP1,N
              W(J,I) = UJ*W(I,J) + VJ*W(J,I)
  260          CONTINUE
           call SCOPY (N-J, W(J, JP1), MDW, W(JP1,J), 1)
  270       CONTINUE
     ENDIF
  end if
!
!     Apply the scaling to the covariance matrix.
!
  if (COV) THEN
     DO 280 I = 1,N
        call SSCAL (N, WS(I+N1-1), W(I,1), MDW)
        call SSCAL (N, WS(I+N1-1), W(1,I), 1)
  280    CONTINUE
  end if
!
!     Rescale solution vector.
!
  290 if (MODE <= 1) THEN
     DO 300 J = 1,N
        X(J) = X(J)*WS(N1+J-1)
  300    CONTINUE
  end if
!
  IP(1) = KRANKE
  IP(3) = IP(3) + 2*KRANKE + N
  return
end
