subroutine SGMRES (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, &
     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, RGWK, LRGW, &
     IGWK, LIGW, RWORK, IWORK)
!
!! SGMRES is a Preconditioned GMRES Iterative Sparse Ax=b Solver.
!
!            This routine uses the generalized minimum residual
!            (GMRES) method with preconditioning to solve
!            non-symmetric linear systems of the form: Ax = b.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A4, D2B4
!***TYPE      SINGLE PRECISION (SGMRES-S, DGMRES-D)
!***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
!             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
!***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
!           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
!           Seager, Mark K., (LLNL), seager@llnl.gov
!             Lawrence Livermore National Laboratory
!             PO Box 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!***DESCRIPTION
!
! *Usage:
!      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
!      INTEGER   ITER, IERR, IUNIT, LRGW, IGWK(LIGW), LIGW
!      INTEGER   IWORK(USER DEFINED)
!      REAL      B(N), X(N), A(NELT), TOL, ERR, SB(N), SX(N)
!      REAL      RGWK(LRGW), RWORK(USER DEFINED)
!      EXTERNAL  MATVEC, MSOLVE
!
!      call SGMRES(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
!     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX,
!     $     RGWK, LRGW, IGWK, LIGW, RWORK, IWORK)
!
! *Arguments:
! N      :IN       Integer.
!         Order of the Matrix.
! B      :IN       Real B(N).
!         Right-hand side vector.
! X      :INOUT    Real X(N).
!         On input X is your initial guess for the solution vector.
!         On output X is the final approximate solution.
! NELT   :IN       Integer.
!         Number of Non-Zeros stored in A.
! IA     :IN       Integer IA(NELT).
! JA     :IN       Integer JA(NELT).
! A      :IN       Real A(NELT).
!         These arrays contain the matrix data structure for A.
!         It could take any form.  See "Description", below,
!         for more details.
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all non-zero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the upper
!         or lower triangle of the matrix is stored.
! MATVEC :EXT      External.
!         Name of a routine which performs the matrix vector multiply
!         Y = A*X given A and X.  The name of the MATVEC routine must
!         be declared external in the calling program.  The calling
!         sequence to MATVEC is:
!             call MATVEC(N, X, Y, NELT, IA, JA, A, ISYM)
!         where N is the number of unknowns, Y is the product A*X
!         upon return, X is an input vector, and NELT is the number of
!         non-zeros in the SLAP IA, JA, A storage for the matrix A.
!         ISYM is a flag which, if non-zero, denotes that A is
!         symmetric and only the lower or upper triangle is stored.
! MSOLVE :EXT      External.
!         Name of the routine which solves a linear system Mz = r for
!         z given r with the preconditioning matrix M (M is supplied via
!         RWORK and IWORK arrays.  The name of the MSOLVE routine must
!         be declared external in the calling program.  The calling
!         sequence to MSOLVE is:
!             call MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!         Where N is the number of unknowns, R is the right-hand side
!         vector and Z is the solution upon return.  NELT, IA, JA, A and
!         ISYM are defined as above.  RWORK is a real array that can
!         be used to pass necessary preconditioning information and/or
!         workspace to MSOLVE.  IWORK is an integer work array for
!         the same purpose as RWORK.
! ITOL   :IN       Integer.
!         Flag to indicate the type of convergence criterion used.
!         ITOL=0  Means the  iteration stops when the test described
!                 below on  the  residual RL  is satisfied.  This is
!                 the  "Natural Stopping Criteria" for this routine.
!                 Other values  of   ITOL  cause  extra,   otherwise
!                 unnecessary, computation per iteration and     are
!                 therefore  much less  efficient.  See  ISSGMR (the
!                 stop test routine) for more information.
!         ITOL=1  Means   the  iteration stops   when the first test
!                 described below on  the residual RL  is satisfied,
!                 and there  is either right  or  no preconditioning
!                 being used.
!         ITOL=2  Implies     that   the  user    is   using    left
!                 preconditioning, and the second stopping criterion
!                 below is used.
!         ITOL=3  Means the  iteration stops   when  the  third test
!                 described below on Minv*Residual is satisfied, and
!                 there is either left  or no  preconditioning being
!                 used.
!         ITOL=11 is    often  useful  for   checking  and comparing
!                 different routines.  For this case, the  user must
!                 supply  the  "exact" solution or  a  very accurate
!                 approximation (one with  an  error much less  than
!                 TOL) through a common block,
!                     COMMON /SSLBLK/ SOLN( )
!                 If ITOL=11, iteration stops when the 2-norm of the
!                 difference between the iterative approximation and
!                 the user-supplied solution  divided by the  2-norm
!                 of the  user-supplied solution  is  less than TOL.
!                 Note that this requires  the  user to  set up  the
!                 "COMMON     /SSLBLK/ SOLN(LENGTH)"  in the calling
!                 routine.  The routine with this declaration should
!                 be loaded before the stop test so that the correct
!                 length is used by  the loader.  This procedure  is
!                 not standard Fortran and may not work correctly on
!                 your   system (although  it  has  worked  on every
!                 system the authors have tried).  If ITOL is not 11
!                 then this common block is indeed standard Fortran.
! TOL    :INOUT    Real.
!         Convergence criterion, as described below.  If TOL is set
!         to zero on input, then a default value of 500*(the smallest
!         positive magnitude, machine epsilon) is used.
! ITMAX  :DUMMY    Integer.
!         Maximum number of iterations in most SLAP routines.  In
!         this routine this does not make sense.  The maximum number
!         of iterations here is given by ITMAX = MAXL*(NRMAX+1).
!         See IGWK for definitions of MAXL and NRMAX.
! ITER   :OUT      Integer.
!         Number of iterations required to reach convergence, or
!         ITMAX if convergence criterion could not be achieved in
!         ITMAX iterations.
! ERR    :OUT      Real.
!         Error estimate of error in final approximate solution, as
!         defined by ITOL.  Letting norm() denote the Euclidean
!         norm, ERR is defined as follows..
!
!         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
!                               for right or no preconditioning, and
!                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
!                                norm(SB*(M-inverse)*B),
!                               for left preconditioning.
!         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
!                               since right or no preconditioning
!                               being used.
!         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
!                                norm(SB*(M-inverse)*B),
!                               since left preconditioning is being
!                               used.
!         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
!                               i=1,n
!         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
! IERR   :OUT      Integer.
!         Return error flag.
!               IERR = 0 => All went well.
!               IERR = 1 => Insufficient storage allocated for
!                           RGWK or IGWK.
!               IERR = 2 => Routine SGMRES failed to reduce the norm
!                           of the current residual on its last call,
!                           and so the iteration has stalled.  In
!                           this case, X equals the last computed
!                           approximation.  The user must either
!                           increase MAXL, or choose a different
!                           initial guess.
!               IERR =-1 => Insufficient length for RGWK array.
!                           IGWK(6) contains the required minimum
!                           length of the RGWK array.
!               IERR =-2 => Illegal value of ITOL, or ITOL and JPRE
!                           values are inconsistent.
!         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the
!         left-hand-side of the relevant stopping test defined
!         below associated with the residual for the current
!         approximation X(L).
! IUNIT  :IN       Integer.
!         Unit number on which to write the error at each iteration,
!         if this is desired for monitoring convergence.  If unit
!         number is 0, no writing will occur.
! SB     :IN       Real SB(N).
!         Array of length N containing scale factors for the right
!         hand side vector B.  If JSCAL.eq.0 (see below), SB need
!         not be supplied.
! SX     :IN       Real SX(N).
!         Array of length N containing scale factors for the solution
!         vector X.  If JSCAL.eq.0 (see below), SX need not be
!         supplied.  SB and SX can be the same array in the calling
!         program if desired.
! RGWK   :INOUT    Real RGWK(LRGW).
!         Real array used for workspace by SGMRES.
!         On return, RGWK(1) = RHOL.  See IERR for definition of RHOL.
! LRGW   :IN       Integer.
!         Length of the real workspace, RGWK.
!         LRGW >= 1 + N*(MAXL+6) + MAXL*(MAXL+3).
!         See below for definition of MAXL.
!         For the default values, RGWK has size at least 131 + 16*N.
! IGWK   :INOUT    Integer IGWK(LIGW).
!         The following IGWK parameters should be set by the user
!         before calling this routine.
!         IGWK(1) = MAXL.  Maximum dimension of Krylov subspace in
!            which X - X0 is to be found (where, X0 is the initial
!            guess).  The default value of MAXL is 10.
!         IGWK(2) = KMP.  Maximum number of previous Krylov basis
!            vectors to which each new basis vector is made orthogonal.
!            The default value of KMP is MAXL.
!         IGWK(3) = JSCAL.  Flag indicating whether the scaling
!            arrays SB and SX are to be used.
!            JSCAL = 0 => SB and SX are not used and the algorithm
!               will perform as if all SB(I) = 1 and SX(I) = 1.
!            JSCAL = 1 =>  Only SX is used, and the algorithm
!               performs as if all SB(I) = 1.
!            JSCAL = 2 =>  Only SB is used, and the algorithm
!               performs as if all SX(I) = 1.
!            JSCAL = 3 =>  Both SB and SX are used.
!         IGWK(4) = JPRE.  Flag indicating whether preconditioning
!            is being used.
!            JPRE = 0  =>  There is no preconditioning.
!            JPRE > 0  =>  There is preconditioning on the right
!               only, and the solver will call routine MSOLVE.
!            JPRE < 0  =>  There is preconditioning on the left
!               only, and the solver will call routine MSOLVE.
!         IGWK(5) = NRMAX.  Maximum number of restarts of the
!            Krylov iteration.  The default value of NRMAX = 10.
!            if IWORK(5) = -1,  then no restarts are performed (in
!            this case, NRMAX is set to zero internally).
!         The following IWORK parameters are diagnostic information
!         made available to the user after this routine completes.
!         IGWK(6) = MLWK.  Required minimum length of RGWK array.
!         IGWK(7) = NMS.  The total number of calls to MSOLVE.
! LIGW   :IN       Integer.
!         Length of the integer workspace, IGWK.  LIGW >= 20.
! RWORK  :WORK     Real RWORK(USER DEFINED).
!         Real array that can be used for workspace in MSOLVE.
! IWORK  :WORK     Integer IWORK(USER DEFINED).
!         Integer array that can be used for workspace in MSOLVE.
!
! *Description:
!       SGMRES solves a linear system A*X = B rewritten in the form:
!
!        (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B,
!
!       with right preconditioning, or
!
!        (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B,
!
!       with left preconditioning, where A is an N-by-N real matrix,
!       X  and  B are N-vectors,   SB and SX   are  diagonal scaling
!       matrices,   and M is  a preconditioning    matrix.   It uses
!       preconditioned  Krylov   subpace  methods  based     on  the
!       generalized minimum residual  method (GMRES).   This routine
!       optionally performs  either  the  full     orthogonalization
!       version of the  GMRES  algorithm or an incomplete variant of
!       it.  Both versions use restarting of the linear iteration by
!       default, although the user can disable this feature.
!
!       The GMRES  algorithm generates a sequence  of approximations
!       X(L) to the  true solution of the above  linear system.  The
!       convergence criteria for stopping the  iteration is based on
!       the size  of the  scaled norm of  the residual  R(L)  =  B -
!       A*X(L).  The actual stopping test is either:
!
!               norm(SB*(B-A*X(L))) .le. TOL*norm(SB*B),
!
!       for right preconditioning, or
!
!               norm(SB*(M-inverse)*(B-A*X(L))) .le.
!                       TOL*norm(SB*(M-inverse)*B),
!
!       for left preconditioning, where norm() denotes the Euclidean
!       norm, and TOL is  a positive scalar less  than one  input by
!       the user.  If TOL equals zero  when SGMRES is called, then a
!       default  value  of 500*(the   smallest  positive  magnitude,
!       machine epsilon) is used.  If the  scaling arrays SB  and SX
!       are used, then  ideally they  should be chosen  so  that the
!       vectors SX*X(or SX*M*X) and  SB*B have all their  components
!       approximately equal  to  one in  magnitude.  If one wants to
!       use the same scaling in X  and B, then  SB and SX can be the
!       same array in the calling program.
!
!       The following is a list of the other routines and their
!       functions used by SGMRES:
!       SPIGMR  Contains the main iteration loop for GMRES.
!       SORTH   Orthogonalizes a new vector against older basis vectors.
!       SHEQR   Computes a QR decomposition of a Hessenberg matrix.
!       SHELS   Solves a Hessenberg least-squares system, using QR
!               factors.
!       SRLCAL  Computes the scaled residual RL.
!       SXLCAL  Computes the solution XL.
!       ISSGMR  User-replaceable stopping routine.
!
!       This routine does  not care  what matrix data   structure is
!       used for  A and M.  It simply   calls  the MATVEC and MSOLVE
!       routines, with  the arguments as  described above.  The user
!       could write any type of structure and the appropriate MATVEC
!       and MSOLVE routines.  It is assumed  that A is stored in the
!       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is
!       stored  in  IWORK  and  RWORK   in  some fashion.   The SLAP
!       routines SSDCG and SSICCG are examples of this procedure.
!
!       Two  examples  of  matrix  data structures  are the: 1) SLAP
!       Triad  format and 2) SLAP Column format.
!
!       =================== S L A P Triad format ===================
!       This routine requires that the  matrix A be   stored in  the
!       SLAP  Triad format.  In  this format only the non-zeros  are
!       stored.  They may appear in  *ANY* order.  The user supplies
!       three arrays of  length NELT, where  NELT is  the number  of
!       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
!       each non-zero the user puts the row and column index of that
!       matrix element  in the IA and  JA arrays.  The  value of the
!       non-zero   matrix  element is  placed  in  the corresponding
!       location of the A array.   This is  an  extremely  easy data
!       structure to generate.  On  the  other hand it   is  not too
!       efficient on vector computers for  the iterative solution of
!       linear systems.  Hence,   SLAP changes   this  input    data
!       structure to the SLAP Column format  for  the iteration (but
!       does not change it back).
!
!       Here is an example of the  SLAP Triad   storage format for a
!       5x5 Matrix.  Recall that the entries may appear in any order.
!
!           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.
!                              1  2  3  4  5  6  7  8  91011
!       |1112  0  015|   A: 5112113315535522354421
!       |2122  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
!       | 0  033  035|  JA:  1  2  1  3  5  3  5  2  5  4  1
!       | 0  0  044  0|
!       |51  053  055|
!
!       =================== S L A P Column format ==================
!
!       This routine  requires that  the matrix A  be stored in  the
!       SLAP Column format.  In this format the non-zeros are stored
!       counting down columns (except for  the diagonal entry, which
!       must appear first in each  "column")  and are stored  in the
!       real array A.  In other words, for each column in the matrix
!       put the diagonal entry in A.  Then put in the other non-zero
!       elements going down   the  column (except  the diagonal)  in
!       order.  The IA array holds the row  index for each non-zero.
!       The JA array holds the offsets into the IA, A arrays for the
!       beginning of   each    column.    That  is,    IA(JA(ICOL)),
!       A(JA(ICOL)) points to the beginning of the ICOL-th column in
!       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
!       end  of   the ICOL-th  column.  Note   that  we  always have
!       JA(N+1) = NELT+1, where  N  is the number of columns in  the
!       matrix and  NELT   is the number of non-zeros in the matrix.
!
!       Here is an example of the  SLAP Column  storage format for a
!       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a
!       column):
!
!           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
!                              1  2  3    4  5    6  7    8    91011
!       |1112  0  015|   A: 112151 | 2212 | 3353 | 44 | 551535
!       |2122  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
!       | 0  033  035|  JA:  1  4  6    8  9   12
!       | 0  0  044  0|
!       |51  053  055|
!
! *Cautions:
!     This routine will attempt to write to the Fortran logical output
!     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
!     this logical unit is attached to a file or terminal before calling
!     this routine with a non-zero value for IUNIT.  This routine does
!     not check for the validity of a non-zero IUNIT unit number.
!
!***REFERENCES  1. Peter N. Brown and A. C. Hindmarsh, Reduced Storage
!                  Matrix Methods in Stiff ODE Systems, Lawrence Liver-
!                  more National Laboratory Report UCRL-95088, Rev. 1,
!                  Livermore, California, June 1987.
!               2. Mark K. Seager, A SLAP for the Masses, in
!                  G. F. Carey, Ed., Parallel Supercomputing: Methods,
!                  Algorithms and Applications, Wiley, 1989, pp.135-155.
!***ROUTINES CALLED  R1MACH, SCOPY, SNRM2, SPIGMR
!***REVISION HISTORY  (YYMMDD)
!   871001  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   891004  Added new reference.
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910506  Corrected errors in C***ROUTINES CALLED list.  (FNF)
!   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
!   920511  Added complete declaration section.  (WRB)
!   920929  Corrected format of references.  (FNF)
!   921019  Changed 500.0 to 500 to reduce SP/DP differences.  (FNF)
!   921026  Added check for valid value of ITOL.  (FNF)
!***END PROLOGUE  SGMRES
!         The following is for optimized compilation on LLNL/LTSS Crays.
!LLL. OPTIMIZE
!     .. Scalar Arguments ..
  REAL ERR, TOL
  INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LIGW, LRGW, N, NELT
!     .. Array Arguments ..
  REAL A(NELT), B(N), RGWK(LRGW), RWORK(*), SB(N), SX(N), X(N)
  INTEGER IA(NELT), IGWK(LIGW), IWORK(*), JA(NELT)
!     .. Subroutine Arguments ..
  EXTERNAL MATVEC, MSOLVE
!     .. Local Scalars ..
  REAL BNRM, RHOL, SUM
  INTEGER I, IFLAG, JPRE, JSCAL, KMP, LDL, LGMR, LHES, LQ, LR, LV, &
          LW, LXL, LZ, LZM1, MAXL, MAXLP1, NMS, NMSL, NRMAX, NRSTS
!     .. External Functions ..
  REAL R1MACH, SNRM2
  EXTERNAL R1MACH, SNRM2
!     .. External Subroutines ..
  EXTERNAL SCOPY, SPIGMR
!     .. Intrinsic Functions ..
  INTRINSIC SQRT
!***FIRST EXECUTABLE STATEMENT  SGMRES
  IERR = 0
!   ------------------------------------------------------------------
!         Load method parameters with user values or defaults.
!   ------------------------------------------------------------------
  MAXL = IGWK(1)
  if (MAXL  ==  0) MAXL = 10
  if (MAXL  >  N) MAXL = N
  KMP = IGWK(2)
  if (KMP  ==  0) KMP = MAXL
  if (KMP  >  MAXL) KMP = MAXL
  JSCAL = IGWK(3)
  JPRE = IGWK(4)
!         Check for valid value of ITOL.
  if (  (ITOL < 0) .OR. ((ITOL > 3).AND.(ITOL /= 11)) ) GOTO 650
!         Check for consistent values of ITOL and JPRE.
  if (  ITOL == 1 .AND. JPRE < 0 ) GOTO 650
  if (  ITOL == 2 .AND. JPRE >= 0 ) GOTO 650
  NRMAX = IGWK(5)
  if (  NRMAX == 0 ) NRMAX = 10
!         If NRMAX .eq. -1, then set NRMAX = 0 to turn off restarting.
  if (  NRMAX == -1 ) NRMAX = 0
!         If input value of TOL is zero, set it to its default value.
  if (  TOL == 0.0E0 ) TOL = 500*R1MACH(3)
!
!         Initialize counters.
  ITER = 0
  NMS = 0
  NRSTS = 0
!   ------------------------------------------------------------------
!         Form work array segment pointers.
!   ------------------------------------------------------------------
  MAXLP1 = MAXL + 1
  LV = 1
  LR = LV + N*MAXLP1
  LHES = LR + N + 1
  LQ = LHES + MAXL*MAXLP1
  LDL = LQ + 2*MAXL
  LW = LDL + N
  LXL = LW + N
  LZ = LXL + N
!
!         Load IGWK(6) with required minimum length of the RGWK array.
  IGWK(6) = LZ + N - 1
  if (  LZ+N-1 > LRGW ) GOTO 640
!   ------------------------------------------------------------------
!         Calculate scaled-preconditioned norm of RHS vector b.
!   ------------------------------------------------------------------
  if (JPRE  <  0) THEN
     call MSOLVE(N, B, RGWK(LR), NELT, IA, JA, A, ISYM, &
          RWORK, IWORK)
     NMS = NMS + 1
  ELSE
     call SCOPY(N, B, 1, RGWK(LR), 1)
  end if
  if (  JSCAL == 2 .OR. JSCAL == 3 ) THEN
     SUM = 0
     DO 10 I = 1,N
        SUM = SUM + (RGWK(LR-1+I)*SB(I))**2
 10      CONTINUE
     BNRM = SQRT(SUM)
  ELSE
     BNRM = SNRM2(N,RGWK(LR),1)
  end if
!   ------------------------------------------------------------------
!         Calculate initial residual.
!   ------------------------------------------------------------------
  call MATVEC(N, X, RGWK(LR), NELT, IA, JA, A, ISYM)
  DO 50 I = 1,N
     RGWK(LR-1+I) = B(I) - RGWK(LR-1+I)
 50   CONTINUE
!   ------------------------------------------------------------------
!         If performing restarting, then load the residual into the
!         correct location in the RGWK array.
!   ------------------------------------------------------------------
 100  CONTINUE
  if (  NRSTS > NRMAX ) GOTO 610
  if (  NRSTS > 0 ) THEN
!         Copy the current residual to a different location in the RGWK
!         array.
     call SCOPY(N, RGWK(LDL), 1, RGWK(LR), 1)
  end if
!   ------------------------------------------------------------------
!         Use the SPIGMR algorithm to solve the linear system A*Z = R.
!   ------------------------------------------------------------------
  call SPIGMR(N, RGWK(LR), SB, SX, JSCAL, MAXL, MAXLP1, KMP, &
         NRSTS, JPRE, MATVEC, MSOLVE, NMSL, RGWK(LZ), RGWK(LV), &
         RGWK(LHES), RGWK(LQ), LGMR, RWORK, IWORK, RGWK(LW), &
         RGWK(LDL), RHOL, NRMAX, B, BNRM, X, RGWK(LXL), ITOL, &
         TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)
  ITER = ITER + LGMR
  NMS = NMS + NMSL
!
!         Increment X by the current approximate solution Z of A*Z = R.
!
  LZM1 = LZ - 1
  DO 110 I = 1,N
     X(I) = X(I) + RGWK(LZM1+I)
 110  CONTINUE
  if (  IFLAG == 0 ) GOTO 600
  if (  IFLAG == 1 ) THEN
     NRSTS = NRSTS + 1
     GOTO 100
  end if
  if (  IFLAG == 2 ) GOTO 620
!   ------------------------------------------------------------------
!         All returns are made through this section.
!   ------------------------------------------------------------------
!         The iteration has converged.
!
 600  CONTINUE
  IGWK(7) = NMS
  RGWK(1) = RHOL
  IERR = 0
  return
!
!         Max number((NRMAX+1)*MAXL) of linear iterations performed.
 610  CONTINUE
  IGWK(7) = NMS
  RGWK(1) = RHOL
  IERR = 1
  return
!
!         GMRES failed to reduce last residual in MAXL iterations.
!         The iteration has stalled.
 620  CONTINUE
  IGWK(7) = NMS
  RGWK(1) = RHOL
  IERR = 2
  return
!         Error return.  Insufficient length for RGWK array.
 640  CONTINUE
  ERR = TOL
  IERR = -1
  return
!         Error return.  Inconsistent ITOL and JPRE values.
 650  CONTINUE
  ERR = TOL
  IERR = -2
  return
!------------- LAST LINE OF SGMRES FOLLOWS ----------------------------
end
