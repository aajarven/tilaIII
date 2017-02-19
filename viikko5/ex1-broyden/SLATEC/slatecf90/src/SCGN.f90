subroutine SCGN (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC, &
     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, ATP, &
     ATZ, DZ, ATDZ, RWORK, IWORK)
!
!! SCGN is the Preconditioned CG Sparse Ax=b Solver for Normal Equations.
!
!            Routine to solve a general linear system  Ax = b  using the
!            Preconditioned Conjugate Gradient method applied to the
!            normal equations  AA'y = b, x=A'y.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A4, D2B4
!***TYPE      SINGLE PRECISION (SCGN-S, DCGN-D)
!***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM SOLVE,
!             NORMAL EQUATIONS., SLAP, SPARSE
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!
! *Usage:
!     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
!     INTEGER  ITER, IERR, IUNIT, IWORK(USER DEFINED)
!     REAL     B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N)
!     REAL     P(N), ATP(N), ATZ(N), DZ(N), ATDZ(N)
!     REAL     RWORK(USER DEFINED)
!     EXTERNAL MATVEC, MTTVEC, MSOLVE
!
!     call SCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC,
!    $     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R,
!    $     Z, P, ATP, ATZ, DZ, ATDZ, RWORK, IWORK)
!
! *Arguments:
! N      :IN       Integer
!         Order of the Matrix.
! B      :IN       Real B(N).
!         Right-hand side vector.
! X      :INOUT    Real X(N).
!         On input X is your initial guess for solution vector.
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
!         y = A*X given A and X.  The name of the MATVEC routine must
!         be declared external in the calling program.  The calling
!         sequence to MATVEC is:
!             call MATVEC( N, X, Y, NELT, IA, JA, A, ISYM )
!         Where N is the number of unknowns, Y is the product A*X
!         upon return X is an input vector, NELT is the number of
!         non-zeros in the SLAP-Column IA, JA, A storage for the matrix
!         A.  ISYM is a flag which, if non-zero, denotes that A is
!         symmetric and only the lower or upper triangle is stored.
! MTTVEC :EXT      External.
!         Name of a routine which performs the matrix transpose vector
!         multiply y = A'*X given A and X (where ' denotes transpose).
!         The name of the MTTVEC routine must be declared external in
!         the calling program.  The calling sequence to MTTVEC is the
!         same as that for MATVEC, viz.:
!             call MTTVEC( N, X, Y, NELT, IA, JA, A, ISYM )
!         Where N is the number of unknowns, Y is the product A'*X
!         upon return X is an input vector, NELT is the number of
!         non-zeros in the SLAP-Column IA, JA, A storage for the matrix
!         A.  ISYM is a flag which, if non-zero, denotes that A is
!         symmetric and only the lower or upper triangle is stored.
! MSOLVE :EXT      External.
!         Name of a routine which solves a linear system MZ = R for
!         Z given R with the preconditioning matrix M (M is supplied via
!         RWORK and IWORK arrays).  The name of the MSOLVE routine must
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
!         Flag to indicate type of convergence criterion.
!         If ITOL=1, iteration stops when the 2-norm of the residual
!         divided by the 2-norm of the right-hand side is less than TOL.
!         If ITOL=2, iteration stops when the 2-norm of M-inv times the
!         residual divided by the 2-norm of M-inv times the right hand
!         side is less than TOL, where M-inv is the inverse of the
!         diagonal of A.
!         ITOL=11 is often useful for checking and comparing different
!         routines.  For this case, the user must supply the "exact"
!         solution or a very accurate approximation (one with an error
!         much less than TOL) through a common block,
!             COMMON /SSLBLK/ SOLN( )
!         If ITOL=11, iteration stops when the 2-norm of the difference
!         between the iterative approximation and the user-supplied
!         solution divided by the 2-norm of the user-supplied solution
!         is less than TOL.  Note that this requires the user to set up
!         the "COMMON /SSLBLK/ SOLN(LENGTH)" in the calling routine.
!         The routine with this declaration should be loaded before the
!         stop test so that the correct length is used by the loader.
!         This procedure is not standard Fortran and may not work
!         correctly on your system (although it has worked on every
!         system the authors have tried).  If ITOL is not 11 then this
!         common block is indeed standard Fortran.
! TOL    :INOUT    Real.
!         Convergence criterion, as described above.  (Reset if IERR=4.)
! ITMAX  :IN       Integer.
!         Maximum number of iterations.
! ITER   :OUT      Integer.
!         Number of iterations required to reach convergence, or
!         ITMAX+1 if convergence criterion could not be achieved in
!         ITMAX iterations.
! ERR    :OUT      Real.
!         Error estimate of error in final approximate solution, as
!         defined by ITOL.
! IERR   :OUT      Integer.
!         Return error flag.
!           IERR = 0 => All went well.
!           IERR = 1 => Insufficient space allocated for WORK or IWORK.
!           IERR = 2 => Method failed to converge in ITMAX steps.
!           IERR = 3 => Error in user input.
!                       Check input values of N, ITOL.
!           IERR = 4 => User error tolerance set too tight.
!                       Reset to 500*R1MACH(3).  Iteration proceeded.
!           IERR = 5 => Preconditioning matrix, M, is not positive
!                       definite.  (r,z) < 0.
!           IERR = 6 => Matrix A is not positive definite.  (p,Ap) < 0.
! IUNIT  :IN       Integer.
!         Unit number on which to write the error at each iteration,
!         if this is desired for monitoring convergence.  If unit
!         number is 0, no writing will occur.
! R      :WORK     Real R(N).
! Z      :WORK     Real Z(N).
! P      :WORK     Real P(N).
! ATP    :WORK     Real ATP(N).
! ATZ    :WORK     Real ATZ(N).
! DZ     :WORK     Real DZ(N).
! ATDZ   :WORK     Real ATDZ(N).
!         Real arrays used for workspace.
! RWORK  :WORK     Real RWORK(USER DEFINED).
!         Real array that can be used by  MSOLVE.
! IWORK  :WORK     Integer IWORK(USER DEFINED).
!         Integer array that can be used by  MSOLVE.
!
! *Description:
!       This  routine applies the  preconditioned conjugate gradient
!       (PCG) method to a non-symmetric system of equations Ax=b. To
!       do this the normal equations are solved:
!               AA' y  = b, where  x  = A'y.
!       In PCG method the iteration count is determined by condition
!                               -1
!       number of the  matrix (M  A).   In the  situation where  the
!       normal equations are  used  to solve a  non-symmetric system
!       the condition number depends on  AA' and should therefore be
!       much worse than that of A.  This is the conventional wisdom.
!       When one has a good preconditioner for AA' this may not hold.
!       The latter is the situation when SCGN should be tried.
!
!       If one is trying to solve  a symmetric system, SCG should be
!       used instead.
!
!       This routine does  not care  what matrix data   structure is
!       used for A and M.  It simply calls MATVEC, MTTVEC and MSOLVE
!       routines, with arguments as described above.  The user could
!       write any type of structure, and  appropriate MATVEC, MTTVEC
!       and MSOLVE routines.  It is assumed  that A is stored in the
!       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is
!       stored  in  IWORK  and  RWORK)  in  some fashion.   The SLAP
!       routines SSDCGN and SSLUCN are examples of this procedure.
!
!       Two  examples  of  matrix  data structures  are the: 1) SLAP
!       Triad  format and 2) SLAP Column format.
!
!       =================== S L A P Triad format ===================
!
!       In  this   format only the  non-zeros are  stored.  They may
!       appear  in *ANY* order.   The user  supplies three arrays of
!       length NELT, where  NELT  is the number  of non-zeros in the
!       matrix:  (IA(NELT), JA(NELT),  A(NELT)).  For each  non-zero
!       the  user puts   the row  and  column index   of that matrix
!       element in the IA and JA arrays.  The  value of the non-zero
!       matrix  element is  placed in  the corresponding location of
!       the A  array.  This is  an extremely easy data  structure to
!       generate.  On  the other hand it  is  not too  efficient  on
!       vector  computers   for the  iterative  solution  of  linear
!       systems.  Hence, SLAP  changes this input  data structure to
!       the SLAP   Column  format for the  iteration (but   does not
!       change it back).
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
!       In  this format   the non-zeros are    stored counting  down
!       columns (except  for the diagonal  entry, which must  appear
!       first in each "column") and are  stored in the real array A.
!       In other words,  for  each column    in the matrix   put the
!       diagonal  entry  in A.   Then   put  in the  other  non-zero
!       elements going   down the  column (except  the  diagonal) in
!       order.  The IA array holds the row index  for each non-zero.
!       The JA array holds the offsets into the IA, A arrays for the
!       beginning   of   each  column.      That is,   IA(JA(ICOL)),
!       A(JA(ICOL)) points to the beginning of the ICOL-th column in
!       IA and  A.  IA(JA(ICOL+1)-1), A(JA(ICOL+1)-1)  points to the
!       end of the ICOL-th column.  Note that we always have JA(N+1)
!       = NELT+1, where N is the number of columns in the matrix and
!       NELT is the number of non-zeros in the matrix.
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
!***SEE ALSO  SSDCGN, SSLUCN, ISSCGN
!***REFERENCES  1. Mark K. Seager, A SLAP for the Masses, in
!                  G. F. Carey, Ed., Parallel Supercomputing: Methods,
!                  Algorithms and Applications, Wiley, 1989, pp.135-155.
!***ROUTINES CALLED  ISSCGN, R1MACH, SAXPY, SCOPY, SDOT
!***REVISION HISTORY  (YYMMDD)
!   871119  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890921  Removed TeX from comments.  (FNF)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   891004  Added new reference.
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Removed MATVEC, MTTVEC and MSOLVE from ROUTINES CALLED
!           list.  (FNF)
!   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
!   920511  Added complete declaration section.  (WRB)
!   920929  Corrected format of reference.  (FNF)
!   921019  Changed 500.0 to 500 to reduce SP/DP differences.  (FNF)
!   921113  Corrected C***CATEGORY line.  (FNF)
!***END PROLOGUE  SCGN
!     .. Scalar Arguments ..
  REAL ERR, TOL
  INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, N, NELT
!     .. Array Arguments ..
  REAL A(N), ATDZ(N), ATP(N), ATZ(N), B(N), DZ(N), P(N), R(N), &
       RWORK(*), X(N), Z(N)
  INTEGER IA(NELT), IWORK(*), JA(NELT)
!     .. Subroutine Arguments ..
  EXTERNAL MATVEC, MSOLVE, MTTVEC
!     .. Local Scalars ..
  REAL AK, AKDEN, BK, BKDEN, BKNUM, BNRM, SOLNRM, TOLMIN
  INTEGER I, K
!     .. External Functions ..
  REAL R1MACH, SDOT
  INTEGER ISSCGN
  EXTERNAL R1MACH, SDOT, ISSCGN
!     .. External Subroutines ..
  EXTERNAL SAXPY, SCOPY
!***FIRST EXECUTABLE STATEMENT  SCGN
!
!         Check user input.
!
  ITER = 0
  IERR = 0
  if (  N < 1 ) THEN
     IERR = 3
     return
  end if
  TOLMIN = 500*R1MACH(3)
  if (  TOL < TOLMIN ) THEN
     TOL = TOLMIN
     IERR = 4
  end if
!         Calculate initial residual and pseudo-residual, and check
!         stopping criterion.
  call MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
  DO 10 I = 1, N
     R(I) = B(I) - R(I)
 10   CONTINUE
  call MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
  call MTTVEC(N, Z, ATZ, NELT, IA, JA, A, ISYM)
!
  if (  ISSCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC, MSOLVE, &
       ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, ATP, ATZ, &
       DZ, ATDZ, RWORK, IWORK, AK, BK, BNRM, SOLNRM)  /=  0 ) &
       go to 200
  if (  IERR /= 0 ) RETURN
!
!         ***** iteration loop *****
!
  DO 100 K=1,ITMAX
     ITER = K
!
!         Calculate coefficient BK and direction vector P.
     BKNUM = SDOT(N, Z, 1, R, 1)
     if (  BKNUM <= 0.0E0 ) THEN
        IERR = 6
        return
     ENDIF
     if ( ITER  ==  1) THEN
        call SCOPY(N, Z, 1, P, 1)
     ELSE
        BK = BKNUM/BKDEN
        DO 20 I = 1, N
           P(I) = Z(I) + BK*P(I)
 20         CONTINUE
     ENDIF
     BKDEN = BKNUM
!
!         Calculate coefficient AK, new iterate X, new residual R,
!         and new pseudo-residual ATZ.
     if ( ITER  /=  1) call SAXPY(N, BK, ATP, 1, ATZ, 1)
     call SCOPY(N, ATZ, 1, ATP, 1)
     AKDEN = SDOT(N, ATP, 1, ATP, 1)
     if (  AKDEN <= 0.0E0 ) THEN
        IERR = 6
        return
     ENDIF
     AK = BKNUM/AKDEN
     call SAXPY(N, AK, ATP, 1, X, 1)
     call MATVEC(N, ATP, Z, NELT, IA, JA, A, ISYM)
     call SAXPY(N, -AK, Z, 1, R, 1)
     call MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
     call MTTVEC(N, Z, ATZ, NELT, IA, JA, A, ISYM)
!
!         check stopping criterion.
     if (  ISSCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC, &
          MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, &
          Z, P, ATP, ATZ, DZ, ATDZ, RWORK, IWORK, AK, BK, BNRM, &
          SOLNRM)  /=  0) GOTO 200
!
 100  CONTINUE
!
!         *****   end of loop  *****
!
!         stopping criterion not satisfied.
  ITER = ITMAX + 1
!
 200  return
!------------- LAST LINE OF SCGN FOLLOWS ----------------------------
end
