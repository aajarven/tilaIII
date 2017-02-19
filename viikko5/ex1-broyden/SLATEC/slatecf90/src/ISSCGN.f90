  INTEGER FUNCTION ISSCGN (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, &
     MTTVEC, MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, &
     P, ATP, ATZ, DZ, ATDZ, RWORK, IWORK, AK, BK, BNRM, SOLNRM)
!
!! ISSCGN is the Preconditioned CG on Normal Equations Stop Test.
!
!            This routine calculates the stop test for the Conjugate
!            Gradient iteration scheme applied to the normal equations.
!            It returns a non-zero if the error estimate (the type of
!            which is determined by ITOL) is less than the user
!            specified tolerance TOL.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A4, D2B4
!***TYPE      SINGLE PRECISION (ISSCGN-S, ISDCGN-D)
!***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM,
!             NORMAL EQUATIONS, SLAP, SPARSE
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!
! *Usage:
!     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
!     INTEGER  IERR, IUNIT, IWORK(USER DEFINED)
!     REAL     B(N), X(N), A(N), TOL, ERR, R(N), Z(N), P(N)
!     REAL     ATP(N), ATZ(N), DZ(N), ATDZ(N)
!     REAL     RWORK(USER DEFINED), AK, BK, BNRM, SOLNRM
!     EXTERNAL MATVEC, MTTVEC, MSOLVE
!
!     if (  ISTPCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC,
!    $     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P,
!    $     ATP, ATZ, DZ, ATDZ, RWORK, IWORK, AK, BK, BNRM, SOLNRM)
!    $      /=  0 ) THEN ITERATION DONE
!
! *Arguments:
! N      :IN       Integer
!         Order of the Matrix.
! B      :IN       Real B(N).
!         Right-hand side vector.
! X      :IN       Real X(N).
!         The current approximate solution vector.
! NELT   :IN       Integer.
!         Number of Non-Zeros stored in A.
! IA     :IN       Integer IA(NELT).
! JA     :IN       Integer JA(NELT).
! A      :IN       Real A(NELT).
!         These arrays contain the matrix data structure for A.
!         It could take any form.  See "Description" in the
!         SCGN routine.
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
! TOL    :IN       Real.
!         Convergence criterion, as described above.
! ITMAX  :IN       Integer.
!         Maximum number of iterations.
! ITER   :IN       Integer.
!         Current iteration count.  (Must be zero on first call.)
! ERR    :OUT      Real.
!         Error estimate of error in the X(N) approximate solution, as
!         defined by ITOL.
! IERR   :OUT      Integer.
!         Error flag.  IERR is set to 3 if ITOL is not one of the
!         acceptable values, see above.
! IUNIT  :IN       Integer.
!         Unit number on which to write the error at each iteration,
!         if this is desired for monitoring convergence.  If unit
!         number is 0, no writing will occur.
! R      :IN       Real R(N).
!         The residual R = B-AX.
! Z      :WORK     Real Z(N).
!         Real array used for workspace.
! P      :IN       Real P(N).
!         The conjugate direction vector.
! ATP    :IN       Real ATP(N).
!         A-transpose times the conjugate direction vector.
! ATZ    :IN       Real ATZ(N).
!         A-transpose times the pseudo-residual.
! DZ     :IN       Real DZ(N).
!         Workspace used to hold temporary vector(s).
! ATDZ   :WORK       Real ATDZ(N).
!         Workspace.
! RWORK  :WORK     Real RWORK(USER DEFINED).
!         Real array that can be used by MSOLVE.
! IWORK  :WORK     Integer IWORK(USER DEFINED).
!         Integer array that can be used by MSOLVE.
! AK     :IN       Real.
! BK     :IN       Real.
!         Current conjugate gradient parameters alpha and beta.
! BNRM   :INOUT    Real.
!         Norm of the right hand side.  Type of norm depends on ITOL.
!         Calculated only on the first call.
! SOLNRM :INOUT    Real.
!         2-Norm of the true solution, SOLN.  Only computed and used
!         if ITOL = 11.
!
! *Function Return Values:
!       0 : Error estimate (determined by ITOL) is *NOT* less than the
!           specified tolerance, TOL.  The iteration must continue.
!       1 : Error estimate (determined by ITOL) is less than the
!           specified tolerance, TOL.  The iteration can be considered
!           complete.
!
! *Cautions:
!     This routine will attempt to write to the Fortran logical output
!     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
!     this logical unit is attached to a file or terminal before calling
!     this routine with a non-zero value for IUNIT.  This routine does
!     not check for the validity of a non-zero IUNIT unit number.
!
!***SEE ALSO  SCGN
!***ROUTINES CALLED  R1MACH, SNRM2
!***COMMON BLOCKS    SSLBLK
!***REVISION HISTORY  (YYMMDD)
!   871119  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   891003  Removed C***REFER TO line, per MKS.
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Removed MATVEC, MTTVEC and MSOLVE from ROUTINES CALLED
!           list.  (FNF)
!   910506  Made subsidiary to SCGN.  (FNF)
!   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
!   920511  Added complete declaration section.  (WRB)
!   920930  Corrected to not print AK,BK when ITER=0.  (FNF)
!   921026  Changed 1.0E10 to R1MACH(2).  (FNF)
!   921113  Corrected C***CATEGORY line.  (FNF)
!***END PROLOGUE  ISSCGN
!     .. Scalar Arguments ..
  REAL AK, BK, BNRM, ERR, SOLNRM, TOL
  INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, N, NELT
!     .. Array Arguments ..
  REAL A(N), ATDZ(N), ATP(N), ATZ(N), B(N), DZ(N), P(N), R(N), &
       RWORK(*), X(N), Z(N)
  INTEGER IA(NELT), IWORK(*), JA(NELT)
!     .. Subroutine Arguments ..
  EXTERNAL MATVEC, MSOLVE, MTTVEC
!     .. Arrays in Common ..
  REAL SOLN(1)
!     .. Local Scalars ..
  INTEGER I
!     .. External Functions ..
  REAL R1MACH, SNRM2
  EXTERNAL R1MACH, SNRM2
!     .. Common blocks ..
  COMMON /SSLBLK/ SOLN
!***FIRST EXECUTABLE STATEMENT  ISSCGN
  ISSCGN = 0
!
  if (  ITOL == 1 ) THEN
!         err = ||Residual||/||RightHandSide|| (2-Norms).
     if ( ITER  ==  0) BNRM = SNRM2(N, B, 1)
     ERR = SNRM2(N, R, 1)/BNRM
  ELSE if (  ITOL == 2 ) THEN
!                  -1              -1
!         err = ||M  Residual||/||M  RightHandSide|| (2-Norms).
     if ( ITER  ==  0) THEN
        call MSOLVE(N, B, DZ, NELT, IA, JA, A, ISYM, RWORK, IWORK)
        call MTTVEC(N, DZ, ATDZ, NELT, IA, JA, A, ISYM)
        BNRM = SNRM2(N, ATDZ, 1)
     ENDIF
     ERR = SNRM2(N, ATZ, 1)/BNRM
  ELSE if (  ITOL == 11 ) THEN
!         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms).
     if ( ITER  ==  0) SOLNRM = SNRM2(N, SOLN, 1)
     DO 10 I = 1, N
        DZ(I) = X(I) - SOLN(I)
 10      CONTINUE
     ERR = SNRM2(N, DZ, 1)/SOLNRM
  ELSE
!
!         If we get here ITOL is not one of the acceptable values.
     ERR = R1MACH(2)
     IERR = 3
  end if
!
  if (  IUNIT /= 0 ) THEN
     if (  ITER == 0 ) THEN
        WRITE(IUNIT,1000) N, ITOL
        WRITE(IUNIT,1010) ITER, ERR
     ELSE
        WRITE(IUNIT,1010) ITER, ERR, AK, BK
     ENDIF
  end if
  if (  ERR <= TOL ) ISSCGN = 1
!
  return
 1000 FORMAT(' PCG Applied to the Normal Equations for ', &
       'N, ITOL = ',I5, I5, &
       /' ITER','   Error Estimate','            Alpha', &
       '             Beta')
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7,1X,E16.7)
!------------- LAST LINE OF ISSCGN FOLLOWS ----------------------------
end
