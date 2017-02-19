  INTEGER FUNCTION ISDIR (N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, &
     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, DZ, RWORK, &
     IWORK, BNRM, SOLNRM)
!
!! ISDIR is the Preconditioned Iterative Refinement Stop Test.
!
!            This routine calculates the stop test for the iterative
!            refinement iteration scheme.  It returns a non-zero if the
!            error estimate (the type of which is determined by ITOL)
!            is less than the user specified tolerance TOL.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A4, D2B4
!***TYPE      DOUBLE PRECISION (ISSIR-S, ISDIR-D)
!***KEYWORDS  LINEAR SYSTEM, SLAP, SPARSE, STOP TEST
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!
! *Usage:
!     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
!     INTEGER IERR, IUNIT, IWORK(USER DEFINED)
!     DOUBLE PRECISION B(N), X(N), A(N), TOL, ERR, R(N), Z(N), DZ(N)
!     DOUBLE PRECISION RWORK(USER DEFINED), BNRM, SOLNRM
!     EXTERNAL MSOLVE
!
!     if (  ISDIR(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,
!    $     ITMAX, ITER, ERR, IERR, IUNIT, R, Z, DZ, RWORK, IWORK,
!    $     BNRM, SOLNRM)  /=  0 ) THEN ITERATION DONE
!
! *Arguments:
! N      :IN       Integer.
!         Order of the Matrix.
! B      :IN       Double Precision B(N).
!         Right-hand side vector.
! X      :IN       Double Precision X(N).
!         The current approximate solution vector.
! NELT   :IN       Integer.
!         Number of Non-Zeros stored in A.
! IA     :IN       Integer IA(NELT).
! JA     :IN       Integer JA(NELT).
! A      :IN       Double Precision A(NELT).
!         These arrays contain the matrix data structure for A.
!         It could take any form.  See "C *Description" in the
!         DIR routine.
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all non-zero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the upper
!         or lower triangle of the matrix is stored.
! MSOLVE :EXT      External.
!         Name of a routine which solves a linear system Mz = r for
!         z given r with the preconditioning matrix M (M is supplied via
!         RWORK and IWORK arrays.  The name of the MSOLVE routine must
!         be declared external in the calling program.  The calling
!         sequence to MSOLVE is:
!             call MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!         Where N is the number of unknowns, R is the right-hand side
!         vector and Z is the solution upon return.  NELT, IA, JA, A and
!         ISYM are defined as above.  RWORK is a double precision array
!         that can be used to pass necessary preconditioning information
!         and/or workspace to MSOLVE.  IWORK is an integer work array
!         for the same purpose as RWORK.
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
!             COMMON /DSLBLK/ SOLN( )
!         If ITOL=11, iteration stops when the 2-norm of the difference
!         between the iterative approximation and the user-supplied
!         solution divided by the 2-norm of the user-supplied solution
!         is less than TOL.  Note that this requires the user to set up
!         the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.
!         The routine with this declaration should be loaded before the
!         stop test so that the correct length is used by the loader.
!         This procedure is not standard Fortran and may not work
!         correctly on your system (although it has worked on every
!         system the authors have tried).  If ITOL is not 11 then this
!         common block is indeed standard Fortran.
! TOL    :IN       Double Precision.
!         Convergence criterion, as described above.
! ITMAX  :IN       Integer.
!         Maximum number of iterations.
! ITER   :IN       Integer.
!         Current iteration count.  (Must be zero on first call.)
! ERR    :OUT      Double Precision.
!         Error estimate of error in the X(N) approximate solution, as
!         defined by ITOL.
! IERR   :OUT      Integer.
!         Error flag.  IERR is set to 3 if ITOL is not one of the
!         acceptable values, see above.
! IUNIT  :IN       Integer.
!         Unit number on which to write the error at each iteration,
!         if this is desired for monitoring convergence.  If unit
!         number is 0, no writing will occur.
! R      :IN       Double Precision R(N).
!         The residual R = B-AX.
! Z      :WORK     Double Precision Z(N).
!         Workspace used to hold the pseudo-residual M z = r.
! DZ     :WORK     Double Precision DZ(N).
!         Workspace used to hold temporary vector(s).
! RWORK  :WORK     Double Precision RWORK(USER DEFINED).
!         Double Precision array that can be used by  MSOLVE.
! IWORK  :WORK     Integer IWORK(USER DEFINED).
!         Integer array that can be used by MSOLVE.
! BNRM   :INOUT    Double Precision.
!         Norm of the right hand side.  Type of norm depends on ITOL.
!         Calculated only on the first call.
! SOLNRM :INOUT    Double Precision.
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
!***SEE ALSO  DIR, DSJAC, DSGS
!***ROUTINES CALLED  D1MACH, DNRM2
!***COMMON BLOCKS    DSLBLK
!***REVISION HISTORY  (YYMMDD)
!   871119  DATE WRITTEN
!   880320  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   891003  Removed C***REFER TO line, per MKS.
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
!   910506  Made subsidiary to DIR.  (FNF)
!   920407  COMMON BLOCK renamed DSLBLK.  (WRB)
!   920511  Added complete declaration section.  (WRB)
!   921026  Changed 1.0E10 to D1MACH(2) and corrected E to D in
!           output format.  (FNF)
!***END PROLOGUE  ISDIR
!     .. Scalar Arguments ..
  DOUBLE PRECISION BNRM, ERR, SOLNRM, TOL
  INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, N, NELT
!     .. Array Arguments ..
  DOUBLE PRECISION A(NELT), B(N), DZ(N), R(N), RWORK(*), X(N), Z(N)
  INTEGER IA(NELT), IWORK(*), JA(NELT)
!     .. Subroutine Arguments ..
  EXTERNAL MSOLVE
!     .. Arrays in Common ..
  DOUBLE PRECISION SOLN(1)
!     .. Local Scalars ..
  INTEGER I
!     .. External Functions ..
  DOUBLE PRECISION D1MACH, DNRM2
  EXTERNAL D1MACH, DNRM2
!     .. Common blocks ..
  COMMON /DSLBLK/ SOLN
!***FIRST EXECUTABLE STATEMENT  ISDIR
  ISDIR = 0
  if (  ITOL == 1 ) THEN
!         err = ||Residual||/||RightHandSide|| (2-Norms).
     if ( ITER  ==  0) BNRM = DNRM2(N, B, 1)
     ERR = DNRM2(N, R, 1)/BNRM
  ELSE if (  ITOL == 2 ) THEN
!                  -1              -1
!         err = ||M  Residual||/||M  RightHandSide|| (2-Norms).
     if ( ITER  ==  0) THEN
        call MSOLVE(N, B, DZ, NELT, IA, JA, A, ISYM, RWORK, IWORK)
        BNRM = DNRM2(N, DZ, 1)
     ENDIF
     ERR = DNRM2(N, Z, 1)/BNRM
  ELSE if (  ITOL == 11 ) THEN
!         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms).
     if (  ITER == 0 ) SOLNRM = DNRM2(N, SOLN, 1)
     DO 10 I = 1, N
        DZ(I) = X(I) - SOLN(I)
 10      CONTINUE
     ERR = DNRM2(N, DZ, 1)/SOLNRM
  ELSE
!
!         If we get here ITOL is not one of the acceptable values.
     ERR = D1MACH(2)
     IERR = 3
  end if
!
  if (  IUNIT /= 0 ) THEN
     WRITE(IUNIT,1000) ITER,ERR
  end if
!
  if (  ERR <= TOL ) ISDIR = 1
!
  return
 1000 FORMAT(5X,'ITER = ',I4,' Error Estimate = ',D16.7)
!------------- LAST LINE OF ISDIR FOLLOWS -----------------------------
end
