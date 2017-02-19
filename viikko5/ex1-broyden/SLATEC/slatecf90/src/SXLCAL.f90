subroutine SXLCAL (N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM, &
     WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR, NELT, IA, JA, A, &
     ISYM)
!
!! SXLCAL is an internal routine for SGMRES.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A4, D2B4
!***TYPE      SINGLE PRECISION (SXLCAL-S, DXLCAL-D)
!***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
!             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
!***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
!           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
!           Seager, Mark K., (LLNL), seager@llnl.gov
!             Lawrence Livermore National Laboratory
!             PO Box 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!***DESCRIPTION
!        This  routine computes the solution  XL,  the current SGMRES
!        iterate, given the  V(I)'s and  the  QR factorization of the
!        Hessenberg  matrix HES.   This routine  is  only called when
!        ITOL=11.
!
! *Usage:
!      INTEGER N, LGMR, MAXLP1, JSCAL, JPRE, NMSL, IPAR(USER DEFINED)
!      INTEGER NELT, IA(NELT), JA(NELT), ISYM
!      REAL X(N), XL(N), ZL(N), HES(MAXLP1,MAXL), Q(2*MAXL),
!     $     V(N,MAXLP1), R0NRM, WK(N), SZ(N), RPAR(USER DEFINED),
!     $     A(NELT)
!      EXTERNAL MSOLVE
!
!      call SXLCAL(N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM,
!     $     WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR,
!     $     NELT, IA, JA, A, ISYM)
!
! *Arguments:
! N      :IN       Integer
!         The order of the matrix A, and the lengths
!         of the vectors SR, SZ, R0 and Z.
! LGMR   :IN       Integer
!         The number of iterations performed and
!         the current order of the upper Hessenberg
!         matrix HES.
! X      :IN       Real X(N)
!         The current approximate solution as of the last restart.
! XL     :OUT      Real XL(N)
!         An array of length N used to hold the approximate
!         solution X(L).
!         Warning: XL and ZL are the same array in the calling routine.
! ZL     :IN       Real ZL(N)
!         An array of length N used to hold the approximate
!         solution Z(L).
! HES    :IN       Real HES(MAXLP1,MAXL)
!         The upper triangular factor of the QR decomposition
!         of the (LGMR+1) by LGMR upper Hessenberg matrix whose
!         entries are the scaled inner-products of A*V(*,i) and V(*,k).
! MAXLP1 :IN       Integer
!         MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
!         MAXL is the maximum allowable order of the matrix HES.
! Q      :IN       Real Q(2*MAXL)
!         A real array of length 2*MAXL containing the components
!         of the Givens rotations used in the QR decomposition
!         of HES.  It is loaded in SHEQR.
! V      :IN       Real V(N,MAXLP1)
!         The N by(LGMR+1) array containing the LGMR
!         orthogonal vectors V(*,1) to V(*,LGMR).
! R0NRM  :IN       Real
!         The scaled norm of the initial residual for the
!         current call to SPIGMR.
! WK     :IN       Real WK(N)
!         A real work array of length N.
! SZ     :IN       Real SZ(N)
!         A vector of length N containing the non-zero
!         elements of the diagonal scaling matrix for Z.
! JSCAL  :IN       Integer
!         A flag indicating whether arrays SR and SZ are used.
!         JSCAL=0 means SR and SZ are not used and the
!                 algorithm will perform as if all
!                 SR(i) = 1 and SZ(i) = 1.
!         JSCAL=1 means only SZ is used, and the algorithm
!                 performs as if all SR(i) = 1.
!         JSCAL=2 means only SR is used, and the algorithm
!                 performs as if all SZ(i) = 1.
!         JSCAL=3 means both SR and SZ are used.
! JPRE   :IN       Integer
!         The preconditioner type flag.
! MSOLVE :EXT      External.
!         Name of the routine which solves a linear system Mz = r for
!         z given r with the preconditioning matrix M (M is supplied via
!         RPAR and IPAR arrays.  The name of the MSOLVE routine must
!         be declared external in the calling program.  The calling
!         sequence to MSOLVE is:
!             call MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
!         Where N is the number of unknowns, R is the right-hand side
!         vector and Z is the solution upon return.  NELT, IA, JA, A and
!         ISYM are defined as below.  RPAR is a real array that can be
!         used to pass necessary preconditioning information and/or
!         workspace to MSOLVE.  IPAR is an integer work array for the
!         same purpose as RPAR.
! NMSL   :IN       Integer
!         The number of calls to MSOLVE.
! RPAR   :IN       Real RPAR(USER DEFINED)
!         Real workspace passed directly to the MSOLVE routine.
! IPAR   :IN       Integer IPAR(USER DEFINED)
!         Integer workspace passed directly to the MSOLVE routine.
! NELT   :IN       Integer
!         The length of arrays IA, JA and A.
! IA     :IN       Integer IA(NELT)
!         An integer array of length NELT containing matrix data.
!         It is passed directly to the MATVEC and MSOLVE routines.
! JA     :IN       Integer JA(NELT)
!         An integer array of length NELT containing matrix data.
!         It is passed directly to the MATVEC and MSOLVE routines.
! A      :IN       Real A(NELT)
!         A real array of length NELT containing matrix data.
!         It is passed directly to the MATVEC and MSOLVE routines.
! ISYM   :IN       Integer
!         A flag to indicate symmetric matrix storage.
!         If ISYM=0, all non-zero entries of the matrix are
!         stored.  If ISYM=1, the matrix is symmetric and
!         only the upper or lower triangular part is stored.
!
!***SEE ALSO  SGMRES
!***ROUTINES CALLED  SAXPY, SCOPY, SHELS
!***REVISION HISTORY  (YYMMDD)
!   871001  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
!   910506  Made subsidiary to SGMRES.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!***END PROLOGUE  SXLCAL
!         The following is for optimized compilation on LLNL/LTSS Crays.
!LLL. OPTIMIZE
!     .. Scalar Arguments ..
  REAL R0NRM
  INTEGER ISYM, JPRE, JSCAL, LGMR, MAXLP1, N, NELT, NMSL
!     .. Array Arguments ..
  REAL A(NELT), HES(MAXLP1,*), Q(*), RPAR(*), SZ(*), V(N,*), WK(N), &
       X(N), XL(N), ZL(N)
  INTEGER IA(NELT), IPAR(*), JA(NELT)
!     .. Subroutine Arguments ..
  EXTERNAL MSOLVE
!     .. Local Scalars ..
  INTEGER I, K, LL, LLP1
!     .. External Subroutines ..
  EXTERNAL SAXPY, SCOPY, SHELS
!***FIRST EXECUTABLE STATEMENT  SXLCAL
  LL = LGMR
  LLP1 = LL + 1
  DO 10 K = 1,LLP1
     WK(K) = 0
 10   CONTINUE
  WK(1) = R0NRM
  call SHELS(HES, MAXLP1, LL, Q, WK)
  DO 20 K = 1,N
     ZL(K) = 0
 20   CONTINUE
  DO 30 I = 1,LL
     call SAXPY(N, WK(I), V(1,I), 1, ZL, 1)
 30   CONTINUE
  if ((JSCAL  ==  1) .OR.(JSCAL  ==  3)) THEN
     DO 40 K = 1,N
        ZL(K) = ZL(K)/SZ(K)
 40      CONTINUE
  end if

  if (JPRE  >  0) THEN
     call SCOPY(N, ZL, 1, WK, 1)
     call MSOLVE(N, WK, ZL, NELT, IA, JA, A, ISYM, RPAR, IPAR)
     NMSL = NMSL + 1
  end if
!
!  Calculate XL from X and ZL.
!
  XL(1:n) = X(1:n) + ZL(1:n)

  return
end
