subroutine SHELS (A, LDA, N, Q, B)
!
!! SHELS is an internal routine for SGMRES.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A4, D2B4
!***TYPE      SINGLE PRECISION (SHELS-S, DHELS-D)
!***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
!             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
!***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
!           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
!           Seager, Mark K., (LLNL), seager@llnl.gov
!             Lawrence Livermore National Laboratory
!             PO Box 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!***DESCRIPTION
!        This routine is extracted from the LINPACK routine SGESL with
!        changes due to the fact that A is an upper Hessenberg matrix.
!
!        SHELS solves the least squares problem:
!
!                   MIN(B-A*X,B-A*X)
!
!        using the factors computed by SHEQR.
!
! *Usage:
!      INTEGER LDA, N
!      REAL A(LDA,N), Q(2*N), B(N+1)
!
!      call SHELS(A, LDA, N, Q, B)
!
! *Arguments:
! A       :IN       Real A(LDA,N)
!          The output from SHEQR which contains the upper
!          triangular factor R in the QR decomposition of A.
! LDA     :IN       Integer
!          The leading dimension of the array A.
! N       :IN       Integer
!          A is originally an (N+1) by N matrix.
! Q       :IN       Real Q(2*N)
!          The coefficients of the N Givens rotations
!          used in the QR factorization of A.
! B       :INOUT    Real B(N+1)
!          On input, B is the right hand side vector.
!          On output, B is the solution vector X.
!
!***SEE ALSO  SGMRES
!***ROUTINES CALLED  SAXPY
!***REVISION HISTORY  (YYMMDD)
!   871001  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Added C***FIRST EXECUTABLE STATEMENT line.  (FNF)
!   910506  Made subsidiary to SGMRES.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!***END PROLOGUE  SHELS
!         The following is for optimized compilation on LLNL/LTSS Crays.
!LLL. OPTIMIZE
!     .. Scalar Arguments ..
  INTEGER LDA, N
!     .. Array Arguments ..
  REAL A(LDA,*), B(*), Q(*)
!     .. Local Scalars ..
  REAL C, S, T, T1, T2
  INTEGER IQ, K, KB, KP1
!     .. External Subroutines ..
  EXTERNAL SAXPY
!***FIRST EXECUTABLE STATEMENT  SHELS
!
!         Minimize(B-A*X,B-A*X).  First form Q*B.
!
  DO K = 1, N
     KP1 = K + 1
     IQ = 2*(K-1) + 1
     C = Q(IQ)
     S = Q(IQ+1)
     T1 = B(K)
     T2 = B(KP1)
     B(K) = C*T1 - S*T2
     B(KP1) = S*T1 + C*T2
  end do
!
!  Now solve  R*X = Q*B.
!
  DO KB = 1, N
     K = N + 1 - KB
     B(K) = B(K)/A(K,K)
     T = -B(K)
     call SAXPY(K-1, T, A(1,K), 1, B(1), 1)
  end do

  return
end
