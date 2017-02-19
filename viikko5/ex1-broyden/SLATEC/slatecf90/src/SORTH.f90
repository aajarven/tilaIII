subroutine SORTH (VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
!
!! SORTH is an internal routine for SGMRES.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A4, D2B4
!***TYPE      SINGLE PRECISION (SORTH-S, DORTH-D)
!***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
!             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
!***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
!           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
!           Seager, Mark K., (LLNL), seager@llnl.gov
!             Lawrence Livermore National Laboratory
!             PO Box 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!***DESCRIPTION
!        This routine  orthogonalizes  the  vector  VNEW  against the
!        previous KMP  vectors in the   V array.  It uses  a modified
!        Gram-Schmidt   orthogonalization procedure with  conditional
!        reorthogonalization.
!
! *Usage:
!      INTEGER N, LL, LDHES, KMP
!      REAL VNEW(N), V(N,LL), HES(LDHES,LL), SNORMW
!
!      call SORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
!
! *Arguments:
! VNEW   :INOUT    Real VNEW(N)
!         On input, the vector of length N containing a scaled
!         product of the Jacobian and the vector V(*,LL).
!         On output, the new vector orthogonal to V(*,i0) to V(*,LL),
!         where i0 = max(1, LL-KMP+1).
! V      :IN       Real V(N,LL)
!         The N x LL array containing the previous LL
!         orthogonal vectors V(*,1) to V(*,LL).
! HES    :INOUT    Real HES(LDHES,LL)
!         On input, an LL x LL upper Hessenberg matrix containing,
!         in HES(I,K), K.lt.LL, the scaled inner products of
!         A*V(*,K) and V(*,i).
!         On return, column LL of HES is filled in with
!         the scaled inner products of A*V(*,LL) and V(*,i).
! N      :IN       Integer
!         The order of the matrix A, and the length of VNEW.
! LL     :IN       Integer
!         The current order of the matrix HES.
! LDHES  :IN       Integer
!         The leading dimension of the HES array.
! KMP    :IN       Integer
!         The number of previous vectors the new vector VNEW
!         must be made orthogonal to (KMP .le. MAXL).
! SNORMW :OUT      REAL
!         Scalar containing the l-2 norm of VNEW.
!
!***SEE ALSO  SGMRES
!***ROUTINES CALLED  SAXPY, SDOT, SNRM2
!***REVISION HISTORY  (YYMMDD)
!   871001  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910506  Made subsidiary to SGMRES.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!***END PROLOGUE  SORTH
!         The following is for optimized compilation on LLNL/LTSS Crays.
!LLL. OPTIMIZE
!     .. Scalar Arguments ..
  REAL SNORMW
  INTEGER KMP, LDHES, LL, N
!     .. Array Arguments ..
  REAL HES(LDHES,*), V(N,*), VNEW(*)
!     .. Local Scalars ..
  REAL ARG, SUMDSQ, TEM, VNRM
  INTEGER I, I0
!     .. External Functions ..
  REAL SDOT, SNRM2
  EXTERNAL SDOT, SNRM2
!     .. External Subroutines ..
  EXTERNAL SAXPY
!     .. Intrinsic Functions ..
  INTRINSIC MAX, SQRT
!***FIRST EXECUTABLE STATEMENT  SORTH
!
!         Get norm of unaltered VNEW for later use.
!
  VNRM = SNRM2(N, VNEW, 1)
!   -------------------------------------------------------------------
!         Perform the modified Gram-Schmidt procedure on VNEW =A*V(LL).
!         Scaled inner products give new column of HES.
!         Projections of earlier vectors are subtracted from VNEW.
!   -------------------------------------------------------------------
  I0 = MAX(1,LL-KMP+1)
  DO 10 I = I0,LL
     HES(I,LL) = SDOT(N, V(1,I), 1, VNEW, 1)
     TEM = -HES(I,LL)
     call SAXPY(N, TEM, V(1,I), 1, VNEW, 1)
 10   CONTINUE
!   -------------------------------------------------------------------
!         Compute SNORMW = norm of VNEW.  If VNEW is small compared
!         to its input value (in norm), then reorthogonalize VNEW to
!         V(*,1) through V(*,LL).  Correct if relative correction
!         exceeds 1000*(unit roundoff).  Finally, correct SNORMW using
!         the dot products involved.
!   -------------------------------------------------------------------
  SNORMW = SNRM2(N, VNEW, 1)
  if (VNRM + 0.001E0*SNORMW  /=  VNRM) RETURN
  SUMDSQ = 0
  DO 30 I = I0,LL
     TEM = -SDOT(N, V(1,I), 1, VNEW, 1)
     if (HES(I,LL) + 0.001E0*TEM  ==  HES(I,LL)) go to 30
     HES(I,LL) = HES(I,LL) - TEM
     call SAXPY(N, TEM, V(1,I), 1, VNEW, 1)
     SUMDSQ = SUMDSQ + TEM**2
 30   CONTINUE
  if (SUMDSQ  ==  0.0E0) RETURN
  ARG = MAX(0.0E0,SNORMW**2 - SUMDSQ)
  SNORMW = SQRT(ARG)
!
  return
!------------- LAST LINE OF SORTH FOLLOWS ----------------------------
end
