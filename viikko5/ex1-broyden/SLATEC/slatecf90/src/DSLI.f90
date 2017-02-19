subroutine DSLI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!
!! DSLI is the SLAP MSOLVE for Lower Triangle Matrix.
!
!            This routine acts as an interface between the SLAP generic
!            MSOLVE calling convention and the routine that actually
!            computes inv(L) * B = X.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A3
!***TYPE      DOUBLE PRECISION (SSLI-S, DSLI-D)
!***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!       It is assumed that RWORK and IWORK have initialized with
!       the information required for DSLI2:
!          IWORK(1) = NEL
!          IWORK(2) = Starting location of IEL in IWORK.
!          IWORK(3) = Starting location of JEL in IWORK.
!          IWORK(4) = Starting location of EL in RWORK.
!       See the DESCRIPTION of DSLI2 for details.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DSLI2
!***REVISION HISTORY  (YYMMDD)
!   871119  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   920511  Added complete declaration section.  (WRB)
!   921113  Corrected C***CATEGORY line.  (FNF)
!   930701  Updated CATEGORY section.  (FNF, WRB)
!***END PROLOGUE  DSLI
!     .. Scalar Arguments ..
  INTEGER ISYM, N, NELT
!     .. Array Arguments ..
  DOUBLE PRECISION A(NELT), B(N), RWORK(*), X(N)
  INTEGER IA(NELT), IWORK(10), JA(NELT)
!     .. Local Scalars ..
  INTEGER LOCEL, LOCIEL, LOCJEL, NEL
!     .. External Subroutines ..
  EXTERNAL DSLI2
!***FIRST EXECUTABLE STATEMENT  DSLI
!
  NEL = IWORK(1)
  LOCIEL = IWORK(2)
  LOCJEL = IWORK(3)
  LOCEL = IWORK(4)
  call DSLI2(N, B, X, NEL, IWORK(LOCIEL), IWORK(LOCJEL), &
       RWORK(LOCEL))
!
  return
!------------- LAST LINE OF DSLI FOLLOWS ----------------------------
end
