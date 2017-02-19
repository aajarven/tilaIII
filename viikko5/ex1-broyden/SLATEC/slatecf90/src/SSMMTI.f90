subroutine SSMMTI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!
!! SSMMTI is the SLAP MSOLVE for LDU Factorization of Normal Equations.
!
!            This routine acts as an interface between the SLAP generic
!            MMTSLV calling convention and the routine that actually
!                                    -1
!            computes  [(LDU)*(LDU)']  B = X.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2E
!***TYPE      SINGLE PRECISION (SSMMTI-S, DSMMTI-D)
!***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!       It is assumed that RWORK and IWORK have initialized with
!       the information required for SSMMI2:
!          IWORK(1) = Starting location of IL in IWORK.
!          IWORK(2) = Starting location of JL in IWORK.
!          IWORK(3) = Starting location of IU in IWORK.
!          IWORK(4) = Starting location of JU in IWORK.
!          IWORK(5) = Starting location of L in RWORK.
!          IWORK(6) = Starting location of DINV in RWORK.
!          IWORK(7) = Starting location of U in RWORK.
!       See the DESCRIPTION of SSMMI2 for details.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  SSMMI2
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
!***END PROLOGUE  SSMMTI
!     .. Scalar Arguments ..
  INTEGER ISYM, N, NELT
!     .. Array Arguments ..
  REAL A(NELT), B(N), RWORK(*), X(N)
  INTEGER IA(NELT), IWORK(10), JA(NELT)
!     .. Local Scalars ..
  INTEGER LOCDIN, LOCIL, LOCIU, LOCJL, LOCJU, LOCL, LOCU
!     .. External Subroutines ..
  EXTERNAL SSMMI2
!***FIRST EXECUTABLE STATEMENT  SSMMTI
!
!         Pull out the locations of the arrays holding the ILU
!         factorization.
!
  LOCIL = IWORK(1)
  LOCJL = IWORK(2)
  LOCIU = IWORK(3)
  LOCJU = IWORK(4)
  LOCL = IWORK(5)
  LOCDIN = IWORK(6)
  LOCU = IWORK(7)
!
  call SSMMI2(N, B, X, IWORK(LOCIL), IWORK(LOCJL), &
       RWORK(LOCL), RWORK(LOCDIN), IWORK(LOCIU), &
       IWORK(LOCJU), RWORK(LOCU))
!
  return
!------------- LAST LINE OF SSMMTI FOLLOWS ----------------------------
end
