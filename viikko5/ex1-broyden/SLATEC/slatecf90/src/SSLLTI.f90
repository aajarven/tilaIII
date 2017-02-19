subroutine SSLLTI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!
!! SSLLTI is the SLAP MSOLVE for LDL' (IC) Factorization.
!
!            This routine acts as an interface between the SLAP generic
!            MSOLVE calling convention and the routine that actually
!                           -1
!            computes (LDL')  B = X.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2E
!***TYPE      SINGLE PRECISION (SSLLTI-S, DSLLTI-D)
!***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!       It is assumed that RWORK and IWORK have initialized with
!       the information required for SLLTI2:
!          IWORK(1) = NEL
!          IWORK(2) = Starting location of IEL in IWORK.
!          IWORK(3) = Starting location of JEL in IWORK.
!          IWORK(4) = Starting location of EL in RWORK.
!          IWORK(5) = Starting location of DINV in RWORK.
!       See the DESCRIPTION of SLLTI2 for details.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  SLLTI2
!***REVISION HISTORY  (YYMMDD)
!   871119  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Corrected conversion error.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!   921113  Corrected C***CATEGORY line.  (FNF)
!   930701  Updated CATEGORY section.  (FNF, WRB)
!***END PROLOGUE  SSLLTI
!     .. Scalar Arguments ..
  INTEGER ISYM, N, NELT
!     .. Array Arguments ..
  REAL A(NELT), B(*), RWORK(*), X(*)
  INTEGER IA(NELT), IWORK(*), JA(NELT)
!     .. Local Scalars ..
  INTEGER LOCDIN, LOCEL, LOCIEL, LOCJEL, NEL
!     .. External Subroutines ..
  EXTERNAL SLLTI2
!***FIRST EXECUTABLE STATEMENT  SSLLTI
  NEL = IWORK(1)
  LOCIEL = IWORK(3)
  LOCJEL = IWORK(2)
  LOCEL  = IWORK(4)
  LOCDIN = IWORK(5)
  call SLLTI2(N, B, X, NEL, IWORK(LOCIEL), IWORK(LOCJEL), &
       RWORK(LOCEL), RWORK(LOCDIN))
!
  return
!------------- LAST LINE OF SSLLTI FOLLOWS ----------------------------
end
