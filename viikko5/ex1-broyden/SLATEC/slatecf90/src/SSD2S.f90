subroutine SSD2S (N, NELT, IA, JA, A, ISYM, DINV)
!
!! SSD2S is the Diagonal Scaling Preconditioner SLAP Normal Eqns Set Up.
!
!            Routine to compute the inverse of the diagonal of the
!            matrix A*A', where A is stored in SLAP-Column format.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2E
!***TYPE      SINGLE PRECISION (SSD2S-S, DSD2S-D)
!***KEYWORDS  DIAGONAL, SLAP SPARSE
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!
! *Usage:
!     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
!     REAL    A(NELT), DINV(N)
!
!     call SSD2S( N, NELT, IA, JA, A, ISYM, DINV )
!
! *Arguments:
! N      :IN       Integer
!         Order of the Matrix.
! NELT   :IN       Integer.
!         Number of elements in arrays IA, JA, and A.
! IA     :IN       Integer IA(NELT).
! JA     :IN       Integer JA(NELT).
! A      :IN       Real A(NELT).
!         These arrays should hold the matrix A in the SLAP Column
!         format.  See "Description", below.
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all non-zero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the upper
!         or lower triangle of the matrix is stored.
! DINV   :OUT      Real DINV(N).
!         Upon return this array holds 1./DIAG(A*A').
!
! *Description
!       =================== S L A P Column format ==================
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
!       With the SLAP  format  all  of  the   "inner  loops" of this
!       routine should vectorize  on  machines with hardware support
!       for vector   gather/scatter  operations.  Your compiler  may
!       require a compiler directive to  convince it that  there are
!       no  implicit  vector  dependencies.  Compiler directives for
!       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
!       supplied with the standard SLAP distribution.
!
!
! *Cautions:
!       This routine assumes that the diagonal of A is all  non-zero
!       and that the operation DINV = 1.0/DIAG(A*A') will not under-
!       flow or overflow. This is done so that the loop  vectorizes.
!       Matrices  with zero or near zero or very  large entries will
!       have numerical difficulties  and  must  be fixed before this
!       routine is called.
!
!***SEE ALSO  SSDCGN
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
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
!***END PROLOGUE  SSD2S
!     .. Scalar Arguments ..
  INTEGER ISYM, N, NELT
!     .. Array Arguments ..
  REAL A(NELT), DINV(N)
  INTEGER IA(NELT), JA(NELT)
!     .. Local Scalars ..
  INTEGER I, K, KBGN, KEND
!***FIRST EXECUTABLE STATEMENT  SSD2S
  DO 10 I = 1, N
     DINV(I) = 0
 10   CONTINUE
!
!         Loop over each column.
!VD$R NOCONCUR
  DO 40 I = 1, N
     KBGN = JA(I)
     KEND = JA(I+1) - 1
!
!         Add in the contributions for each row that has a non-zero
!         in this column.
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NODEPCHK
     DO 20 K = KBGN, KEND
        DINV(IA(K)) = DINV(IA(K)) + A(K)**2
 20      CONTINUE
     if (  ISYM == 1 ) THEN
!
!         Lower triangle stored by columns => upper triangle stored by
!         rows with Diagonal being the first entry.  Loop across the
!         rest of the row.
        KBGN = KBGN + 1
        if (  KBGN <= KEND ) THEN
           DO 30 K = KBGN, KEND
              DINV(I) = DINV(I) + A(K)**2
 30            CONTINUE
        ENDIF
     ENDIF
 40   CONTINUE
  DO 50 I=1,N
     DINV(I) = 1.0E0/DINV(I)
 50   CONTINUE
!
  return
!------------- LAST LINE OF SSD2S FOLLOWS ----------------------------
end
