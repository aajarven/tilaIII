subroutine SSDSCL (N, NELT, IA, JA, A, ISYM, X, B, DINV, JOB, ITOL)
!
!! SSDSCL is the Diagonal Scaling of system Ax = b.
!
!            This routine scales (and unscales) the system  Ax = b
!            by symmetric diagonal scaling.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2E
!***TYPE      SINGLE PRECISION (SSDSCL-S, DSDSCL-D)
!***KEYWORDS  DIAGONAL, SLAP SPARSE
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!
!    This routine scales (and unscales) the system Ax = b by symmetric
!    diagonal scaling.  The new system is:
!         -1/2  -1/2  1/2      -1/2
!        D    AD    (D   x) = D    b
!    when scaling is selected with the JOB parameter.  When unscaling
!    is selected this process is reversed.  The true solution is also
!    scaled or unscaled if ITOL is set appropriately, see below.
!
! *Usage:
!     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, JOB, ITOL
!     REAL    A(NELT), X(N), B(N), DINV(N)
!
!     call SSDSCL( N, NELT, IA, JA, A, ISYM, X, B, DINV, JOB, ITOL )
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
! X      :INOUT    Real X(N).
!         Initial guess that will be later used in the iterative
!         solution.
!         of the scaled system.
! B      :INOUT    Real B(N).
!         Right hand side vector.
! DINV   :INOUT    Real DINV(N).
!         Upon return this array holds 1./DIAG(A).
!         This is an input if JOB = 0.
! JOB    :IN       Integer.
!         Flag indicating whether to scale or not.
!         JOB non-zero means do scaling.
!         JOB = 0 means do unscaling.
! ITOL   :IN       Integer.
!         Flag indicating what type of error estimation to do in the
!         iterative method.  When ITOL = 11 the exact solution from
!         common block SSLBLK will be used.  When the system is scaled
!         then the true solution must also be scaled.  If ITOL is not
!         11 then this vector is not referenced.
!
! *Common Blocks:
! SOLN    :INOUT   Real SOLN(N).  COMMON BLOCK /SSLBLK/
!         The true solution, SOLN, is scaled (or unscaled) if ITOL is
!         set to 11, see above.
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
!       and that the operation DINV = 1.0/DIAG(A)  will  not  under-
!       flow or overflow. This is done so that the loop  vectorizes.
!       Matrices  with zero or near zero or very  large entries will
!       have numerical difficulties  and  must  be fixed before this
!       routine is called.
!
!***SEE ALSO  SSDCG
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    SSLBLK
!***REVISION HISTORY  (YYMMDD)
!   871119  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Added C***FIRST EXECUTABLE STATEMENT line.  (FNF)
!   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
!   920511  Added complete declaration section.  (WRB)
!   921113  Corrected C***CATEGORY line.  (FNF)
!   930701  Updated CATEGORY section.  (FNF, WRB)
!***END PROLOGUE  SSDSCL
!     .. Scalar Arguments ..
  INTEGER ISYM, ITOL, JOB, N, NELT
!     .. Array Arguments ..
  REAL A(NELT), B(N), DINV(N), X(N)
  INTEGER IA(NELT), JA(NELT)
!     .. Arrays in Common ..
  REAL SOLN(1)
!     .. Local Scalars ..
  REAL DI
  INTEGER ICOL, J, JBGN, JEND
!     .. Intrinsic Functions ..
  INTRINSIC SQRT
!     .. Common blocks ..
  COMMON /SSLBLK/ SOLN
!***FIRST EXECUTABLE STATEMENT  SSDSCL
!
!         SCALING...
!
  if (  JOB /= 0 ) THEN
     DO 10 ICOL = 1, N
        DINV(ICOL) = 1.0E0/SQRT( A(JA(ICOL)) )
 10      CONTINUE
  ELSE
!
!         UNSCALING...
!
     DO 15 ICOL = 1, N
        DINV(ICOL) = 1.0E0/DINV(ICOL)
 15      CONTINUE
  end if
!
  DO 30 ICOL = 1, N
     JBGN = JA(ICOL)
     JEND = JA(ICOL+1)-1
     DI = DINV(ICOL)
     DO 20 J = JBGN, JEND
        A(J) = DINV(IA(J))*A(J)*DI
 20      CONTINUE
 30   CONTINUE
!
  DO 40 ICOL = 1, N
     B(ICOL) = B(ICOL)*DINV(ICOL)
     X(ICOL) = X(ICOL)/DINV(ICOL)
 40   CONTINUE
!
!         Check to see if we need to scale the "true solution" as well.
!
  if (  ITOL == 11 ) THEN
     DO 50 ICOL = 1, N
        SOLN(ICOL) = SOLN(ICOL)/DINV(ICOL)
 50      CONTINUE
  end if
!
  return
end
