subroutine DTIN (N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB)
!
!! DTIN reads in SLAP Triad Format Linear System.
!            Routine to read in a SLAP Triad format matrix and right
!            hand side and solution to the system, if known.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  N1
!***TYPE      DOUBLE PRECISION (STIN-S, DTIN-D)
!***KEYWORDS  DIAGNOSTICS, LINEAR SYSTEM, SLAP SPARSE
!***AUTHOR  Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!
! *Usage:
!     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT, JOB
!     DOUBLE PRECISION A(NELT), SOLN(N), RHS(N)
!
!     call DTIN( N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB )
!
! *Arguments:
! N      :OUT      Integer
!         Order of the Matrix.
! NELT   :INOUT    Integer.
!         On input NELT is the maximum number of non-zeros that
!         can be stored in the IA, JA, A arrays.
!         On output NELT is the number of non-zeros stored in A.
! IA     :OUT      Integer IA(NELT).
! JA     :OUT      Integer JA(NELT).
! A      :OUT      Double Precision A(NELT).
!         On output these arrays hold the matrix A in the SLAP
!         Triad format.  See "Description", below.
! ISYM   :OUT      Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all non-zero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the lower
!         triangle of the matrix is stored.
! SOLN   :OUT      Double Precision SOLN(N).
!         The solution to the linear system, if present.  This array
!         is accessed if and only if JOB to read it in, see below.
!         If the user requests that SOLN be read in, but it is not in
!         the file, then it is simply zeroed out.
! RHS    :OUT      Double Precision RHS(N).
!         The right hand side vector.  This array is accessed if and
!         only if JOB is set to read it in, see below.
!         If the user requests that RHS be read in, but it is not in
!         the file, then it is simply zeroed out.
! IUNIT  :IN       Integer.
!         Fortran logical I/O device unit number to write the matrix
!         to.  This unit must be connected in a system dependent fashion
!         to a file or the console or you will get a nasty message
!         from the Fortran I/O libraries.
! JOB    :INOUT    Integer.
!         Flag indicating what I/O operations to perform.
!         On input JOB indicates what Input operations to try to
!         perform.
!         JOB = 0 => Read only the matrix.
!         JOB = 1 => Read matrix and RHS (if present).
!         JOB = 2 => Read matrix and SOLN (if present).
!         JOB = 3 => Read matrix, RHS and SOLN (if present).
!         On output JOB indicates what operations were actually
!         performed.
!         JOB = 0 => Read in only the matrix.
!         JOB = 1 => Read in the matrix and RHS.
!         JOB = 2 => Read in the matrix and SOLN.
!         JOB = 3 => Read in the matrix, RHS and SOLN.
!
! *Description:
!       The format for the  input is as follows.  On  the first line
!       are counters and flags: N, NELT, ISYM, IRHS, ISOLN.  N, NELT
!       and ISYM are described above.  IRHS is  a flag indicating if
!       the RHS was  written out (1 is  yes, 0 is  no).  ISOLN  is a
!       flag indicating if the SOLN was written out  (1 is yes, 0 is
!       no).  The format for the fist line is: 5i10.  Then comes the
!       NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format
!       for  these lines is   :  1X,I5,1X,I5,1X,D16.7.   Then  comes
!       RHS(I), I = 1, N, if IRHS = 1.  Then  comes SOLN(I), I  = 1,
!       N, if ISOLN = 1.  The format for these lines is: 1X,D16.7.
!
!       =================== S L A P Triad format ===================
!       This routine requires that the  matrix A be   stored in  the
!       SLAP  Triad format.  In  this format only the non-zeros  are
!       stored.  They may appear in  *ANY* order.  The user supplies
!       three arrays of  length NELT, where  NELT is  the number  of
!       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
!       each non-zero the user puts the row and column index of that
!       matrix element  in the IA and  JA arrays.  The  value of the
!       non-zero  matrix  element is  placed   in  the corresponding
!       location of the A array.   This is  an  extremely  easy data
!       structure to generate.  On  the  other hand it   is  not too
!       efficient on vector computers for  the iterative solution of
!       linear systems.  Hence,   SLAP changes   this  input    data
!       structure to the SLAP Column format  for  the iteration (but
!       does not change it back).
!
!       Here is an example of the  SLAP Triad   storage format for a
!       5x5 Matrix.  Recall that the entries may appear in any order.
!
!           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.
!                              1  2  3  4  5  6  7  8  91011
!       |1112  0  015|   A: 5112113315535522354421
!       |2122  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
!       | 0  033  035|  JA:  1  2  1  3  5  3  5  2  5  4  1
!       | 0  0  044  0|
!       |51  053  055|
!
! *Cautions:
!     This routine will attempt to write to the Fortran logical output
!     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
!     this logical unit is attached to a file or terminal before calling
!     this routine with a non-zero value for IUNIT.  This routine does
!     not check for the validity of a non-zero IUNIT unit number.
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
!   921007  Changed E's to D's in formats.  (FNF)
!   930701  Updated CATEGORY section.  (FNF, WRB)
!***END PROLOGUE  DTIN
!     .. Scalar Arguments ..
  INTEGER ISYM, IUNIT, JOB, N, NELT
!     .. Array Arguments ..
  DOUBLE PRECISION A(NELT), RHS(N), SOLN(N)
  INTEGER IA(NELT), JA(NELT)
!     .. Local Scalars ..
  INTEGER I, IRHS, ISOLN, JOBRET, NELTMX
!     .. Intrinsic Functions ..
  INTRINSIC MIN
!***FIRST EXECUTABLE STATEMENT  DTIN
!
!         Read in the information heading.
!
  NELTMX = NELT
  READ(IUNIT,1000) N, NELT, ISYM, IRHS, ISOLN
  NELT = MIN( NELT, NELTMX )
!
!         Read in the matrix non-zeros in Triad format.
  DO 10 I = 1, NELT
     READ(IUNIT,1010) IA(I), JA(I), A(I)
 10   CONTINUE
!
!         If requested, read in the rhs.
  JOBRET = 0
  if (  JOB == 1 .OR. JOB == 3 ) THEN
!
!         Check to see if rhs is in the file.
     if (  IRHS == 1 ) THEN
        JOBRET = 1
        READ(IUNIT,1020) (RHS(I),I=1,N)
     ELSE
        DO 20 I = 1, N
           RHS(I) = 0
 20         CONTINUE
     ENDIF
  end if
!
!         If requested, read in the solution.
  if (  JOB > 1 ) THEN
!
!         Check to see if solution is in the file.
     if (  ISOLN == 1 ) THEN
        JOBRET = JOBRET + 2
        READ(IUNIT,1020) (SOLN(I),I=1,N)
     ELSE
        DO 30 I = 1, N
           SOLN(I) = 0
 30         CONTINUE
     ENDIF
  end if
!
  JOB = JOBRET
  return
 1000 FORMAT(5I10)
 1010 FORMAT(1X,I5,1X,I5,1X,D16.7)
 1020 FORMAT(1X,D16.7)
!------------- LAST LINE OF DTIN FOLLOWS ----------------------------
end
