subroutine SSICS (N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL, D, &
     R, IWARN)
!
!! SSICS is the Incomplete Cholesky Decomposition Preconditioner SLAP Set Up.
!
!            Routine to generate the Incomplete Cholesky decomposition,
!            L*D*L-trans, of a symmetric positive definite matrix, A,
!            which is stored in SLAP Column format.  The unit lower
!            triangular matrix L is stored by rows, and the inverse of
!            the diagonal matrix D is stored.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2E
!***TYPE      SINGLE PRECISION (SSICS-S, DSICS-D)
!***KEYWORDS  INCOMPLETE CHOLESKY FACTORIZATION,
!             ITERATIVE PRECONDITION, LINEAR SYSTEM, SLAP SPARSE
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
!     INTEGER NEL, IEL(NEL), JEL(NEL), IWARN
!     REAL    A(NELT), EL(NEL), D(N), R(N)
!
!     call SSICS( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL, D, R,
!    $    IWARN )
!
! *Arguments:
! N      :IN       Integer.
!         Order of the Matrix.
! NELT   :IN       Integer.
!         Number of elements in arrays IA, JA, and A.
! IA     :INOUT    Integer IA(NELT).
! JA     :INOUT    Integer JA(NELT).
! A      :INOUT    Real A(NELT).
!         These arrays should hold the matrix A in the SLAP Column
!         format.  See "Description", below.
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all non-zero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the lower
!         triangle of the matrix is stored.
! NEL    :OUT      Integer.
!         Number of non-zeros in the lower triangle of A.   Also
!         corresponds to the length of the IEL, JEL, EL arrays.
! IEL    :OUT      Integer IEL(NEL).
! JEL    :OUT      Integer JEL(NEL).
! EL     :OUT      Real     EL(NEL).
!         IEL, JEL, EL contain the unit lower triangular factor  of the
!         incomplete decomposition   of the A  matrix  stored  in  SLAP
!         Row format.   The Diagonal of   ones   *IS*   stored.     See
!         "Description", below for more details about the SLAP Row fmt.
! D      :OUT      Real D(N)
!         Upon return this array holds D(I) = 1./DIAG(A).
! R      :WORK     Real R(N).
!         Temporary real workspace needed for the factorization.
! IWARN  :OUT      Integer.
!         This is a warning variable and is zero if the IC factoriza-
!         tion goes well.  It is set to the row index corresponding to
!         the last zero pivot found.  See "Description", below.
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
!       ==================== S L A P Row format ====================
!
!       This routine requires  that the matrix A  be  stored  in the
!       SLAP  Row format.   In this format  the non-zeros are stored
!       counting across  rows (except for the diagonal  entry, which
!       must appear first in each "row") and  are stored in the real
!       array A.  In other words, for each row in the matrix put the
!       diagonal entry in  A.   Then   put  in the   other  non-zero
!       elements   going  across the  row (except   the diagonal) in
!       order.   The  JA array  holds   the column   index for  each
!       non-zero.   The IA  array holds the  offsets into  the JA, A
!       arrays  for   the   beginning  of   each  row.   That    is,
!       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
!       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
!       points to the  end of the  IROW-th row.  Note that we always
!       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
!       the matrix  and NELT  is the  number   of  non-zeros in  the
!       matrix.
!
!       Here is an example of the SLAP Row storage format for a  5x5
!       Matrix (in the A and JA arrays '|' denotes the end of a row):
!
!           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
!                              1  2  3    4  5    6  7    8    91011
!       |1112  0  015|   A: 111215 | 2221 | 3335 | 44 | 555153
!       |2122  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
!       | 0  033  035|  IA:  1  4  6    8  9   12
!       | 0  0  044  0|
!       |51  053  055|
!
!       With the SLAP  format some  of  the   "inner  loops" of this
!       routine should vectorize  on  machines with hardware support
!       for vector   gather/scatter  operations.  Your compiler  may
!       require a compiler directive to  convince it that  there are
!       no  implicit  vector  dependencies.  Compiler directives for
!       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
!       supplied with the standard SLAP distribution.
!
!       The IC factorization does not always exist for SPD matrices.
!       In the event that a zero pivot is found it is set  to be 1.0
!       and the factorization proceeds.   The integer variable IWARN
!       is set to the last row where the Diagonal was fudged.  This
!       eventuality hardly ever occurs in practice.
!
!***SEE ALSO  SCG, SSICCG
!***REFERENCES  1. Gene Golub and Charles Van Loan, Matrix Computations,
!                  Johns Hopkins University Press, Baltimore, Maryland,
!                  1983.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   871119  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   900805  Changed XERRWV calls to calls to XERMSG.  (RWC)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   920511  Added complete declaration section.  (WRB)
!   920929  Corrected format of reference.  (FNF)
!   930701  Updated CATEGORY section.  (FNF, WRB)
!***END PROLOGUE  SSICS
!     .. Scalar Arguments ..
  INTEGER ISYM, IWARN, N, NEL, NELT
!     .. Array Arguments ..
  REAL A(NELT), D(N), EL(NEL), R(N)
  INTEGER IA(NELT), IEL(NEL), JA(NELT), JEL(NEL)
!     .. Local Scalars ..
  REAL ELTMP
  INTEGER I, IBGN, IC, ICBGN, ICEND, ICOL, IEND, IR, IRBGN, IREND, &
          IROW, IRR, J, JBGN, JELTMP, JEND
  CHARACTER XERN1*8
!     .. External Subroutines ..
  EXTERNAL XERMSG
!***FIRST EXECUTABLE STATEMENT  SSICS
!
!         Set the lower triangle in IEL, JEL, EL
!
  IWARN = 0
!
!         All matrix elements stored in IA, JA, A.  Pick out the lower
!         triangle (making sure that the Diagonal of EL is one) and
!         store by rows.
!
  NEL = 1
  IEL(1) = 1
  JEL(1) = 1
  EL(1) = 1
  D(1) = A(1)
!VD$R NOCONCUR
  DO 30 IROW = 2, N
!         Put in the Diagonal.
     NEL = NEL + 1
     IEL(IROW) = NEL
     JEL(NEL) = IROW
     EL(NEL) = 1
     D(IROW) = A(JA(IROW))
!
!         Look in all the lower triangle columns for a matching row.
!         Since the matrix is symmetric, we can look across the
!         ITOW-th row by looking down the IROW-th column (if it is
!         stored ISYM=0)...
     if (  ISYM == 0 ) THEN
        ICBGN = JA(IROW)
        ICEND = JA(IROW+1)-1
     ELSE
        ICBGN = 1
        ICEND = IROW-1
     ENDIF
     DO 20 IC = ICBGN, ICEND
        if (  ISYM == 0 ) THEN
           ICOL = IA(IC)
           if (  ICOL >= IROW ) GOTO 20
        ELSE
           ICOL = IC
        ENDIF
        JBGN = JA(ICOL)+1
        JEND = JA(ICOL+1)-1
        if (  JBGN <= JEND .AND. IA(JEND) >= IROW ) THEN
!VD$ NOVECTOR
           DO 10 J = JBGN, JEND
              if (  IA(J) == IROW ) THEN
                 NEL = NEL + 1
                 JEL(NEL) = ICOL
                 EL(NEL)  = A(J)
                 GOTO 20
              ENDIF
 10            CONTINUE
        ENDIF
 20      CONTINUE
 30   CONTINUE
  IEL(N+1) = NEL+1
!
!         Sort ROWS of lower triangle into descending order (count out
!         along rows out from Diagonal).
!
  DO 60 IROW = 2, N
     IBGN = IEL(IROW)+1
     IEND = IEL(IROW+1)-1
     if (  IBGN < IEND ) THEN
        DO 50 I = IBGN, IEND-1
!VD$ NOVECTOR
           DO 40 J = I+1, IEND
              if (  JEL(I) > JEL(J) ) THEN
                 JELTMP = JEL(J)
                 JEL(J) = JEL(I)
                 JEL(I) = JELTMP
                 ELTMP = EL(J)
                 EL(J) = EL(I)
                 EL(I) = ELTMP
              ENDIF
 40            CONTINUE
 50         CONTINUE
     ENDIF
 60   CONTINUE
!
!         Perform the Incomplete Cholesky decomposition by looping
!         over the rows.
!         Scale the first column.  Use the structure of A to pick out
!         the rows with something in column 1.
!
  IRBGN = JA(1)+1
  IREND = JA(2)-1
  DO 65 IRR = IRBGN, IREND
     IR = IA(IRR)
!         Find the index into EL for EL(1,IR).
!         Hint: it's the second entry.
     I = IEL(IR)+1
     EL(I) = EL(I)/D(1)
 65   CONTINUE
!
  DO 110 IROW = 2, N
!
!         Update the IROW-th diagonal.
!
     DO 66 I = 1, IROW-1
        R(I) = 0
 66      CONTINUE
     IBGN = IEL(IROW)+1
     IEND = IEL(IROW+1)-1
     if (  IBGN <= IEND ) THEN
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NODEPCHK
        DO 70 I = IBGN, IEND
           R(JEL(I)) = EL(I)*D(JEL(I))
           D(IROW) = D(IROW) - EL(I)*R(JEL(I))
 70         CONTINUE
!
!         Check to see if we have a problem with the diagonal.
!
        if (  D(IROW) <= 0.0E0 ) THEN
           if (  IWARN == 0 ) IWARN = IROW
           D(IROW) = 1
        ENDIF
     ENDIF
!
!         Update each EL(IROW+1:N,IROW), if there are any.
!         Use the structure of A to determine the Non-zero elements
!         of the IROW-th column of EL.
!
     IRBGN = JA(IROW)
     IREND = JA(IROW+1)-1
     DO 100 IRR = IRBGN, IREND
        IR = IA(IRR)
        if (  IR <= IROW ) GOTO 100
!         Find the index into EL for EL(IR,IROW)
        IBGN = IEL(IR)+1
        IEND = IEL(IR+1)-1
        if (  JEL(IBGN) > IROW ) GOTO 100
        DO 90 I = IBGN, IEND
           if (  JEL(I) == IROW ) THEN
              ICEND = IEND
 91               if (  JEL(ICEND) >= IROW ) THEN
                 ICEND = ICEND - 1
                 GOTO 91
              ENDIF
!         Sum up the EL(IR,1:IROW-1)*R(1:IROW-1) contributions.
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NODEPCHK
              DO 80 IC = IBGN, ICEND
                 EL(I) = EL(I) - EL(IC)*R(JEL(IC))
 80               CONTINUE
              EL(I) = EL(I)/D(IROW)
              GOTO 100
           ENDIF
 90         CONTINUE
!
!         If we get here, we have real problems...
        WRITE (XERN1, '(I8)') IROW
        call XERMSG ('SLATEC', 'SSICS', &
           'A and EL data structure mismatch in row '// XERN1, 1, 2)
 100     CONTINUE
 110  CONTINUE
!
!         Replace diagonals by their inverses.
!
  D(1:n) = 1.0E0/D(1:n)

  return
end
