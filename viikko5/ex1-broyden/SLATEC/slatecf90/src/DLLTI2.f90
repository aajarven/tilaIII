subroutine DLLTI2 (N, B, X, NEL, IEL, JEL, EL, DINV)
!
!! DLLTI2 is the SLAP Backsolve routine for LDL' Factorization.
!
!            Routine to solve a system of the form  L*D*L' X = B,
!            where L is a unit lower triangular matrix and D is a
!            diagonal matrix and ' means transpose.
!
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2E
!***TYPE      DOUBLE PRECISION (SLLTI2-S, DLLTI2-D)
!***KEYWORDS  INCOMPLETE FACTORIZATION, ITERATIVE PRECONDITION, SLAP,
!             SPARSE, SYMMETRIC LINEAR SYSTEM SOLVE
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!
! *Usage:
!     INTEGER N, NEL, IEL(NEL), JEL(NEL)
!     DOUBLE PRECISION B(N), X(N), EL(NEL), DINV(N)
!
!     call DLLTI2( N, B, X, NEL, IEL, JEL, EL, DINV )
!
! *Arguments:
! N      :IN       Integer
!         Order of the Matrix.
! B      :IN       Double Precision B(N).
!         Right hand side vector.
! X      :OUT      Double Precision X(N).
!         Solution to L*D*L' x = b.
! NEL    :IN       Integer.
!         Number of non-zeros in the EL array.
! IEL    :IN       Integer IEL(NEL).
! JEL    :IN       Integer JEL(NEL).
! EL     :IN       Double Precision     EL(NEL).
!         IEL, JEL, EL contain the unit lower triangular factor   of
!         the incomplete decomposition   of the A  matrix  stored in
!         SLAP Row format.   The diagonal of ones *IS* stored.  This
!         structure can be set  up  by  the DS2LT routine.  See  the
!         "Description", below for more details about the  SLAP  Row
!         format.
! DINV   :IN       Double Precision DINV(N).
!         Inverse of the diagonal matrix D.
!
! *Description:
!       This routine is supplied with  the SLAP package as a routine
!       to perform the MSOLVE operation in the SCG iteration routine
!       for  the driver  routine DSICCG.   It must be called via the
!       SLAP  MSOLVE calling sequence  convention  interface routine
!       DSLLI.
!         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****
!               **** SLAP MSOLVE CALLING CONVENTION ****
!
!       IEL, JEL, EL should contain the unit lower triangular factor
!       of  the incomplete decomposition of  the A matrix  stored in
!       SLAP Row format.   This IC factorization  can be computed by
!       the  DSICS routine.  The  diagonal  (which is all one's) is
!       stored.
!
!       ==================== S L A P Row format ====================
!
!       This routine requires  that the matrix A  be  stored  in the
!       SLAP  Row format.   In this format  the non-zeros are stored
!       counting across  rows (except for the diagonal  entry, which
!       must  appear first  in each  "row")  and  are stored  in the
!       double precision  array A.  In other words, for each row  in
!       the matrix  put the diagonal  entry in A.   Then put in  the
!       other  non-zero elements  going across  the row  (except the
!       diagonal) in order.  The JA array holds the column index for
!       each non-zero.  The IA array holds the offsets  into the JA,
!       A  arrays  for  the   beginning  of  each  row.    That  is,
!       JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW-
!       th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
!       are  the last elements  of the  IROW-th row.   Note  that we
!       always have  IA(N+1) = NELT+1, where N is the number of rows
!       in the matrix  and  NELT is the  number of non-zeros  in the
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
!       With  the SLAP  Row format  the "inner loop" of this routine
!       should vectorize   on machines with   hardware  support  for
!       vector gather/scatter operations.  Your compiler may require
!       a  compiler directive  to  convince   it that there  are  no
!       implicit vector  dependencies.  Compiler directives  for the
!       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
!       with the standard SLAP distribution.
!
!***SEE ALSO  DSICCG, DSICS
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
!***END PROLOGUE  DLLTI2
!     .. Scalar Arguments ..
  INTEGER N, NEL
!     .. Array Arguments ..
  DOUBLE PRECISION B(N), DINV(N), EL(NEL), X(N)
  INTEGER IEL(NEL), JEL(NEL)
!     .. Local Scalars ..
  INTEGER I, IBGN, IEND, IROW
!***FIRST EXECUTABLE STATEMENT  DLLTI2
!
!         Solve  L*y = b,  storing result in x.
!
  DO 10 I=1,N
     X(I) = B(I)
 10   CONTINUE
  DO 30 IROW = 1, N
     IBGN = IEL(IROW) + 1
     IEND = IEL(IROW+1) - 1
     if (  IBGN <= IEND ) THEN
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NOCONCUR
!VD$ NODEPCHK
        DO 20 I = IBGN, IEND
           X(IROW) = X(IROW) - EL(I)*X(JEL(I))
 20         CONTINUE
     ENDIF
 30   CONTINUE
!
!         Solve  D*Z = Y,  storing result in X.
!
  DO 40 I=1,N
     X(I) = X(I)*DINV(I)
 40   CONTINUE
!
!         Solve  L-trans*X = Z.
!
  DO 60 IROW = N, 2, -1
     IBGN = IEL(IROW) + 1
     IEND = IEL(IROW+1) - 1
     if (  IBGN <= IEND ) THEN
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP
!VD$ NOCONCUR
!VD$ NODEPCHK
        DO 50 I = IBGN, IEND
           X(JEL(I)) = X(JEL(I)) - EL(I)*X(IROW)
 50         CONTINUE
     ENDIF
 60   CONTINUE
!
  return
!------------- LAST LINE OF DLLTI2 FOLLOWS ----------------------------
end
