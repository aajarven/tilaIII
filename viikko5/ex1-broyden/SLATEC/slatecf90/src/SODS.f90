subroutine SODS (A, X, B, NEQ, NUK, NRDA, IFLAG, WORK, IWORK)
!
!! SODS solves an overdetermined linear system for BVSUP.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to BVSUP
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (SODS-S)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     SODS solves the overdetermined system of linear equations A X = B,
!     where A is NEQ by NUK and NEQ  >=  NUK. If rank A = NUK,
!     X is the UNIQUE least squares solution vector. That is,
!              R(1)**2 + ..... + R(NEQ)**2 = minimum
!     where R is the residual vector  R = B - A X.
!     If rank A  <  NUK , the least squares solution of minimal
!     length can be provided.
!     SODS is an interfacing routine which calls subroutine LSSODS
!     for the solution. LSSODS in turn calls subroutine ORTHOL and
!     possibly subroutine OHTROR for the decomposition of A by
!     orthogonal transformations. In the process, ORTHOL calls upon
!     subroutine CSCALE for scaling.
!
! **********************************************************************
!   Input
! **********************************************************************
!
!     A -- Contains the matrix of NEQ equations in NUK unknowns and must
!          be dimensioned NRDA by NUK. The original A is destroyed
!     X -- Solution array of length at least NUK
!     B -- Given constant vector of length NEQ, B is destroyed
!     NEQ -- Number of equations, NEQ greater or equal to 1
!     NUK -- Number of columns in the matrix (which is also the number
!            of unknowns), NUK not larger than NEQ
!     NRDA -- Row dimension of A, NRDA greater or equal to NEQ
!     IFLAG -- Status indicator
!            =0 For the first call (and for each new problem defined by
!               a new matrix A) when the matrix data is treated as exact
!           =-K For the first call (and for each new problem defined by
!               a new matrix A) when the matrix data is assumed to be
!               accurate to about K digits
!            =1 For subsequent calls whenever the matrix A has already
!               been decomposed (problems with new vectors B but
!               same matrix a can be handled efficiently)
!     WORK(*),IWORK(*) -- Arrays for storage of internal information,
!                     WORK must be dimensioned at least  2 + 5*NUK
!                     IWORK must be dimensioned at least NUK+2
!     IWORK(2) -- Scaling indicator
!                 =-1 If the matrix A is to be pre-scaled by
!                 columns when appropriate
!                 If the scaling indicator is not equal to -1
!                 no scaling will be attempted
!              For most problems scaling will probably not be necessary
!
! **********************************************************************
!   OUTPUT
! **********************************************************************
!
!     IFLAG -- Status indicator
!            =1 If solution was obtained
!            =2 If improper input is detected
!            =3 If rank of matrix is less than NUK
!               If the minimal length least squares solution is
!               desired, simply reset IFLAG=1 and call the code again
!     X -- Least squares solution of  A X = B
!     A -- Contains the strictly upper triangular part of the reduced
!           matrix and the transformation information
!     WORK(*),IWORK(*) -- Contains information needed on subsequent
!                         Calls (IFLAG=1 case on input) which must not
!                         be altered
!                         WORK(1) contains the Euclidean norm of
!                         the residual vector
!                         WORK(2) contains the Euclidean norm of
!                         the solution vector
!                         IWORK(1) contains the numerically determined
!                         rank of the matrix A
!
! **********************************************************************
!
!***SEE ALSO  BVSUP
!***REFERENCES  G. Golub, Numerical methods for solving linear least
!                 squares problems, Numerische Mathematik 7, (1965),
!                 pp. 206-216.
!               P. Businger and G. Golub, Linear least squares
!                 solutions by Householder transformations, Numerische
!                 Mathematik  7, (1965), pp. 269-276.
!               H. A. Watts, Solving linear least squares problems
!                 using SODS/SUDS/CODS, Sandia Report SAND77-0683,
!                 Sandia Laboratories, 1977.
!***ROUTINES CALLED  LSSODS
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SODS
  DIMENSION A(NRDA,*),X(*),B(*),WORK(*),IWORK(*)
!
!***FIRST EXECUTABLE STATEMENT  SODS
  ITER=0
  IS=2
  IP=3
  KS=2
  KD=3
  KZ=KD+NUK
  KV=KZ+NUK
  KT=KV+NUK
  KC=KT+NUK
!
  call LSSODS(A,X,B,NEQ,NUK,NRDA,IFLAG,IWORK(1),IWORK(IS),A, &
              WORK(KD),IWORK(IP),ITER,WORK(1),WORK(KS), &
              WORK(KZ),B,WORK(KV),WORK(KT),WORK(KC))
!
  return
end
