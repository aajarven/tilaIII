subroutine SSIEV (A, LDA, N, E, WORK, JOB, INFO)
!
!! SSIEV computes the eigenvalues and eigenvectors of a real symmetric matrix.
!
!***LIBRARY   SLATEC
!***CATEGORY  D4A1
!***TYPE      SINGLE PRECISION (SSIEV-S, CHIEV-C)
!***KEYWORDS  COMPLEX HERMITIAN, EIGENVALUES, EIGENVECTORS, MATRIX,
!             SYMMETRIC
!***AUTHOR  Kahaner, D. K., (NBS)
!           Moler, C. B., (U. of New Mexico)
!           Stewart, G. W., (U. of Maryland)
!***DESCRIPTION
!
!     Abstract
!      SSIEV computes the eigenvalues and, optionally, the eigenvectors
!      of a real symmetric matrix.
!
!     Call Sequence Parameters-
!       (The values of parameters marked with * (star) will be  changed
!         by SSIEV.)
!
!       A*      REAL (LDA,N)
!               real symmetric input matrix.
!               Only the diagonal and upper triangle of A must be input,
!               as SSIEV copies the upper triangle to the lower.
!               That is, the user must define A(I,J), I=1,..N, and J=I,.
!               ..,N.
!               On return from SSIEV, if the user has set JOB
!               = 0        the lower triangle of A has been altered.
!               = nonzero  the N eigenvectors of A are stored in its
!               first N columns.  See also INFO below.
!
!       LDA     INTEGER
!               set by the user to
!               the leading dimension of the array A.
!
!       N       INTEGER
!               set by the user to
!               the order of the matrix A and
!               the number of elements in E.
!
!       E*      REAL (N)
!               on return from SSIEV, E contains the N
!               eigenvalues of A.  See also INFO below.
!
!       WORK*   REAL (2*N)
!               temporary storage vector.  Contents changed by SSIEV.
!
!       JOB     INTEGER
!               set by user on input
!               = 0         only calculate eigenvalues of A.
!               = nonzero   calculate eigenvalues and eigenvectors of A.
!
!       INFO*   INTEGER
!               on return from SSIEV, the value of INFO is
!               = 0 for normal return.
!               = K if the eigenvalue iteration fails to converge.
!                   eigenvalues and vectors 1 through K-1 are correct.
!
!
!     Error Messages-
!          No. 1   recoverable  N is greater than LDA
!          No. 2   recoverable  N is less than one
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  IMTQL2, TQLRAT, TRED1, TRED2, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800808  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  SSIEV
  INTEGER INFO,JOB,LDA,N
  REAL A(LDA,*),E(*),WORK(*)
!***FIRST EXECUTABLE STATEMENT  SSIEV
   if (N  >  LDA) call XERMSG ('SLATEC', 'SSIEV', 'N  >  LDA.', &
         1, 1)
   if ( N  >  LDA) RETURN
   if (N  <  1) call XERMSG ('SLATEC', 'SSIEV', 'N  <  1', 2, 1)
   if ( N  <  1) RETURN
!
!       CHECK N=1 CASE
!
  E(1) = A(1,1)
  INFO = 0
  if ( N  ==  1) RETURN
!
!     COPY UPPER TRIANGLE TO LOWER
!
  DO 10 J=1,N
  DO 10 I=1,J
     A(J,I)=A(I,J)
   10 CONTINUE
!
  if ( JOB /= 0) go to 20
!
!     EIGENVALUES ONLY
!
  call TRED1(LDA,N,A,E,WORK(1),WORK(N+1))
  call TQLRAT(N,E,WORK(N+1),INFO)
  return
!
!     EIGENVALUES AND EIGENVECTORS
!
   20 call TRED2(LDA,N,A,E,WORK,A)
  call IMTQL2(LDA,N,E,WORK,A,INFO)
  return
end
