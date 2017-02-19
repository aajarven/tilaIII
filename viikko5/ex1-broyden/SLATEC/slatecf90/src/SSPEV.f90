subroutine SSPEV (A, N, E, V, LDV, WORK, JOB, INFO)
!
!! SSPEV computes eigenvalues and eigenvectors of a real symmetric matrix ...
!  stored in packed form.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A1
!***TYPE      SINGLE PRECISION (SSPEV-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK, PACKED, SYMMETRIC
!***AUTHOR  Kahaner, D. K., (NBS)
!           Moler, C. B., (U. of New Mexico)
!           Stewart, G. W., (U. of Maryland)
!***DESCRIPTION
!
!     Abstract
!      SSPEV computes the eigenvalues and, optionally, the eigenvectors
!      of a real symmetric matrix stored in packed form.
!
!     Call Sequence Parameters-
!       (The values of parameters marked with * (star) will be  changed
!         by SSPEV.)
!
!        A*      REAL(N*(N+1)/2)
!                real symmetric packed input matrix.  Contains upper
!                triangle and diagonal of A, by column (elements
!                11, 12, 22, 13, 23, 33, ...).
!
!        N       INTEGER
!                set by the user to
!                the order of the matrix A.
!
!        E*      REAL(N)
!                on return from SSPEV, E contains the eigenvalues of A.
!                See also INFO below.
!
!        V*      REAL(LDV,N)
!                on return from SSPEV, if the user has set JOB
!                = 0        V is not referenced.
!                = nonzero  the N eigenvectors of A are stored in the
!                first N columns of V.  See also INFO below.
!
!        LDV     INTEGER
!                set by the user to
!                the leading dimension of the array V if JOB is also
!                set nonzero.  In that case, N must be  <=  LDV.
!                If JOB is set to zero, LDV is not referenced.
!
!        WORK*   REAL(2N)
!                temporary storage vector.  Contents changed by SSPEV.
!
!        JOB     INTEGER
!                set by the user to
!                = 0        eigenvalues only to be calculated by SSPEV.
!                           Neither V nor LDV are referenced.
!                = nonzero  eigenvalues and vectors to be calculated.
!                           In this case, A & V must be distinct arrays.
!                           Also, if LDA  >  LDV, SSPEV changes all the
!                           elements of A thru column N.  If LDA < LDV,
!                           SSPEV changes all the elements of V through
!                           column N.  If LDA=LDV, only A(I,J) and V(I,
!                           J) for I,J = 1,...,N are changed by SSPEV.
!
!       INFO*   INTEGER
!               on return from SSPEV, the value of INFO is
!               = 0 for normal return.
!               = K if the eigenvalue iteration fails to converge.
!                   Eigenvalues and vectors 1 through K-1 are correct.
!
!
!     Error Messages-
!          No. 1   recoverable  N is greater than LDV and JOB is nonzero
!          No. 2   recoverable  N is less than one
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  IMTQL2, TQLRAT, TRBAK3, TRED3, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800808  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  SSPEV
  INTEGER I,INFO,J,LDV,M,N
  REAL A(*),E(*),V(LDV,*),WORK(*)
!***FIRST EXECUTABLE STATEMENT  SSPEV
   if (N  >  LDV) call XERMSG ('SLATEC', 'SSPEV', 'N  >  LDV.', &
      1, 1)
   if ( N  >  LDV) RETURN
   if (N  <  1) call XERMSG ('SLATEC', 'SSPEV', 'N  <  1', 2, 1)
   if ( N  <  1) RETURN
!
!       CHECK N=1 CASE
!
  E(1) = A(1)
  INFO = 0
  if ( N  ==  1) RETURN
!
  if ( JOB /= 0) go to 20
!
!     EIGENVALUES ONLY
!
  call TRED3(N,1,A,E,WORK(1),WORK(N+1))
  call TQLRAT(N,E,WORK(N+1),INFO)
  return
!
!     EIGENVALUES AND EIGENVECTORS
!
   20 call TRED3(N,1,A,E,WORK(1),WORK(1))
  DO 30 I = 1, N
    DO 25 J = 1, N
   25     V(I,J) = 0.
   30   V(I,I) = 1.
  call IMTQL2(LDV,N,E,WORK,V,INFO)
  M = N
  if ( INFO  /=  0) M = INFO - 1
  call TRBAK3(LDV,N,1,A,M,V)
  return
end
