subroutine CGEEV (A, LDA, N, E, V, LDV, WORK, JOB, INFO)
!
!! CGEEV computes the eigenvalues and, optionally, the eigenvectors ...
!            of a complex general matrix.
!
!***LIBRARY   SLATEC
!***CATEGORY  D4A4
!***TYPE      COMPLEX (SGEEV-S, CGEEV-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, GENERAL MATRIX
!***AUTHOR  Kahaner, D. K., (NBS)
!           Moler, C. B., (U. of New Mexico)
!           Stewart, G. W., (U. of Maryland)
!***DESCRIPTION
!
!     Abstract
!      CGEEV computes the eigenvalues and, optionally,
!      the eigenvectors of a general complex matrix.
!
!     Call Sequence Parameters-
!       (The values of parameters marked with * (star) will be changed
!         by CGEEV.)
!
!        A*      COMPLEX(LDA,N)
!                complex nonsymmetric input matrix.
!
!        LDA     INTEGER
!                set by the user to
!                the leading dimension of the complex array A.
!
!        N       INTEGER
!                set by the user to
!                the order of the matrices A and V, and
!                the number of elements in E.
!
!        E*      COMPLEX(N)
!                on return from CGEEV E contains the eigenvalues of A.
!                See also INFO below.
!
!        V*      COMPLEX(LDV,N)
!                on return from CGEEV if the user has set JOB
!                = 0        V is not referenced.
!                = nonzero  the N eigenvectors of A are stored in the
!                first N columns of V.  See also INFO below.
!                (If the input matrix A is nearly degenerate, V
!                 will be badly conditioned, i.e. have nearly
!                 dependent columns.)
!
!        LDV     INTEGER
!                set by the user to
!                the leading dimension of the array V if JOB is also
!                set nonzero.  In that case N must be  <=  LDV.
!                If JOB is set to zero LDV is not referenced.
!
!        WORK*   REAL(3N)
!                temporary storage vector.  Contents changed by CGEEV.
!
!        JOB     INTEGER
!                set by the user to
!                = 0        eigenvalues only to be calculated by CGEEV.
!                           neither V nor LDV are referenced.
!                = nonzero  eigenvalues and vectors to be calculated.
!                           In this case A & V must be distinct arrays.
!                           Also,  if LDA > LDV,  CGEEV changes all the
!                           elements of A thru column N.  If LDA < LDV,
!                           CGEEV changes all the elements of V through
!                           column N.  If LDA = LDV only A(I,J) and V(I,
!                           J) for I,J = 1,...,N are changed by CGEEV.
!
!        INFO*   INTEGER
!                on return from CGEEV the value of INFO is
!                = 0  normal return, calculation successful.
!                = K  if the eigenvalue iteration fails to converge,
!                     eigenvalues K+1 through N are correct, but
!                     no eigenvectors were computed even if they were
!                     requested (JOB nonzero).
!
!      Error Messages
!           No. 1  recoverable  N is greater than LDA
!           No. 2  recoverable  N is less than one.
!           No. 3  recoverable  JOB is nonzero and N is greater than LDV
!           No. 4  warning      LDA > LDV,  elements of A other than the
!                               N by N input elements have been changed
!           No. 5  warning      LDA < LDV,  elements of V other than the
!                               N by N output elements have been changed
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CBABK2, CBAL, COMQR, COMQR2, CORTH, SCOPY, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800808  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  CGEEV
  INTEGER I,IHI,ILO,INFO,J,K,L,LDA,LDV,MDIM,N
  REAL A(*),E(*),WORK(*),V(*)
!***FIRST EXECUTABLE STATEMENT  CGEEV
  if (N  >  LDA) call XERMSG ('SLATEC', 'CGEEV', 'N  >  LDA.', 1, &
     1)
  if ( N  >  LDA) RETURN
  if (N  <  1) call XERMSG ('SLATEC', 'CGEEV', 'N  <  1', 2, 1)
  if ( N  <  1) RETURN
  if ( N  ==  1 .AND. JOB  ==  0) go to 35
  MDIM = 2 * LDA
  if ( JOB  ==  0) go to 5
  if (N  >  LDV) call XERMSG ('SLATEC', 'CGEEV', &
     'JOB  /=  0 AND N  >  LDV.', 3, 1)
  if ( N  >  LDV) RETURN
  if ( N  ==  1) go to 35
!
!       REARRANGE A if NECESSARY WHEN LDA > LDV AND JOB  /= 0
!
  MDIM = MIN(MDIM,2 * LDV)
  if (LDA  <  LDV) call XERMSG ('SLATEC', 'CGEEV', &
     'LDA < LDV,  ELEMENTS OF V OTHER THAN THE N BY N OUTPUT ' // &
     'ELEMENTS HAVE BEEN CHANGED.', 5, 0)
  if ( LDA <= LDV) go to 5
  call XERMSG ('SLATEC', 'CGEEV', &
     'LDA > LDV, ELEMENTS OF A OTHER THAN THE N BY N INPUT ' // &
     'ELEMENTS HAVE BEEN CHANGED.', 4, 0)
  L = N - 1
  DO 4 J=1,L
      I = 2 * N
     M = 1+J*2*LDV
     K = 1+J*2*LDA
     call SCOPY(I,A(K),1,A(M),1)
    4 CONTINUE
    5 CONTINUE
!
!     SEPARATE REAL AND IMAGINARY PARTS
!
  DO 6 J = 1, N
   K = (J-1) * MDIM +1
   L = K + N
   call SCOPY(N,A(K+1),2,WORK(1),1)
   call SCOPY(N,A(K),2,A(K),1)
   call SCOPY(N,WORK(1),1,A(L),1)
    6 CONTINUE
!
!     SCALE AND ORTHOGONAL REDUCTION TO HESSENBERG.
!
  call CBAL(MDIM,N,A(1),A(N+1),ILO,IHI,WORK(1))
  call CORTH(MDIM,N,ILO,IHI,A(1),A(N+1),WORK(N+1),WORK(2*N+1))
  if ( JOB  /=  0) go to 10
!
!     EIGENVALUES ONLY
!
  call COMQR(MDIM,N,ILO,IHI,A(1),A(N+1),E(1),E(N+1),INFO)
  go to 30
!
!     EIGENVALUES AND EIGENVECTORS.
!
   10 call COMQR2(MDIM,N,ILO,IHI,WORK(N+1),WORK(2*N+1),A(1),A(N+1), &
    E(1),E(N+1),V(1),V(N+1),INFO)
  if (INFO  /=  0) go to 30
  call CBABK2(MDIM,N,ILO,IHI,WORK(1),N,V(1),V(N+1))
!
!     CONVERT EIGENVECTORS TO COMPLEX STORAGE.
!
  DO 20 J = 1,N
   K = (J-1) * MDIM + 1
   I = (J-1) * 2 * LDV + 1
   L = K + N
   call SCOPY(N,V(K),1,WORK(1),1)
   call SCOPY(N,V(L),1,V(I+1),2)
   call SCOPY(N,WORK(1),1,V(I),2)
   20 CONTINUE
!
!     CONVERT EIGENVALUES TO COMPLEX STORAGE.
!
   30 call SCOPY(N,E(1),1,WORK(1),1)
  call SCOPY(N,E(N+1),1,E(2),2)
  call SCOPY(N,WORK(1),1,E(1),2)
  return
!
!     TAKE CARE OF N=1 CASE
!
   35 E(1) = A(1)
  E(2) = A(2)
  INFO = 0
  if ( JOB  ==  0) RETURN
  V(1) = A(1)
  V(2) = A(2)
  return
end
