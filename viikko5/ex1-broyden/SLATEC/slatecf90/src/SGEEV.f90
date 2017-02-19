subroutine SGEEV (A, LDA, N, E, V, LDV, WORK, JOB, INFO)
!
!! SGEEV computes the eigenvalues and, optionally, the eigenvectors ...
!            of a real general matrix.
!
!***LIBRARY   SLATEC
!***CATEGORY  D4A2
!***TYPE      SINGLE PRECISION (SGEEV-S, CGEEV-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, GENERAL MATRIX
!***AUTHOR  Kahaner, D. K., (NBS)
!           Moler, C. B., (U. of New Mexico)
!           Stewart, G. W., (U. of Maryland)
!***DESCRIPTION
!
!     Abstract
!      SGEEV computes the eigenvalues and, optionally,
!      the eigenvectors of a general real matrix.
!
!     Call Sequence Parameters-
!       (The values of parameters marked with * (star) will be changed
!         by SGEEV.)
!
!        A*      REAL(LDA,N)
!                real nonsymmetric input matrix.
!
!        LDA     INTEGER
!                set by the user to
!                the leading dimension of the real array A.
!
!        N       INTEGER
!                set by the user to
!                the order of the matrices A and V, and
!                the number of elements in E.
!
!        E*      COMPLEX(N)
!                on return from SGEEV, E contains the eigenvalues of A.
!                See also INFO below.
!
!        V*      COMPLEX(LDV,N)
!                on return from SGEEV, if the user has set JOB
!                = 0        V is not referenced.
!                = nonzero  the N eigenvectors of A are stored in the
!                first N columns of V.  See also INFO below.
!                (Note that if the input matrix A is nearly degenerate,
!                 V may be badly conditioned, i.e., may have nearly
!                 dependent columns.)
!
!        LDV     INTEGER
!                set by the user to
!                the leading dimension of the array V if JOB is also
!                set nonzero.  In that case, N must be  <=  LDV.
!                If JOB is set to zero, LDV is not referenced.
!
!        WORK*   REAL(2N)
!                temporary storage vector.  Contents changed by SGEEV.
!
!        JOB     INTEGER
!                set by the user to
!                = 0        eigenvalues only to be calculated by SGEEV.
!                           Neither V nor LDV is referenced.
!                = nonzero  eigenvalues and vectors to be calculated.
!                           In this case, A & V must be distinct arrays.
!                           Also, if LDA  >  LDV, SGEEV changes all the
!                           elements of A thru column N.  If LDA < LDV,
!                           SGEEV changes all the elements of V through
!                           column N. If LDA = LDV, only A(I,J) and V(I,
!                           J) for I,J = 1,...,N are changed by SGEEV.
!
!        INFO*   INTEGER
!                on return from SGEEV the value of INFO is
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
!           No. 4  warning      LDA > LDV, elements of A other than the
!                               N by N input elements have been changed.
!           No. 5  warning      LDA < LDV, elements of V other than the
!                               N x N output elements have been changed.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  BALANC, BALBAK, HQR, HQR2, ORTHES, ORTRAN, SCOPY,
!                    SCOPYM, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800808  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  SGEEV
  INTEGER I,IHI,ILO,INFO,J,JB,JOB,K,KM,KP,L,LDA,LDV, &
          MDIM,N
  REAL A(*),E(*),WORK(*),V(*)
!***FIRST EXECUTABLE STATEMENT  SGEEV
  if (N  >  LDA) call XERMSG ('SLATEC', 'SGEEV', 'N  >  LDA.', 1, &
     1)
  if (N  >  LDA) RETURN
  if (N  <  1) call XERMSG ('SLATEC', 'SGEEV', 'N  <  1', 2, 1)
  if ( N  <  1) RETURN
  if ( N  ==  1 .AND. JOB  ==  0) go to 35
  MDIM = LDA
  if ( JOB  ==  0) go to 5
  if (N  >  LDV) call XERMSG ('SLATEC', 'SGEEV', &
     'JOB  /=  0 AND N  >  LDV.', 3, 1)
  if ( N  >  LDV) RETURN
  if ( N  ==  1) go to 35
!
!       REARRANGE A if NECESSARY WHEN LDA > LDV AND JOB  /= 0
!
  MDIM = MIN(LDA,LDV)
  if (LDA  <  LDV) call XERMSG ('SLATEC', 'SGEEV', &
     'LDA < LDV,  ELEMENTS OF V OTHER THAN THE N BY N OUTPUT ' // &
     'ELEMENTS HAVE BEEN CHANGED.', 5, 0)
  if ( LDA <= LDV) go to 5
  call XERMSG ('SLATEC', 'SGEEV', &
     'LDA > LDV, ELEMENTS OF A OTHER THAN THE N BY N INPUT ' // &
     'ELEMENTS HAVE BEEN CHANGED.', 4, 0)
  L = N - 1
  DO 4 J=1,L
     M = 1+J*LDV
     K = 1+J*LDA
     call SCOPY(N,A(K),1,A(M),1)
    4 CONTINUE
    5 CONTINUE
!
!     SCALE AND ORTHOGONAL REDUCTION TO HESSENBERG.
!
  call BALANC(MDIM,N,A,ILO,IHI,WORK(1))
  call ORTHES(MDIM,N,ILO,IHI,A,WORK(N+1))
  if ( JOB  /=  0) go to 10
!
!     EIGENVALUES ONLY
!
  call HQR(LDA,N,ILO,IHI,A,E(1),E(N+1),INFO)
  go to 30
!
!     EIGENVALUES AND EIGENVECTORS.
!
   10 call ORTRAN(MDIM,N,ILO,IHI,A,WORK(N+1),V)
  call HQR2(MDIM,N,ILO,IHI,A,E(1),E(N+1),V,INFO)
  if (INFO  /=  0) go to 30
  call BALBAK(MDIM,N,ILO,IHI,WORK(1),N,V)
!
!     CONVERT EIGENVECTORS TO COMPLEX STORAGE.
!
  DO 20 JB = 1,N
     J=N+1-JB
     I=N+J
     K=(J-1)*MDIM+1
     KP=K+MDIM
     KM=K-MDIM
     if ( E(I) >= 0.0E0) call SCOPY(N,V(K),1,WORK(1),2)
     if ( E(I) < 0.0E0) call SCOPY(N,V(KM),1,WORK(1),2)
     if ( E(I) == 0.0E0) call sinit ( N,0.0E0,WORK(2),2)
     if ( E(I) > 0.0E0) call SCOPY(N,V(KP),1,WORK(2),2)
     if ( E(I) < 0.0E0) call SCOPYM(N,V(K),1,WORK(2),2)
     L=2*(J-1)*LDV+1
     call SCOPY(2*N,WORK(1),1,V(L),1)
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
  E(2) = 0.E0
  INFO = 0
  if ( JOB  ==  0) RETURN
  V(1) = A(1)
  V(2) = 0.E0
  return
end
