subroutine CHIEV (A, LDA, N, E, V, LDV, WORK, JOB, INFO)
!
!! CHIEV computes the eigenvalues and, optionally, the eigenvectors ...
!            of a complex Hermitian matrix.
!
!***LIBRARY   SLATEC
!***CATEGORY  D4A3
!***TYPE      COMPLEX (SSIEV-S, CHIEV-C)
!***KEYWORDS  COMPLEX HERMITIAN, EIGENVALUES, EIGENVECTORS, MATRIX,
!             SYMMETRIC
!***AUTHOR  Kahaner, D. K., (NBS)
!           Moler, C. B., (U. of New Mexico)
!           Stewart, G. W., (U. of Maryland)
!***DESCRIPTION
!
!     David Kahaner, Cleve Moler, G. W. Stewart,
!       N.B.S.         U.N.M.      N.B.S./U.MD.
!
!     Abstract
!      CHIEV computes the eigenvalues and, optionally,
!      the eigenvectors of a complex Hermitian matrix.
!
!     Call Sequence Parameters-
!       (the values of parameters marked with * (star) will be changed
!         by CHIEV.)
!
!        A*      COMPLEX(LDA,N)
!                complex Hermitian input matrix.
!                Only the upper triangle of A need be
!                filled in.  Elements on diagonal must be real.
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
!        E*      REAL(N)
!                on return from CHIEV E contains the eigenvalues of A.
!                See also INFO below.
!
!        V*      COMPLEX(LDV,N)
!                on return from CHIEV if the user has set JOB
!                = 0        V is not referenced.
!                = nonzero  the N eigenvectors of A are stored in the
!                first N columns of V.  See also INFO below.
!
!        LDV     INTEGER
!                set by the user to
!                the leading dimension of the array V if JOB is also
!                set nonzero.  In that case N must be  <=  LDV.
!                If JOB is set to zero LDV is not referenced.
!
!        WORK*   REAL(4N)
!                temporary storage vector.  Contents changed by CHIEV.
!
!        JOB     INTEGER
!                set by the user to
!                = 0        eigenvalues only to be calculated by CHIEV.
!                           Neither V nor LDV are referenced.
!                = nonzero  eigenvalues and vectors to be calculated.
!                           In this case A and V must be distinct arrays
!                           also if LDA  >  LDV CHIEV changes all the
!                           elements of A thru column N.  If LDA < LDV
!                           CHIEV changes all the elements of V through
!                           column N.  If LDA = LDV only A(I,J) and V(I,
!                           J) for I,J = 1,...,N are changed by CHIEV.
!
!        INFO*   INTEGER
!                on return from CHIEV the value of INFO is
!                = 0  normal return, calculation successful.
!                = K  if the eigenvalue iteration fails to converge,
!                     eigenvalues (and eigenvectors if requested)
!                     1 through K-1 are correct.
!
!      Error Messages
!           No. 1  recoverable  N is greater than LDA
!           No. 2  recoverable  N is less than one.
!           No. 3  recoverable  JOB is nonzero and N is greater than LDV
!           No. 4  warning      LDA > LDV,  elements of A other than the
!                               N by N input elements have been changed
!           No. 5  warning      LDA < LDV,  elements of V other than the
!                               N by N output elements have been changed
!           No. 6  recoverable  nonreal element on diagonal of A.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  HTRIBK, HTRIDI, IMTQL2, SCOPY, SCOPYM, TQLRAT,
!                    XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800808  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  CHIEV
  INTEGER I,INFO,J,JOB,K,L,LDA,LDV,M,MDIM,N
  REAL A(*),E(*),WORK(*),V(*)
!***FIRST EXECUTABLE STATEMENT  CHIEV
  if (N  >  LDA) call XERMSG ('SLATEC', 'CHIEV', 'N  >  LDA.', 1, &
     1)
  if ( N  >  LDA) RETURN
  if (N  <  1) call XERMSG ('SLATEC', 'CHIEV', 'N  <  1', 2, 1)
  if ( N  <  1) RETURN
  if ( N  ==  1 .AND. JOB  ==  0) go to 35
  MDIM = 2 * LDA
  if ( JOB  ==  0) go to 5
  if (N  >  LDV) call XERMSG ('SLATEC', 'CHIEV', &
     'JOB  /=  0 AND N  >  LDV.', 3, 1)
  if ( N  >  LDV) RETURN
  if ( N  ==  1) go to 35
!
!       REARRANGE A if NECESSARY WHEN LDA > LDV AND JOB  /= 0
!
  MDIM = MIN(MDIM,2 * LDV)
  if (LDA  <  LDV) call XERMSG ('SLATEC', 'CHIEV', &
     'LDA < LDV,  ELEMENTS OF V OTHER THAN THE N BY N OUTPUT ' // &
     'ELEMENTS HAVE BEEN CHANGED.', 5, 0)
  if ( LDA <= LDV) go to 5
  call XERMSG ('SLATEC', 'CHIEV', &
     'LDA > LDV, ELEMENTS OF A OTHER THAN THE N BY N INPUT ' // &
     'ELEMENTS HAVE BEEN CHANGED.', 4, 0)
  L = N - 1
  DO 4 J=1,L
     M = 1+J*2*LDV
     K = 1+J*2*LDA
     call SCOPY(2*N,A(K),1,A(M),1)
    4 CONTINUE
    5 CONTINUE
!
!     FILL IN LOWER TRIANGLE OF A, COLUMN BY COLUMN.
!
  DO 6 J = 1,N
   K = (J-1)*(MDIM+2)+1
   if (A(K+1)  /=  0.0) call XERMSG ('SLATEC', 'CHIEV', &
      'NONREAL ELEMENT ON DIAGONAL OF A', 6, 1)
  if ( A(K+1)  /= 0.0) RETURN
   call SCOPY(N-J+1,A(K),MDIM,A(K),2)
   call SCOPYM(N-J+1,A(K+1),MDIM,A(K+1),2)
    6 CONTINUE
!
!     SEPARATE REAL AND IMAGINARY PARTS
!
  DO 10 J = 1, N
   K = (J-1) * MDIM +1
   L = K + N
   call SCOPY(N,A(K+1),2,WORK(1),1)
   call SCOPY(N,A(K),2,A(K),1)
   call SCOPY(N,WORK(1),1,A(L),1)
   10 CONTINUE
!
!    REDUCE A TO TRIDIAGONAL MATRIX.
!
  call HTRIDI(MDIM,N,A(1),A(N+1),E,WORK(1),WORK(N+1), &
              WORK(2*N+1))
  if ( JOB  /=  0) GOTO 15
!
!     EIGENVALUES ONLY.
!
  call TQLRAT(N,E,WORK(N+1),INFO)
  return
!
!     EIGENVALUES AND EIGENVECTORS.
!
   15 DO 17 J = 1,N
   K = (J-1) * MDIM + 1
   M = K + N - 1
   DO 16 I = K,M
   16   V(I) = 0.
   I = K + J - 1
   V(I) = 1.
   17 CONTINUE
  call IMTQL2(MDIM,N,E,WORK(1),V,INFO)
  if ( INFO  /=  0) RETURN
  call HTRIBK(MDIM,N,A(1),A(N+1),WORK(2*N+1),N,V(1),V(N+1))
!
!    CONVERT EIGENVECTORS TO COMPLEX STORAGE.
!
  DO 20 J = 1,N
   K = (J-1) * MDIM + 1
   I = (J-1) * 2 * LDV + 1
   L = K + N
   call SCOPY(N,V(K),1,WORK(1),1)
   call SCOPY(N,V(L),1,V(I+1),2)
   call SCOPY(N,WORK(1),1,V(I),2)
   20 CONTINUE
  return
!
!     TAKE CARE OF N=1 CASE.
!
   35 if (A(2)  /=  0.) call XERMSG ('SLATEC', 'CHIEV', &
     'NONREAL ELEMENT ON DIAGONAL OF A', 6, 1)
  if ( A(2)  /=  0.) RETURN
  E(1) = A(1)
  INFO = 0
  if ( JOB  ==  0) RETURN
  V(1) = A(1)
  V(2) = 0.
  return
end
