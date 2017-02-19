subroutine CPZERO (IN, A, R, T, IFLG, S)
!
!! CPZERO finds the zeros of a polynomial with complex coefficients.
!
!***LIBRARY   SLATEC
!***CATEGORY  F1A1B
!***TYPE      COMPLEX (RPZERO-S, CPZERO-C)
!***KEYWORDS  POLYNOMIAL ROOTS, POLYNOMIAL ZEROS, REAL ROOTS
!***AUTHOR  Kahaner, D. K., (NBS)
!***DESCRIPTION
!
!      Find the zeros of the complex polynomial
!         P(Z)= A(1)*Z**N + A(2)*Z**(N-1) +...+ A(N+1)
!
!    Input...
!       IN = degree of P(Z)
!       A = complex vector containing coefficients of P(Z),
!            A(I) = coefficient of Z**(N+1-i)
!       R = N word complex vector containing initial estimates for zeros
!            if these are known.
!       T = 4(N+1) word array used for temporary storage
!       IFLG = flag to indicate if initial estimates of
!              zeros are input.
!            If IFLG  ==  0, no estimates are input.
!            If IFLG  /=  0, the vector R contains estimates of
!               the zeros
!       ** WARNING ****** If estimates are input, they must
!                         be separated, that is, distinct or
!                         not repeated.
!       S = an N word array
!
!    Output...
!       R(I) = Ith zero,
!       S(I) = bound for R(I) .
!       IFLG = error diagnostic
!    Error Diagnostics...
!       If IFLG  ==  0 on return, all is well
!       If IFLG  ==  1 on return, A(1)=0.0 or N=0 on input
!       If IFLG  ==  2 on return, the program failed to converge
!                after 25*N iterations.  Best current estimates of the
!                zeros are in R(I).  Error bounds are not calculated.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CPEVL
!***REVISION HISTORY  (YYMMDD)
!   810223  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CPZERO
!
  REAL  S(*)
  COMPLEX R(*),T(*),A(*),PN,TEMP
!***FIRST EXECUTABLE STATEMENT  CPZERO
  if (  IN  <=  0 .OR. ABS(A(1))  ==  0.0 ) go to 30
!
!       CHECK FOR EASILY OBTAINED ZEROS
!
  N=IN
  N1=N+1
  if ( IFLG  /=  0) go to 14
    1 N1=N+1
  if ( N  >  1) go to 2
     R(1)=-A(2)/A(1)
     S(1)=0.0
     return
    2 if (  ABS(A(N1))  /=  0.0 ) go to 3
     R(N)=0.0
     S(N)=0.0
     N=N-1
     go to 1
!
!          if INITIAL ESTIMATES FOR ZEROS NOT GIVEN, FIND SOME
!
    3 TEMP=-A(2)/(A(1)*N)
  call CPEVL(N,N,A,TEMP,T,T,.FALSE.)
  IMAX=N+2
  T(N1)=ABS(T(N1))
  DO 6 I=2,N1
     T(N+I)=-ABS(T(N+2-I))
     if ( REAL(T(N+I))  <  REAL(T(IMAX))) IMAX=N+I
    6 CONTINUE
  X=(-REAL(T(IMAX))/REAL(T(N1)))**(1./(IMAX-N1))
    7 X=2.*X
     call CPEVL(N,0,T(N1),CMPLX(X,0.0),PN,PN,.FALSE.)
  if (REAL(PN) < 0.) go to 7
  U=.5*X
  V=X
   10 X=.5*(U+V)
     call CPEVL(N,0,T(N1),CMPLX(X,0.0),PN,PN,.FALSE.)
     if (REAL(PN) > 0.) V=X
     if (REAL(PN) <= 0.) U=X
     if ( (V-U)  >  .001*(1.+V)) go to 10
  DO 13 I=1,N
     U=(3.14159265/N)*(2*I-1.5)
   13    R(I)=MAX(X,.001*ABS(TEMP))*CMPLX(COS(U),SIN(U))+TEMP
!
!          MAIN ITERATION LOOP STARTS HERE
!
   14 NR=0
  NMAX=25*N
  DO 19 NIT=1,NMAX
     DO 18 I=1,N
        if ( NIT  /=  1 .AND. ABS(T(I))  ==  0.) go to 18
           call CPEVL(N,0,A,R(I),PN,TEMP,.TRUE.)
           if ( ABS(REAL(PN))+ABS(AIMAG(PN))  >  REAL(TEMP)+ &
                AIMAG(TEMP)) go to 16
              T(I)=0.0
              NR=NR+1
              go to 18
   16          TEMP=A(1)
           DO 17 J=1,N
   17             if ( J  /=  I) TEMP=TEMP*(R(I)-R(J))
           T(I)=PN/TEMP
   18    CONTINUE
     DO 15 I=1,N
   15       R(I)=R(I)-T(I)
     if ( NR  ==  N) go to 21
   19 CONTINUE
  go to 26
!
!          CALCULATE ERROR BOUNDS FOR ZEROS
!
   21 DO 25 NR=1,N
     call CPEVL(N,N,A,R(NR),T,T(N+2),.TRUE.)
     X=ABS(CMPLX(ABS(REAL(T(1))),ABS(AIMAG(T(1))))+T(N+2))
     S(NR)=0.0
     DO 23 I=1,N
        X=X*REAL(N1-I)/I
        TEMP=CMPLX(MAX(ABS(REAL(T(I+1)))-REAL(T(N1+I)),0.0), &
             MAX(ABS(AIMAG(T(I+1)))-AIMAG(T(N1+I)),0.0))
   23       S(NR)=MAX(S(NR),(ABS(TEMP)/X)**(1./I))
   25    S(NR)=1./S(NR)
  return
!        ERROR EXITS
   26 IFLG=2
  return
   30 IFLG=1
  return
end
