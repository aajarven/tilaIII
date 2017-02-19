subroutine CPEVL (N, M, A, Z, C, B, KBD)
!
!! CPEVL is subsidiary to CPZERO.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (CPEVL-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!        Evaluate a complex polynomial and its derivatives.
!        Optionally compute error bounds for these values.
!
!   INPUT...
!        N = Degree of the polynomial
!        M = Number of derivatives to be calculated,
!            M=0 evaluates only the function
!            M=1 evaluates the function and first derivative, etc.
!             if M  >  N+1 function and all N derivatives will be
!                calculated.
!       A = Complex vector containing the N+1 coefficients of polynomial
!               A(I)= coefficient of Z**(N+1-I)
!        Z = Complex point at which the evaluation is to take place.
!        C = Array of 2(M+1) words into which values are placed.
!        B = Array of 2(M+1) words only needed if bounds are to be
!              calculated.  It is not used otherwise.
!        KBD = A logical variable, e.g. .TRUE. or .FALSE. which is
!              to be set .TRUE. if bounds are to be computed.
!
!  OUTPUT...
!        C =  C(I+1) contains the complex value of the I-th
!              derivative at Z, I=0,...,M
!        B =  B(I) contains the bounds on the real and imaginary parts
!              of C(I) if they were requested.
!
!***SEE ALSO  CPZERO
!***ROUTINES CALLED  I1MACH
!***REVISION HISTORY  (YYMMDD)
!   810223  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  CPEVL
!
  COMPLEX A(*),C(*),Z,CI,CIM1,B(*),BI,BIM1,T,ZA,Q
  LOGICAL KBD
  SAVE D1
  DATA D1 /0.0/
  ZA(Q)=CMPLX(ABS(REAL(Q)),ABS(AIMAG(Q)))
!***FIRST EXECUTABLE STATEMENT  CPEVL
  if (D1  ==  0.0) D1 = REAL(I1MACH(10))**(1-I1MACH(11))
  NP1=N+1
  DO 1 J=1,NP1
     CI=0.0
     CIM1=A(J)
     BI=0.0
     BIM1=0.0
     MINI=MIN(M+1,N+2-J)
        DO 1 I=1,MINI
           if ( J  /=  1) CI=C(I)
           if ( I  /=  1) CIM1=C(I-1)
           C(I)=CIM1+Z*CI
           if ( .NOT. KBD) go to 1
           if ( J  /=  1) BI=B(I)
           if ( I  /=  1) BIM1=B(I-1)
           T=BI+(3.*D1+4.*D1*D1)*ZA(CI)
           R=REAL(ZA(Z)*CMPLX(REAL(T),-AIMAG(T)))
           S=AIMAG(ZA(Z)*T)
           B(I)=(1.+8.*D1)*(BIM1+D1*ZA(CIM1)+CMPLX(R,S))
           if ( J  ==  1) B(I)=0.0
    1 CONTINUE
  return
end
