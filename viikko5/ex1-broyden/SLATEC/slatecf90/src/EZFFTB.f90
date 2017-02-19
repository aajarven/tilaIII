subroutine EZFFTB (N, R, AZERO, A, B, WSAVE)
!
!! EZFFTB is a simplified real, periodic, backward fast Fourier transform.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A1
!***TYPE      SINGLE PRECISION (EZFFTB-S)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine EZFFTB computes a real periodic sequence from its
!  Fourier coefficients (Fourier synthesis).  The transform is
!  defined below at Output Parameter R.  EZFFTB is a simplified
!  but slower version of RFFTB.
!
!  Input Parameters
!
!  N       the length of the output array R.  The method is most
!          efficient when N is the product of small primes.
!
!  AZERO   the constant Fourier coefficient
!
!  A,B     arrays which contain the remaining Fourier coefficients.
!          These arrays are not destroyed.
!
!          The length of these arrays depends on whether N is even or
!          odd.
!
!          If N is even, N/2    locations are required.
!          If N is odd, (N-1)/2 locations are required
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15
!          in the program that calls EZFFTB.  The WSAVE array must be
!          initialized by calling subroutine EZFFTI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!          The same WSAVE array can be used by EZFFTF and EZFFTB.
!
!  Output Parameters
!
!  R       if N is even, define KMAX=N/2
!          if N is odd,  define KMAX=(N-1)/2
!
!          Then for I=1,...,N
!
!               R(I)=AZERO plus the sum from K=1 to K=KMAX of
!
!               A(K)*COS(K*(I-1)*2*PI/N)+B(K)*SIN(K*(I-1)*2*PI/N)
!
!  ********************* Complex Notation **************************
!
!          For J=1,...,N
!
!          R(J) equals the sum from K=-KMAX to K=KMAX of
!
!               C(K)*EXP(I*K*(J-1)*2*PI/N)
!
!          where
!
!               C(K) = .5*CMPLX(A(K),-B(K))   for K=1,...,KMAX
!
!               C(-K) = CONJG(C(K))
!
!               C(0) = AZERO
!
!                    and I=SQRT(-1)
!
!  *************** Amplitude - Phase Notation ***********************
!
!          For I=1,...,N
!
!          R(I) equals AZERO plus the sum from K=1 to K=KMAX of
!
!               ALPHA(K)*COS(K*(I-1)*2*PI/N+BETA(K))
!
!          where
!
!               ALPHA(K) = SQRT(A(K)*A(K)+B(K)*B(K))
!
!               COS(BETA(K))=A(K)/ALPHA(K)
!
!               SIN(BETA(K))=-B(K)/ALPHA(K)
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  RFFTB
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*)
!   861211  REVISION DATE from Version 3.2
!   881128  Modified by Dick Valent to meet prologue standards.
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  EZFFTB
  DIMENSION R(*), A(*), B(*), WSAVE(*)
!***FIRST EXECUTABLE STATEMENT  EZFFTB
  if (N-2) 101,102,103
  101 R(1) = AZERO
  return
  102 R(1) = AZERO+A(1)
  R(2) = AZERO-A(1)
  return
  103 NS2 = (N-1)/2
  DO 104 I=1,NS2
     R(2*I) = .5*A(I)
     R(2*I+1) = -.5*B(I)
  104 CONTINUE
  R(1) = AZERO
  if (MOD(N,2)  ==  0) R(N) = A(NS2+1)
  call RFFTB (N,R,WSAVE(N+1))
  return
end
