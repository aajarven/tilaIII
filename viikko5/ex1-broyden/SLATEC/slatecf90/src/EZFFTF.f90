subroutine EZFFTF (N, R, AZERO, A, B, WSAVE)
!
!! EZFFTF computes a simplified real, periodic, fast Fourier forward transform.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A1
!***TYPE      SINGLE PRECISION (EZFFTF-S)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine EZFFTF computes the Fourier coefficients of a real
!  periodic sequence (Fourier analysis).  The transform is defined
!  below at Output Parameters AZERO, A and B.  EZFFTF is a simplified
!  but slower version of RFFTF.
!
!  Input Parameters
!
!  N       the length of the array R to be transformed.  The method
!          is most efficient when N is the product of small primes.
!
!  R       a real array of length N which contains the sequence
!          to be transformed.  R is not destroyed.
!
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15
!          in the program that calls EZFFTF.  The WSAVE array must be
!          initialized by calling subroutine EZFFTI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!          The same WSAVE array can be used by EZFFTF and EZFFTB.
!
!  Output Parameters
!
!  AZERO   the sum from I=1 to I=N of R(I)/N
!
!  A,B     for N even B(N/2)=0. and A(N/2) is the sum from I=1 to
!          I=N of (-1)**(I-1)*R(I)/N
!
!          for N even define KMAX=N/2-1
!          for N odd  define KMAX=(N-1)/2
!
!          then for  K=1,...,KMAX
!
!               A(K) equals the sum from I=1 to I=N of
!
!                    2./N*R(I)*COS(K*(I-1)*2*PI/N)
!
!               B(K) equals the sum from I=1 to I=N of
!
!                    2./N*R(I)*SIN(K*(I-1)*2*PI/N)
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  RFFTF
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           (a) changing dummy array size declarations (1) to (*),
!           (b) changing references to intrinsic function FLOAT
!               to REAL.
!   881128  Modified by Dick Valent to meet prologue standards.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  EZFFTF
  DIMENSION R(*), A(*), B(*), WSAVE(*)
!***FIRST EXECUTABLE STATEMENT  EZFFTF
  if (N-2) 101,102,103
  101 AZERO = R(1)
  return
  102 AZERO = .5*(R(1)+R(2))
  A(1) = .5*(R(1)-R(2))
  return
  103 DO 104 I=1,N
     WSAVE(I) = R(I)
  104 CONTINUE
  call RFFTF (N,WSAVE,WSAVE(N+1))
  CF = 2./N
  CFM = -CF
  AZERO = .5*CF*WSAVE(1)
  NS2 = (N+1)/2
  NS2M = NS2-1
  DO 105 I=1,NS2M
     A(I) = CF*WSAVE(2*I)
     B(I) = CFM*WSAVE(2*I+1)
  105 CONTINUE
  if (MOD(N,2)  ==  0) A(NS2) = .5*CF*WSAVE(N)
  return
end
