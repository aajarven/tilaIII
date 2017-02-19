subroutine SINQB (N, X, WSAVE)
!
!! SINQB computes the unnormalized inverse of SINQF.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A3
!***TYPE      SINGLE PRECISION (SINQB-S)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine SINQB computes the fast Fourier transform of quarter
!  wave data.  That is, SINQB computes a sequence from its
!  representation in terms of a sine series with odd wave numbers.
!  the transform is defined below at output parameter X.
!
!  SINQF is the unnormalized inverse of SINQB since a call of SINQB
!  followed by a call of SINQF will multiply the input sequence X
!  by 4*N.
!
!  The array WSAVE which is used by subroutine SINQB must be
!  initialized by calling subroutine SINQI(N,WSAVE).
!
!  Input Parameters
!
!  N       the length of the array X to be transformed.  The method
!          is most efficient when N is a product of small primes.
!
!  X       an array which contains the sequence to be transformed
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15
!          in the program that calls SINQB.  The WSAVE array must be
!          initialized by calling subroutine SINQI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!
!  Output Parameters
!
!  X       For I=1,...,N
!
!               X(I)= the sum from K=1 to K=N of
!
!                 4*X(K)*SIN((2*K-1)*I*PI/(2*N))
!
!               a call of SINQB followed by a call of
!               SINQF will multiply the sequence X by 4*N.
!               Therefore SINQF is the unnormalized inverse
!               of SINQB.
!
!  WSAVE   contains initialization calculations which must not
!          be destroyed between calls of SINQB or SINQF.
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  COSQB
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*).
!   861211  REVISION DATE from Version 3.2
!   881128  Modified by Dick Valent to meet prologue standards.
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SINQB
  DIMENSION X(*), WSAVE(*)
!***FIRST EXECUTABLE STATEMENT  SINQB
  if (N  >  1) go to 101
  X(1) = 4.*X(1)
  return
  101 NS2 = N/2
  DO 102 K=2,N,2
     X(K) = -X(K)
  102 CONTINUE
  call COSQB (N,X,WSAVE)
  DO 103 K=1,NS2
     KC = N-K
     XHOLD = X(K)
     X(K) = X(KC+1)
     X(KC+1) = XHOLD
  103 CONTINUE
  return
end
