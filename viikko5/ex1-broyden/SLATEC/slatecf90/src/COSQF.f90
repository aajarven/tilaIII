subroutine COSQF (N, X, WSAVE)
!
!! COSQF computes the forward cosine transform with odd wave numbers.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A3
!***TYPE      SINGLE PRECISION (COSQF-S)
!***KEYWORDS  COSINE FOURIER TRANSFORM, FFTPACK
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine COSQF computes the fast Fourier transform of quarter
!  wave data. That is, COSQF computes the coefficients in a cosine
!  series representation with only odd wave numbers.  The transform
!  is defined below at Output Parameter X
!
!  COSQF is the unnormalized inverse of COSQB since a call of COSQF
!  followed by a call of COSQB will multiply the input sequence X
!  by 4*N.
!
!  The array WSAVE which is used by subroutine COSQF must be
!  initialized by calling subroutine COSQI(N,WSAVE).
!
!
!  Input Parameters
!
!  N       the length of the array X to be transformed.  The method
!          is most efficient when N is a product of small primes.
!
!  X       an array which contains the sequence to be transformed
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15
!          in the program that calls COSQF.  The WSAVE array must be
!          initialized by calling subroutine COSQI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!
!  Output Parameters
!
!  X       For I=1,...,N
!
!               X(I) = X(1) plus the sum from K=2 to K=N of
!
!                  2*X(K)*COS((2*I-1)*(K-1)*PI/(2*N))
!
!               A call of COSQF followed by a call of
!               COSQB will multiply the sequence X by 4*N.
!               Therefore COSQB is the unnormalized inverse
!               of COSQF.
!
!  WSAVE   contains initialization calculations which must not
!          be destroyed between calls of COSQF or COSQB.
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  COSQF1
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           (a) changing dummy array size declarations (1) to (*),
!           (b) changing definition of variable SQRT2 by using
!               FORTRAN intrinsic function SQRT instead of a DATA
!               statement.
!   861211  REVISION DATE from Version 3.2
!   881128  Modified by Dick Valent to meet prologue standards.
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  COSQF
  DIMENSION X(*), WSAVE(*)
!***FIRST EXECUTABLE STATEMENT  COSQF
  SQRT2 = SQRT(2.)
  if (N-2) 102,101,103
  101 TSQX = SQRT2*X(2)
  X(2) = X(1)-TSQX
  X(1) = X(1)+TSQX
  102 RETURN
  103 call COSQF1 (N,X,WSAVE,WSAVE(N+1))
  return
end
