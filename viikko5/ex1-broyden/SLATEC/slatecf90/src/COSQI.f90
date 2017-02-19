subroutine COSQI (N, WSAVE)
!
!! COSQI initializes a work array for COSQF and COSQB.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A3
!***TYPE      SINGLE PRECISION (COSQI-S)
!***KEYWORDS  COSINE FOURIER TRANSFORM, FFTPACK
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine COSQI initializes the work array WSAVE which is used in
!  both COSQF1 and COSQB1.  The prime factorization of N together with
!  a tabulation of the trigonometric functions are computed and
!  stored in WSAVE.
!
!  Input Parameter
!
!  N       the length of the array to be transformed.  The method
!          is most efficient when N is a product of small primes.
!
!  Output Parameter
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15.
!          The same work array can be used for both COSQF1 and COSQB1
!          as long as N remains unchanged.  Different WSAVE arrays
!          are required for different values of N.  The contents of
!          WSAVE must not be changed between calls of COSQF1 or COSQB1.
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  RFFTI
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           (a) changing dummy array size declarations (1) to (*),
!           (b) changing references to intrinsic function FLOAT
!               to REAL, and
!           (c) changing definition of variable PIH by using
!               FORTRAN intrinsic function ATAN instead of a DATA
!               statement.
!   881128  Modified by Dick Valent to meet prologue standards.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  COSQI
  DIMENSION WSAVE(*)
!***FIRST EXECUTABLE STATEMENT  COSQI
  PIH = 2.*ATAN(1.)
  DT = PIH/N
  FK = 0.
  DO 101 K=1,N
     FK = FK+1.
     WSAVE(K) = COS(FK*DT)
  101 CONTINUE
  call RFFTI (N,WSAVE(N+1))
  return
end
