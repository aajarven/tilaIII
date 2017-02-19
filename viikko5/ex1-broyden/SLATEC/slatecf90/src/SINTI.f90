subroutine SINTI (N, WSAVE)
!
!! SINTI initializes a work array for SINT.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A3
!***TYPE      SINGLE PRECISION (SINTI-S)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine SINTI initializes the array WSAVE which is used in
!  subroutine SINT.  The prime factorization of N together with
!  a tabulation of the trigonometric functions are computed and
!  stored in WSAVE.
!
!  Input Parameter
!
!  N       the length of the sequence to be transformed.  The method
!          is most efficient when N+1 is a product of small primes.
!
!  Output Parameter
!
!  WSAVE   a work array with at least INT(3.5*N+16) locations.
!          Different WSAVE arrays are required for different values
!          of N.  The contents of WSAVE must not be changed between
!          calls of SINT.
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
!           (c) changing definition of variable PI by using
!               FORTRAN intrinsic function ATAN instead of a DATA
!               statement.
!   881128  Modified by Dick Valent to meet prologue standards.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SINTI
  DIMENSION WSAVE(*)
!***FIRST EXECUTABLE STATEMENT  SINTI
  if (N  <=  1) RETURN
  PI = 4.*ATAN(1.)
  NP1 = N+1
  NS2 = N/2
  DT = PI/NP1
  KS = N+2
  KF = KS+NS2-1
  FK = 0.
  DO 101 K=KS,KF
     FK = FK+1.
     WSAVE(K) = 2.*SIN(FK*DT)
  101 CONTINUE
  call RFFTI (NP1,WSAVE(KF+1))
  return
end
