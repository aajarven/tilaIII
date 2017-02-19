subroutine SINT (N, X, WSAVE)
!
!! SINT computes the sine transform of a real, odd sequence.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A3
!***TYPE      SINGLE PRECISION (SINT-S)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine SINT computes the discrete Fourier sine transform
!  of an odd sequence X(I).  The transform is defined below at
!  output parameter X.
!
!  SINT is the unnormalized inverse of itself since a call of SINT
!  followed by another call of SINT will multiply the input sequence
!  X by 2*(N+1).
!
!  The array WSAVE which is used by subroutine SINT must be
!  initialized by calling subroutine SINTI(N,WSAVE).
!
!  Input Parameters
!
!  N       the length of the sequence to be transformed.  The method
!          is most efficient when N+1 is the product of small primes.
!
!  X       an array which contains the sequence to be transformed
!
!
!  WSAVE   a work array with dimension at least INT(3.5*N+16)
!          in the program that calls SINT.  The WSAVE array must be
!          initialized by calling subroutine SINTI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!
!  Output Parameters
!
!  X       For I=1,...,N
!
!               X(I)= the sum from K=1 to K=N
!
!                    2*X(K)*SIN(K*I*PI/(N+1))
!
!               A call of SINT followed by another call of
!               SINT will multiply the sequence X by 2*(N+1).
!               Hence SINT is the unnormalized inverse
!               of itself.
!
!  WSAVE   contains initialization calculations which must not be
!          destroyed between calls of SINT.
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
!           (b) changing definition of variable SQRT3 by using
!               FORTRAN intrinsic function SQRT instead of a DATA
!               statement.
!   881128  Modified by Dick Valent to meet prologue standards.
!   891009  Removed unreferenced statement label.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SINT
  DIMENSION X(*), WSAVE(*)
!***FIRST EXECUTABLE STATEMENT  SINT
  SQRT3 = SQRT(3.)
  if (N-2) 101,102,103
  101 X(1) = X(1)+X(1)
  return
  102 XH = SQRT3*(X(1)+X(2))
  X(2) = SQRT3*(X(1)-X(2))
  X(1) = XH
  return
  103 NP1 = N+1
  NS2 = N/2
  WSAVE(1) = 0.
  KW = NP1
  DO 104 K=1,NS2
     KW = KW+1
     KC = NP1-K
     T1 = X(K)-X(KC)
     T2 = WSAVE(KW)*(X(K)+X(KC))
     WSAVE(K+1) = T1+T2
     WSAVE(KC+1) = T2-T1
  104 CONTINUE
  MODN = MOD(N,2)
  if (MODN  /=  0) WSAVE(NS2+2) = 4.*X(NS2+1)
  NF = NP1+NS2+1
  call RFFTF (NP1,WSAVE,WSAVE(NF))
  X(1) = .5*WSAVE(1)
  DO 105 I=3,N,2
     X(I-1) = -WSAVE(I)
     X(I) = X(I-2)+WSAVE(I-1)
  105 CONTINUE
  if (MODN  /=  0) RETURN
  X(N) = -WSAVE(N+1)
  return
end
