subroutine CFFTF1 (N, C, CH, WA, IFAC)
!
!! CFFTF1 computes the forward transform of a complex, periodic sequence.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A2
!***TYPE      COMPLEX (RFFTF1-S, CFFTF1-C)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine CFFTF1 computes the forward complex discrete Fourier
!  transform (the Fourier analysis).  Equivalently, CFFTF1 computes
!  the Fourier coefficients of a complex periodic sequence.
!  The transform is defined below at output parameter C.
!
!  The transform is not normalized.  To obtain a normalized transform
!  the output must be divided by N.  Otherwise a call of CFFTF1
!  followed by a call of CFFTB1 will multiply the sequence by N.
!
!  The arrays WA and IFAC which are used by subroutine CFFTB1 must be
!  initialized by calling subroutine CFFTI1 (N, WA, IFAC).
!
!  Input Parameters
!
!  N       the length of the complex sequence C.  The method is
!          more efficient when N is the product of small primes.
!
!  C       a complex array of length N which contains the sequence
!
!  CH      a real work array of length at least 2*N
!
!  WA      a real work array which must be dimensioned at least 2*N.
!
!  IFAC    an integer work array which must be dimensioned at least 15.
!
!          The WA and IFAC arrays must be initialized by calling
!          subroutine CFFTI1 (N, WA, IFAC), and different WA and IFAC
!          arrays must be used for each different value of N.  This
!          initialization does not have to be repeated so long as N
!          remains unchanged.  Thus subsequent transforms can be
!          obtained faster than the first.  The same WA and IFAC arrays
!          can be used by CFFTF1 and CFFTB1.
!
!  Output Parameters
!
!  C       For J=1,...,N
!
!              C(J)=the sum from K=1,...,N of
!
!                 C(K)*EXP(-I*(J-1)*(K-1)*2*PI/N)
!
!                         where I=SQRT(-1)
!
!  NOTE:   WA and IFAC contain initialization calculations which must
!          not be destroyed between calls of subroutine CFFTF1 or CFFTB1
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  PASSF, PASSF2, PASSF3, PASSF4, PASSF5
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*).
!   881128  Modified by Dick Valent to meet prologue standards.
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900131  Routine changed from subsidiary to user-callable.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CFFTF1
  DIMENSION CH(*), C(*), WA(*), IFAC(*)
!***FIRST EXECUTABLE STATEMENT  CFFTF1
  NF = IFAC(2)
  NA = 0
  L1 = 1
  IW = 1
  DO 116 K1=1,NF
     IP = IFAC(K1+2)
     L2 = IP*L1
     IDO = N/L2
     IDOT = IDO+IDO
     IDL1 = IDOT*L1
     if (IP  /=  4) go to 103
     IX2 = IW+IDOT
     IX3 = IX2+IDOT
     if (NA  /=  0) go to 101
     call PASSF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
     go to 102
  101    call PASSF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
     go to 115
  103    if (IP  /=  2) go to 106
     if (NA  /=  0) go to 104
     call PASSF2 (IDOT,L1,C,CH,WA(IW))
     go to 105
  104    call PASSF2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
     go to 115
  106    if (IP  /=  3) go to 109
     IX2 = IW+IDOT
     if (NA  /=  0) go to 107
     call PASSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
     go to 108
  107    call PASSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
     go to 115
  109    if (IP  /=  5) go to 112
     IX2 = IW+IDOT
     IX3 = IX2+IDOT
     IX4 = IX3+IDOT
     if (NA  /=  0) go to 110
     call PASSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
     go to 111
  110    call PASSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
     go to 115
  112    if (NA  /=  0) go to 113
     call PASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
     go to 114
  113    call PASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    if (NAC  /=  0) NA = 1-NA
  115    L1 = L2
     IW = IW+(IP-1)*IDOT
  116 CONTINUE
  if (NA  ==  0) RETURN
  N2 = N+N
  DO 117 I=1,N2
     C(I) = CH(I)
  117 CONTINUE
  return
end
