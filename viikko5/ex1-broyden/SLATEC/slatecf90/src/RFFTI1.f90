subroutine RFFTI1 (N, WA, IFAC)
!
!! RFFTI1 initializes a real and an integer work array for RFFTF1 and RFFTB1.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A1
!***TYPE      SINGLE PRECISION (RFFTI1-S, CFFTI1-C)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!   Subroutine RFFTI1 initializes the work arrays WA and IFAC which are
!   used in both RFFTF1 and RFFTB1.  The prime factorization of N and a
!   tabulation of the trigonometric functions are computed and stored in
!   IFAC and WA, respectively.
!
!   Input Argument
!
!   N       the length of the sequence to be transformed.
!
!   Output Arguments
!
!   WA      a real work array which must be dimensioned at least N.
!
!   IFAC    an integer work array which must be dimensioned at least 15.
!
!   The same work arrays can be used for both RFFTF1 and RFFTB1 as long
!   as N remains unchanged.  Different WA and IFAC arrays are required
!   for different values of N.  The contents of WA and IFAC must not be
!   changed between calls of RFFTF1 or RFFTB1.
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           (a) changing dummy array size declarations (1) to (*),
!           (b) changing references to intrinsic function FLOAT
!               to REAL, and
!           (c) changing definition of variable TPI by using
!               FORTRAN intrinsic functions instead of DATA
!               statements.
!   881128  Modified by Dick Valent to meet prologue standards.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900131  Routine changed from subsidiary to user-callable.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  RFFTI1
  DIMENSION WA(*), IFAC(*), NTRYH(4)
  SAVE NTRYH
  DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
!***FIRST EXECUTABLE STATEMENT  RFFTI1
  NL = N
  NF = 0
  J = 0
  101 J = J+1
  if (J-4) 102,102,103
  102 NTRY = NTRYH(J)
  go to 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
  NR = NL-NTRY*NQ
  if (NR) 101,105,101
  105 NF = NF+1
  IFAC(NF+2) = NTRY
  NL = NQ
  if (NTRY  /=  2) go to 107
  if (NF  ==  1) go to 107
  DO 106 I=2,NF
     IB = NF-I+2
     IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
  IFAC(3) = 2
  107 if (NL  /=  1) go to 104
  IFAC(1) = N
  IFAC(2) = NF
  TPI = 8.*ATAN(1.)
  ARGH = TPI/N
  IS = 0
  NFM1 = NF-1
  L1 = 1
  if (NFM1  ==  0) RETURN
  DO 110 K1=1,NFM1
     IP = IFAC(K1+2)
     LD = 0
     L2 = L1*IP
     IDO = N/L2
     IPM = IP-1
     DO 109 J=1,IPM
        LD = LD+L1
        I = IS
        ARGLD = LD*ARGH
        FI = 0.
        DO 108 II=3,IDO,2
           I = I+2
           FI = FI+1.
           ARG = FI*ARGLD
           WA(I-1) = COS(ARG)
           WA(I) = SIN(ARG)
  108       CONTINUE
        IS = IS+IDO
  109    CONTINUE
     L1 = L2
  110 CONTINUE
  return
end
