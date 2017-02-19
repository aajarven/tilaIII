subroutine PASSB2 (IDO, L1, CC, CH, WA1)
!
!! PASSB2 calculates the fast Fourier transform of subvectors of length two.
!
!***LIBRARY   SLATEC (FFTPACK)
!***TYPE      SINGLE PRECISION (PASSB2-S)
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*).
!   881128  Modified by Dick Valent to meet prologue standards.
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PASSB2
  DIMENSION CC(IDO,2,*), CH(IDO,L1,2), WA1(*)
!***FIRST EXECUTABLE STATEMENT  PASSB2
  if (IDO  >  2) go to 102
  DO 101 K=1,L1
     CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
     CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
     CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
     CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
  return
  102 if ( IDO/2 < L1) go to 105
  DO 104 K=1,L1
!DIR$ IVDEP
     DO 103 I=2,IDO,2
        CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
        TR2 = CC(I-1,1,K)-CC(I-1,2,K)
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
        TI2 = CC(I,1,K)-CC(I,2,K)
        CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
        CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
  return
  105 DO 107 I=2,IDO,2
!DIR$ IVDEP
     DO 106 K=1,L1
        CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
        TR2 = CC(I-1,1,K)-CC(I-1,2,K)
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
        TI2 = CC(I,1,K)-CC(I,2,K)
        CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
        CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  106    CONTINUE
  107 CONTINUE
  return
end
