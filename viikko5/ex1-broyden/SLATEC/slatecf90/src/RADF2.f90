subroutine RADF2 (IDO, L1, CC, CH, WA1)
!
!! RADF2 calculates the fast Fourier transform of subvectors of length two.
!
!***LIBRARY   SLATEC (FFTPACK)
!***TYPE      SINGLE PRECISION (RADF2-S)
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
!***END PROLOGUE  RADF2
  DIMENSION CH(IDO,2,*), CC(IDO,L1,2), WA1(*)
!***FIRST EXECUTABLE STATEMENT  RADF2
  DO 101 K=1,L1
     CH(1,1,K) = CC(1,K,1)+CC(1,K,2)
     CH(IDO,2,K) = CC(1,K,1)-CC(1,K,2)
  101 CONTINUE
  if (IDO-2) 107,105,102
  102 IDP2 = IDO+2
  if ( (IDO-1)/2 < L1) go to 108
  DO 104 K=1,L1
!DIR$ IVDEP
     DO 103 I=3,IDO,2
        IC = IDP2-I
        TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        CH(I,1,K) = CC(I,K,1)+TI2
        CH(IC,2,K) = TI2-CC(I,K,1)
        CH(I-1,1,K) = CC(I-1,K,1)+TR2
        CH(IC-1,2,K) = CC(I-1,K,1)-TR2
  103    CONTINUE
  104 CONTINUE
  go to 111
  108 DO 110 I=3,IDO,2
     IC = IDP2-I
!DIR$ IVDEP
     DO 109 K=1,L1
        TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        CH(I,1,K) = CC(I,K,1)+TI2
        CH(IC,2,K) = TI2-CC(I,K,1)
        CH(I-1,1,K) = CC(I-1,K,1)+TR2
        CH(IC-1,2,K) = CC(I-1,K,1)-TR2
  109    CONTINUE
  110 CONTINUE
  111 if (MOD(IDO,2)  ==  1) RETURN
  105 DO 106 K=1,L1
     CH(1,2,K) = -CC(IDO,K,2)
     CH(IDO,1,K) = CC(IDO,K,1)
  106 CONTINUE
  107 RETURN
end
