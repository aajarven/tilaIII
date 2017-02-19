subroutine RADF3 (IDO, L1, CC, CH, WA1, WA2)
!
!! RADF3 calculates the fast Fourier transform of subvectors of length three.
!
!***LIBRARY   SLATEC (FFTPACK)
!***TYPE      SINGLE PRECISION (RADF3-S)
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           (a) changing dummy array size declarations (1) to (*),
!           (b) changing definition of variable TAUI by using
!               FORTRAN intrinsic function SQRT instead of a DATA
!               statement.
!   881128  Modified by Dick Valent to meet prologue standards.
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  RADF3
  DIMENSION CH(IDO,3,*), CC(IDO,L1,3), WA1(*), WA2(*)
!***FIRST EXECUTABLE STATEMENT  RADF3
  TAUR = -.5
  TAUI = .5*SQRT(3.)
  DO 101 K=1,L1
     CR2 = CC(1,K,2)+CC(1,K,3)
     CH(1,1,K) = CC(1,K,1)+CR2
     CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
     CH(IDO,2,K) = CC(1,K,1)+TAUR*CR2
  101 CONTINUE
  if (IDO  ==  1) RETURN
  IDP2 = IDO+2
  if ( (IDO-1)/2 < L1) go to 104
  DO 103 K=1,L1
!DIR$ IVDEP
     DO 102 I=3,IDO,2
        IC = IDP2-I
        DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
        DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
        CR2 = DR2+DR3
        CI2 = DI2+DI3
        CH(I-1,1,K) = CC(I-1,K,1)+CR2
        CH(I,1,K) = CC(I,K,1)+CI2
        TR2 = CC(I-1,K,1)+TAUR*CR2
        TI2 = CC(I,K,1)+TAUR*CI2
        TR3 = TAUI*(DI2-DI3)
        TI3 = TAUI*(DR3-DR2)
        CH(I-1,3,K) = TR2+TR3
        CH(IC-1,2,K) = TR2-TR3
        CH(I,3,K) = TI2+TI3
        CH(IC,2,K) = TI3-TI2
  102    CONTINUE
  103 CONTINUE
  return
  104 DO 106 I=3,IDO,2
     IC = IDP2-I
!DIR$ IVDEP
     DO 105 K=1,L1
        DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
        DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
        CR2 = DR2+DR3
        CI2 = DI2+DI3
        CH(I-1,1,K) = CC(I-1,K,1)+CR2
        CH(I,1,K) = CC(I,K,1)+CI2
        TR2 = CC(I-1,K,1)+TAUR*CR2
        TI2 = CC(I,K,1)+TAUR*CI2
        TR3 = TAUI*(DI2-DI3)
        TI3 = TAUI*(DR3-DR2)
        CH(I-1,3,K) = TR2+TR3
        CH(IC-1,2,K) = TR2-TR3
        CH(I,3,K) = TI2+TI3
        CH(IC,2,K) = TI3-TI2
  105    CONTINUE
  106 CONTINUE
  return
end
