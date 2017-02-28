subroutine RADF5 (IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
!
!! RADF5 calculates the fast Fourier transform of subvectors of length five.
!
!***LIBRARY   SLATEC (FFTPACK)
!***TYPE      SINGLE PRECISION (RADF5-S)
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           (a) changing dummy array size declarations (1) to (*),
!           (b) changing definition of variables PI, TI11, TI12,
!               TR11, TR12 by using FORTRAN intrinsic functions ATAN
!               and SIN instead of DATA statements.
!   881128  Modified by Dick Valent to meet prologue standards.
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  RADF5
  DIMENSION CC(IDO,L1,5), CH(IDO,5,*), WA1(*), WA2(*), WA3(*), &
            WA4(*)
!***FIRST EXECUTABLE STATEMENT  RADF5
  PI = 4.*ATAN(1.)
  TR11 = SIN(.1*PI)
  TI11 = SIN(.4*PI)
  TR12 = -SIN(.3*PI)
  TI12 = SIN(.2*PI)
  DO 101 K=1,L1
     CR2 = CC(1,K,5)+CC(1,K,2)
     CI5 = CC(1,K,5)-CC(1,K,2)
     CR3 = CC(1,K,4)+CC(1,K,3)
     CI4 = CC(1,K,4)-CC(1,K,3)
     CH(1,1,K) = CC(1,K,1)+CR2+CR3
     CH(IDO,2,K) = CC(1,K,1)+TR11*CR2+TR12*CR3
     CH(1,3,K) = TI11*CI5+TI12*CI4
     CH(IDO,4,K) = CC(1,K,1)+TR12*CR2+TR11*CR3
     CH(1,5,K) = TI12*CI5-TI11*CI4
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
        DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
        DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
        DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
        DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
        CR2 = DR2+DR5
        CI5 = DR5-DR2
        CR5 = DI2-DI5
        CI2 = DI2+DI5
        CR3 = DR3+DR4
        CI4 = DR4-DR3
        CR4 = DI3-DI4
        CI3 = DI3+DI4
        CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
        CH(I,1,K) = CC(I,K,1)+CI2+CI3
        TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
        TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
        TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
        TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
        TR5 = TI11*CR5+TI12*CR4
        TI5 = TI11*CI5+TI12*CI4
        TR4 = TI12*CR5-TI11*CR4
        TI4 = TI12*CI5-TI11*CI4
        CH(I-1,3,K) = TR2+TR5
        CH(IC-1,2,K) = TR2-TR5
        CH(I,3,K) = TI2+TI5
        CH(IC,2,K) = TI5-TI2
        CH(I-1,5,K) = TR3+TR4
        CH(IC-1,4,K) = TR3-TR4
        CH(I,5,K) = TI3+TI4
        CH(IC,4,K) = TI4-TI3
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
        DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
        DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
        DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
        DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
        CR2 = DR2+DR5
        CI5 = DR5-DR2
        CR5 = DI2-DI5
        CI2 = DI2+DI5
        CR3 = DR3+DR4
        CI4 = DR4-DR3
        CR4 = DI3-DI4
        CI3 = DI3+DI4
        CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
        CH(I,1,K) = CC(I,K,1)+CI2+CI3
        TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
        TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
        TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
        TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
        TR5 = TI11*CR5+TI12*CR4
        TI5 = TI11*CI5+TI12*CI4
        TR4 = TI12*CR5-TI11*CR4
        TI4 = TI12*CI5-TI11*CI4
        CH(I-1,3,K) = TR2+TR5
        CH(IC-1,2,K) = TR2-TR5
        CH(I,3,K) = TI2+TI5
        CH(IC,2,K) = TI5-TI2
        CH(I-1,5,K) = TR3+TR4
        CH(IC-1,4,K) = TR3-TR4
        CH(I,5,K) = TI3+TI4
        CH(IC,4,K) = TI4-TI3
  105    CONTINUE
  106 CONTINUE
  return
end