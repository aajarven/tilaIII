subroutine PASSF (NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
!
!! PASSF calculates fast Fourier transforms of subvectors of arbitrary length.
!
!***LIBRARY   SLATEC (FFTPACK)
!***TYPE      SINGLE PRECISION (PASSF-S)
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*).
!   881128  Modified by Dick Valent to meet prologue standards.
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced variable.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PASSF
  DIMENSION CH(IDO,L1,*), CC(IDO,IP,*), C1(IDO,L1,*), WA(*), &
            C2(IDL1,*), CH2(IDL1,*)
!***FIRST EXECUTABLE STATEMENT  PASSF
  IDOT = IDO/2
  IPP2 = IP+2
  IPPH = (IP+1)/2
  IDP = IP*IDO
!
  if (IDO  <  L1) go to 106
  DO 103 J=2,IPPH
     JC = IPP2-J
     DO 102 K=1,L1
!DIR$ IVDEP
        DO 101 I=1,IDO
           CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
           CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
  DO 105 K=1,L1
!DIR$ IVDEP
     DO 104 I=1,IDO
        CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
  go to 112
  106 DO 109 J=2,IPPH
     JC = IPP2-J
     DO 108 I=1,IDO
!DIR$ IVDEP
        DO 107 K=1,L1
           CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
           CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
  DO 111 I=1,IDO
!DIR$ IVDEP
     DO 110 K=1,L1
        CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
  INC = 0
  DO 116 L=2,IPPH
     LC = IPP2-L
     IDL = IDL+IDO
!DIR$ IVDEP
     DO 113 IK=1,IDL1
        C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
        C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
  113    CONTINUE
     IDLJ = IDL
     INC = INC+IDO
     DO 115 J=3,IPPH
        JC = IPP2-J
        IDLJ = IDLJ+INC
        if (IDLJ  >  IDP) IDLJ = IDLJ-IDP
        WAR = WA(IDLJ-1)
        WAI = WA(IDLJ)
!DIR$ IVDEP
        DO 114 IK=1,IDL1
           C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
           C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
  DO 118 J=2,IPPH
!DIR$ IVDEP
     DO 117 IK=1,IDL1
        CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
  DO 120 J=2,IPPH
     JC = IPP2-J
!DIR$ IVDEP
     DO 119 IK=2,IDL1,2
        CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
        CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
        CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
        CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
  NAC = 1
  if (IDO  ==  2) RETURN
  NAC = 0
!DIR$ IVDEP
  DO 121 IK=1,IDL1
     C2(IK,1) = CH2(IK,1)
  121 CONTINUE
  DO 123 J=2,IP
!DIR$ IVDEP
     DO 122 K=1,L1
        C1(1,K,J) = CH(1,K,J)
        C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
  if (IDOT  >  L1) go to 127
  IDIJ = 0
  DO 126 J=2,IP
     IDIJ = IDIJ+2
     DO 125 I=4,IDO,2
        IDIJ = IDIJ+2
!DIR$ IVDEP
        DO 124 K=1,L1
           C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
           C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
  return
  127 IDJ = 2-IDO
  DO 130 J=2,IP
     IDJ = IDJ+IDO
     DO 129 K=1,L1
        IDIJ = IDJ
!DIR$ IVDEP
        DO 128 I=4,IDO,2
           IDIJ = IDIJ+2
           C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
           C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
  return
end
