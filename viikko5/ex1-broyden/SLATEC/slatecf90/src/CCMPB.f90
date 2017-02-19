subroutine CCMPB (N, IERROR, AN, BN, CN, B, AH, BH)
!
!! CCMPB is subsidiary to CBLKTR.
!
!***LIBRARY   SLATEC
!***TYPE      COMPLEX (COMPB-S, CCMPB-C)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     CCMPB computes the roots of the B polynomials using subroutine
!     TEVLC which is a modification the EISPACK program TQLRAT.
!     IERROR is set to 4 if either TEVLC fails or if A(J+1)*C(J) is
!     less than zero for some J.  AH,BH are temporary work arrays.
!
!***SEE ALSO  CBLKTR
!***ROUTINES CALLED  CPADD, INXCB, R1MACH, TEVLC
!***COMMON BLOCKS    CCBLK
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  CCMPB
!
  DIMENSION       AN(*)      ,BN(*)      ,CN(*)      ,B(*)       , &
                  AH(*)      ,BH(*)
  COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        , &
                  NM         ,NCMPLX     ,IK
!***FIRST EXECUTABLE STATEMENT  CCMPB
  EPS = R1MACH(4)
  BNORM = ABS(BN(1))
  DO 102 J=2,NM
     BNORM = MAX(BNORM,ABS(BN(J)))
     ARG = AN(J)*CN(J-1)
     if (ARG) 119,101,101
  101    B(J) = SIGN(SQRT(ARG),AN(J))
  102 CONTINUE
  CNV = EPS*BNORM
  if = 2**K
  KDO = K-1
  DO 108 L=1,KDO
     IR = L-1
     I2 = 2**IR
     I4 = I2+I2
     IPL = I4-1
     IFD = IF-I4
     DO 107 I=I4,IFD,I4
        call INXCB (I,L,IB,NB)
        if (NB) 108,108,103
  103       JS = I-IPL
        JF = JS+NB-1
        LS = 0
        DO 104 J=JS,JF
           LS = LS+1
           BH(LS) = BN(J)
           AH(LS) = B(J)
  104       CONTINUE
        call TEVLC (NB,BH,AH,IERROR)
        if (IERROR) 118,105,118
  105       LH = IB-1
        DO 106 J=1,NB
           LH = LH+1
           B(LH) = -BH(J)
  106       CONTINUE
  107    CONTINUE
  108 CONTINUE
  DO 109 J=1,NM
     B(J) = -BN(J)
  109 CONTINUE
  if (NPP) 117,110,117
  110 NMP = NM+1
  NB = NM+NMP
  DO 112 J=1,NB
     L1 = MOD(J-1,NMP)+1
     L2 = MOD(J+NM-1,NMP)+1
     ARG = AN(L1)*CN(L2)
     if (ARG) 119,111,111
  111    BH(J) = SIGN(SQRT(ARG),-AN(L1))
     AH(J) = -BN(L1)
  112 CONTINUE
  call TEVLC (NB,AH,BH,IERROR)
  if (IERROR) 118,113,118
  113 call INXCB (IF,K-1,J2,LH)
  call INXCB (IF/2,K-1,J1,LH)
  J2 = J2+1
  LH = J2
  N2M2 = J2+NM+NM-2
  114 D1 = ABS(B(J1)-B(J2-1))
  D2 = ABS(B(J1)-B(J2))
  D3 = ABS(B(J1)-B(J2+1))
  if ((D2  <  D1) .AND. (D2  <  D3)) go to 115
  B(LH) = B(J2)
  J2 = J2+1
  LH = LH+1
  if (J2-N2M2) 114,114,116
  115 J2 = J2+1
  J1 = J1+1
  if (J2-N2M2) 114,114,116
  116 B(LH) = B(N2M2+1)
  call INXCB (IF,K-1,J1,J2)
  J2 = J1+NMP+NMP
  call CPADD (NM+1,IERROR,AN,CN,B(J1),B(J1),B(J2))
  117 RETURN
  118 IERROR = 4
  return
  119 IERROR = 5
  return
end
