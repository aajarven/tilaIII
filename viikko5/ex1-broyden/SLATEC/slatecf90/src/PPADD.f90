subroutine PPADD (N, IERROR, A, C, CBP, BP, BH)
!
!! PPADD is subsidiary to BLKTRI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PPADD-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   PPADD computes the eigenvalues of the periodic tridiagonal matrix
!   with coefficients AN,BN,CN.
!
!   N    is the order of the BH and BP polynomials.
!   BP   contains the eigenvalues on output.
!   CBP  is the same as BP except type complex.
!   BH   is used to temporarily store the roots of the B HAT polynomial
!        which enters through BP.
!
!***SEE ALSO  BLKTRI
!***ROUTINES CALLED  BSRH, PPSGF, PPSPF, PSGF
!***COMMON BLOCKS    CBLKT
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PPADD
!
  COMPLEX         CX         ,FSG        ,HSG        , &
                  DD         ,F          ,FP         ,FPP        , &
                  CDIS       ,R1         ,R2         ,R3         , &
                  CBP
  DIMENSION       A(*)       ,C(*)       ,BP(*)      ,BH(*)      , &
                  CBP(*)
  COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        , &
                  NM         ,NCMPLX     ,IK
  EXTERNAL        PSGF       ,PPSPF      ,PPSGF
!***FIRST EXECUTABLE STATEMENT  PPADD
  SCNV = SQRT(CNV)
  IZ = N
  if (BP(N)-BP(1)) 101,142,103
  101 DO 102 J=1,N
     NT = N-J
     BH(J) = BP(NT+1)
  102 CONTINUE
  go to 105
  103 DO 104 J=1,N
     BH(J) = BP(J)
  104 CONTINUE
  105 NCMPLX = 0
  MODIZ = MOD(IZ,2)
  IS = 1
  if (MODIZ) 106,107,106
  106 if (A(1)) 110,142,107
  107 XL = BH(1)
  DB = BH(3)-BH(1)
  108 XL = XL-DB
  if (PSGF(XL,IZ,C,A,BH)) 108,108,109
  109 SGN = -1.
  CBP(1) = CMPLX(BSRH(XL,BH(1),IZ,C,A,BH,PSGF,SGN),0.)
  IS = 2
  110 if = IZ-1
  if (MODIZ) 111,112,111
  111 if (A(1)) 112,142,115
  112 XR = BH(IZ)
  DB = BH(IZ)-BH(IZ-2)
  113 XR = XR+DB
  if (PSGF(XR,IZ,C,A,BH)) 113,114,114
  114 SGN = 1.
  CBP(IZ) = CMPLX(BSRH(BH(IZ),XR,IZ,C,A,BH,PSGF,SGN),0.)
  if = IZ-2
  115 DO 136 IG=IS,IF,2
     XL = BH(IG)
     XR = BH(IG+1)
     SGN = -1.
     XM = BSRH(XL,XR,IZ,C,A,BH,PPSPF,SGN)
     PSG = PSGF(XM,IZ,C,A,BH)
     if (ABS(PSG)-EPS) 118,118,116
  116    if (PSG*PPSGF(XM,IZ,C,A,BH)) 117,118,119
!
!     CASE OF A REAL ZERO
!
  117    SGN = 1.
     CBP(IG) = CMPLX(BSRH(BH(IG),XM,IZ,C,A,BH,PSGF,SGN),0.)
     SGN = -1.
     CBP(IG+1) = CMPLX(BSRH(XM,BH(IG+1),IZ,C,A,BH,PSGF,SGN),0.)
     go to 136
!
!     CASE OF A MULTIPLE ZERO
!
  118    CBP(IG) = CMPLX(XM,0.)
     CBP(IG+1) = CMPLX(XM,0.)
     go to 136
!
!     CASE OF A COMPLEX ZERO
!
  119    IT = 0
     ICV = 0
     CX = CMPLX(XM,0.)
  120    FSG = (1.,0.)
     HSG = (1.,0.)
     FP = (0.,0.)
     FPP = (0.,0.)
     DO 121 J=1,IZ
        DD = 1./(CX-BH(J))
        FSG = FSG*A(J)*DD
        HSG = HSG*C(J)*DD
        FP = FP+DD
        FPP = FPP-DD*DD
  121    CONTINUE
     if (MODIZ) 123,122,123
  122    F = (1.,0.)-FSG-HSG
     go to 124
  123    F = (1.,0.)+FSG+HSG
  124    I3 = 0
     if (ABS(FP)) 126,126,125
  125    I3 = 1
     R3 = -F/FP
  126    if (ABS(FPP)) 132,132,127
  127    CDIS = SQRT(FP**2-2.*F*FPP)
     R1 = CDIS-FP
     R2 = -FP-CDIS
     if (ABS(R1)-ABS(R2)) 129,129,128
  128    R1 = R1/FPP
     go to 130
  129    R1 = R2/FPP
  130    R2 = 2.*F/FPP/R1
     if (ABS(R2)  <  ABS(R1)) R1 = R2
     if (I3) 133,133,131
  131    if (ABS(R3)  <  ABS(R1)) R1 = R3
     go to 133
  132    R1 = R3
  133    CX = CX+R1
     IT = IT+1
     if (IT  >  50) go to 142
     if (ABS(R1)  >  SCNV) go to 120
     if (ICV) 134,134,135
  134    ICV = 1
     go to 120
  135    CBP(IG) = CX
     CBP(IG+1) = CONJG(CX)
  136 CONTINUE
  if (ABS(CBP(N))-ABS(CBP(1))) 137,142,139
  137 NHALF = N/2
  DO 138 J=1,NHALF
     NT = N-J
     CX = CBP(J)
     CBP(J) = CBP(NT+1)
     CBP(NT+1) = CX
  138 CONTINUE
  139 NCMPLX = 1
  DO 140 J=2,IZ
     if (AIMAG(CBP(J))) 143,140,143
  140 CONTINUE
  NCMPLX = 0
  DO 141 J=2,IZ
     BP(J) = REAL(CBP(J))
  141 CONTINUE
  go to 143
  142 IERROR = 4
  143 CONTINUE
  return
end
