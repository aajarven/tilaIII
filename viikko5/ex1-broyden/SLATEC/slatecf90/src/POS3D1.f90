subroutine POS3D1 (LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, XRT, &
     YRT, T, D, WX, WY, C1, C2, BB)
!
!! POS3D1 is subsidiary to POIS3D.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (POS3D1-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  POIS3D
!***ROUTINES CALLED  COSQB, COSQF, COSQI, COST, COSTI, PIMACH, RFFTB,
!                    RFFTF, RFFTI, SINQB, SINQF, SINQI, SINT, SINTI,
!                    TRIDQ
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891009  Removed unreferenced variable.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900308  Changed call to TRID to call to TRIDQ.  (WRB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  POS3D1
  DIMENSION       A(*)       ,B(*)       ,C(*)       , &
                  F(LDIMF,MDIMF,*)       ,XRT(*)     ,YRT(*)     , &
                  T(*)       ,D(*)       ,WX(*)      ,WY(*)      , &
                  BB(*)
!***FIRST EXECUTABLE STATEMENT  POS3D1
  PI = PIMACH(DUM)
  LR = L
  MR = M
  NR = N
!
!     GENERATE TRANSFORM ROOTS
!
  LRDEL = ((LP-1)*(LP-3)*(LP-5))/3
  SCALX = LR+LRDEL
  DX = PI/(2.*SCALX)
  go to (108,103,101,102,101),LP
  101 DI = 0.5
  SCALX = 2.*SCALX
  go to 104
  102 DI = 1.0
  go to 104
  103 DI = 0.0
  104 DO 105 I=1,LR
     XRT(I) = -4.*C1*(SIN((I-DI)*DX))**2
  105 CONTINUE
  SCALX = 2.*SCALX
  go to (112,106,110,107,111),LP
  106 call SINTI (LR,WX)
  go to 112
  107 call COSTI (LR,WX)
  go to 112
  108 XRT(1) = 0.
  XRT(LR) = -4.*C1
  DO 109 I=3,LR,2
     XRT(I-1) = -4.*C1*(SIN((I-1)*DX))**2
     XRT(I) = XRT(I-1)
  109 CONTINUE
  call RFFTI (LR,WX)
  go to 112
  110 call SINQI (LR,WX)
  go to 112
  111 call COSQI (LR,WX)
  112 CONTINUE
  MRDEL = ((MP-1)*(MP-3)*(MP-5))/3
  SCALY = MR+MRDEL
  DY = PI/(2.*SCALY)
  go to (120,115,113,114,113),MP
  113 DJ = 0.5
  SCALY = 2.*SCALY
  go to 116
  114 DJ = 1.0
  go to 116
  115 DJ = 0.0
  116 DO 117 J=1,MR
     YRT(J) = -4.*C2*(SIN((J-DJ)*DY))**2
  117 CONTINUE
  SCALY = 2.*SCALY
  go to (124,118,122,119,123),MP
  118 call SINTI (MR,WY)
  go to 124
  119 call COSTI (MR,WY)
  go to 124
  120 YRT(1) = 0.
  YRT(MR) = -4.*C2
  DO 121 J=3,MR,2
     YRT(J-1) = -4.*C2*(SIN((J-1)*DY))**2
     YRT(J) = YRT(J-1)
  121 CONTINUE
  call RFFTI (MR,WY)
  go to 124
  122 call SINQI (MR,WY)
  go to 124
  123 call COSQI (MR,WY)
  124 CONTINUE
  IFWRD = 1
  125 CONTINUE
!
!     TRANSFORM X
!
  DO 141 J=1,MR
     DO 140 K=1,NR
        DO 126 I=1,LR
           T(I) = F(I,J,K)
  126       CONTINUE
        go to (127,130,131,134,135),LP
  127       go to (128,129),IFWRD
  128       call RFFTF (LR,T,WX)
        go to 138
  129       call RFFTB (LR,T,WX)
        go to 138
  130       call SINT (LR,T,WX)
        go to 138
  131       go to (132,133),IFWRD
  132       call SINQF (LR,T,WX)
        go to 138
  133       call SINQB (LR,T,WX)
        go to 138
  134       call COST (LR,T,WX)
        go to 138
  135       go to (136,137),IFWRD
  136       call COSQF (LR,T,WX)
        go to 138
  137       call COSQB (LR,T,WX)
  138       CONTINUE
        DO 139 I=1,LR
           F(I,J,K) = T(I)
  139       CONTINUE
  140    CONTINUE
  141 CONTINUE
  go to (142,164),IFWRD
!
!     TRANSFORM Y
!
  142 CONTINUE
  DO 158 I=1,LR
     DO 157 K=1,NR
        DO 143 J=1,MR
           T(J) = F(I,J,K)
  143       CONTINUE
        go to (144,147,148,151,152),MP
  144       go to (145,146),IFWRD
  145       call RFFTF (MR,T,WY)
        go to 155
  146       call RFFTB (MR,T,WY)
        go to 155
  147       call SINT (MR,T,WY)
        go to 155
  148       go to (149,150),IFWRD
  149       call SINQF (MR,T,WY)
        go to 155
  150       call SINQB (MR,T,WY)
        go to 155
  151       call COST (MR,T,WY)
        go to 155
  152       go to (153,154),IFWRD
  153       call COSQF (MR,T,WY)
        go to 155
  154       call COSQB (MR,T,WY)
  155       CONTINUE
        DO 156 J=1,MR
           F(I,J,K) = T(J)
  156       CONTINUE
  157    CONTINUE
  158 CONTINUE
  go to (159,125),IFWRD
  159 CONTINUE
!
!     SOLVE TRIDIAGONAL SYSTEMS IN Z
!
  DO 163 I=1,LR
     DO 162 J=1,MR
        DO 160 K=1,NR
           BB(K) = B(K)+XRT(I)+YRT(J)
           T(K) = F(I,J,K)
  160       CONTINUE
        call TRIDQ (NR,A,BB,C,T,D)
        DO 161 K=1,NR
           F(I,J,K) = T(K)
  161       CONTINUE
  162    CONTINUE
  163 CONTINUE
  IFWRD = 2
  go to 142
  164 CONTINUE
  DO 167 I=1,LR
     DO 166 J=1,MR
        DO 165 K=1,NR
           F(I,J,K) = F(I,J,K)/(SCALX*SCALY)
  165       CONTINUE
  166    CONTINUE
  167 CONTINUE
  return
end
