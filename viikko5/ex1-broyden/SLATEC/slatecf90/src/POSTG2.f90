subroutine POSTG2 (NPEROD, N, M, A, BB, C, IDIMQ, Q, B, B2, B3, W, &
     W2, W3, D, TCOS, P)
!
!! POSTG2 is subsidiary to POISTG.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (POSTG2-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Subroutine to solve Poisson's equation on a staggered grid.
!
!***SEE ALSO  POISTG
!***ROUTINES CALLED  COSGEN, S1MERG, TRI3, TRIX
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   920130  Modified to use merge routine S1MERG rather than deleted
!           routine MERGE.  (WRB)
!***END PROLOGUE  POSTG2
!
  DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,*) , &
                  B(*)       ,B2(*)      ,B3(*)      ,W(*)       , &
                  W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    , &
                  K(4)       ,P(*)
  EQUIVALENCE     (K(1),K1)  ,(K(2),K2)  ,(K(3),K3)  ,(K(4),K4)
!***FIRST EXECUTABLE STATEMENT  POSTG2
  NP = NPEROD
  FNUM = 0.5*(NP/3)
  FNUM2 = 0.5*(NP/2)
  MR = M
  IP = -MR
  IPSTOR = 0
  I2R = 1
  JR = 2
  NR = N
  NLAST = N
  KR = 1
  LR = 0
  if (NR  <=  3) go to 142
  101 CONTINUE
  JR = 2*I2R
  NROD = 1
  if ((NR/2)*2  ==  NR) NROD = 0
  JSTART = 1
  JSTOP = NLAST-JR
  if (NROD  ==  0) JSTOP = JSTOP-I2R
  I2RBY2 = I2R/2
  if (JSTOP  >=  JSTART) go to 102
  J = JR
  go to 115
  102 CONTINUE
!
!     REGULAR REDUCTION.
!
  IJUMP = 1
  DO 114 J=JSTART,JSTOP,JR
     JP1 = J+I2RBY2
     JP2 = J+I2R
     JP3 = JP2+I2RBY2
     JM1 = J-I2RBY2
     JM2 = J-I2R
     JM3 = JM2-I2RBY2
     if (J  /=  1) go to 106
     call COSGEN (I2R,1,FNUM,0.5,TCOS)
     if (I2R  /=  1) go to 104
     DO 103 I=1,MR
        B(I) = Q(I,1)
        Q(I,1) = Q(I,2)
  103    CONTINUE
     go to 112
  104    DO 105 I=1,MR
        B(I) = Q(I,1)+0.5*(Q(I,JP2)-Q(I,JP1)-Q(I,JP3))
        Q(I,1) = Q(I,JP2)+Q(I,1)-Q(I,JP1)
  105    CONTINUE
     go to 112
  106    CONTINUE
     go to (107,108),IJUMP
  107    CONTINUE
     IJUMP = 2
     call COSGEN (I2R,1,0.5,0.0,TCOS)
  108    CONTINUE
     if (I2R  /=  1) go to 110
     DO 109 I=1,MR
        B(I) = 2.*Q(I,J)
        Q(I,J) = Q(I,JM2)+Q(I,JP2)
  109    CONTINUE
     go to 112
  110    DO 111 I=1,MR
        FI = Q(I,J)
        Q(I,J) = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
        B(I) = FI+Q(I,J)-Q(I,JM3)-Q(I,JP3)
  111    CONTINUE
  112    CONTINUE
     call TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
     DO 113 I=1,MR
        Q(I,J) = Q(I,J)+B(I)
  113    CONTINUE
!
!     END OF REDUCTION FOR REGULAR UNKNOWNS.
!
  114 CONTINUE
!
!     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
!
  J = JSTOP+JR
  115 NLAST = J
  JM1 = J-I2RBY2
  JM2 = J-I2R
  JM3 = JM2-I2RBY2
  if (NROD  ==  0) go to 125
!
!     ODD NUMBER OF UNKNOWNS
!
  if (I2R  /=  1) go to 117
  DO 116 I=1,MR
     B(I) = Q(I,J)
     Q(I,J) = Q(I,JM2)
  116 CONTINUE
  go to 123
  117 DO 118 I=1,MR
     B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  118 CONTINUE
  if (NRODPR  /=  0) go to 120
  DO 119 I=1,MR
     II = IP+I
     Q(I,J) = Q(I,JM2)+P(II)
  119 CONTINUE
  IP = IP-MR
  go to 122
  120 CONTINUE
  DO 121 I=1,MR
     Q(I,J) = Q(I,J)-Q(I,JM1)+Q(I,JM2)
  121 CONTINUE
  122 if (LR  ==  0) go to 123
  call COSGEN (LR,1,FNUM2,0.5,TCOS(KR+1))
  123 CONTINUE
  call COSGEN (KR,1,FNUM2,0.5,TCOS)
  call TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
  DO 124 I=1,MR
     Q(I,J) = Q(I,J)+B(I)
  124 CONTINUE
  KR = KR+I2R
  go to 141
  125 CONTINUE
!
!     EVEN NUMBER OF UNKNOWNS
!
  JP1 = J+I2RBY2
  JP2 = J+I2R
  if (I2R  /=  1) go to 129
  DO 126 I=1,MR
     B(I) = Q(I,J)
  126 CONTINUE
  TCOS(1) = 0.
  call TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  IP = 0
  IPSTOR = MR
  DO 127 I=1,MR
     P(I) = B(I)
     B(I) = B(I)+Q(I,N)
  127 CONTINUE
  TCOS(1) = -1.+2*(NP/2)
  TCOS(2) = 0.
  call TRIX (1,1,MR,A,BB,C,B,TCOS,D,W)
  DO 128 I=1,MR
     Q(I,J) = Q(I,JM2)+P(I)+B(I)
  128 CONTINUE
  go to 140
  129 CONTINUE
  DO 130 I=1,MR
     B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  130 CONTINUE
  if (NRODPR  /=  0) go to 132
  DO 131 I=1,MR
     II = IP+I
     B(I) = B(I)+P(II)
  131 CONTINUE
  go to 134
  132 CONTINUE
  DO 133 I=1,MR
     B(I) = B(I)+Q(I,JP2)-Q(I,JP1)
  133 CONTINUE
  134 CONTINUE
  call COSGEN (I2R,1,0.5,0.0,TCOS)
  call TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
  IP = IP+MR
  IPSTOR = MAX(IPSTOR,IP+MR)
  DO 135 I=1,MR
     II = IP+I
     P(II) = B(I)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
     B(I) = P(II)+Q(I,JP2)
  135 CONTINUE
  if (LR  ==  0) go to 136
  call COSGEN (LR,1,FNUM2,0.5,TCOS(I2R+1))
  call S1MERG (TCOS,0,I2R,I2R,LR,KR)
  go to 138
  136 DO 137 I=1,I2R
     II = KR+I
     TCOS(II) = TCOS(I)
  137 CONTINUE
  138 call COSGEN (KR,1,FNUM2,0.5,TCOS)
  call TRIX (KR,KR,MR,A,BB,C,B,TCOS,D,W)
  DO 139 I=1,MR
     II = IP+I
     Q(I,J) = Q(I,JM2)+P(II)+B(I)
  139 CONTINUE
  140 CONTINUE
  LR = KR
  KR = KR+JR
  141 CONTINUE
  NR = (NLAST-1)/JR+1
  if (NR  <=  3) go to 142
  I2R = JR
  NRODPR = NROD
  go to 101
  142 CONTINUE
!
!      BEGIN SOLUTION
!
  J = 1+JR
  JM1 = J-I2R
  JP1 = J+I2R
  JM2 = NLAST-I2R
  if (NR  ==  2) go to 180
  if (LR  /=  0) go to 167
  if (N  /=  3) go to 156
!
!     CASE N = 3.
!
  go to (143,148,143),NP
  143 DO 144 I=1,MR
     B(I) = Q(I,2)
     B2(I) = Q(I,1)+Q(I,3)
     B3(I) = 0.
  144 CONTINUE
  go to (146,146,145),NP
  145 TCOS(1) = -1.
  TCOS(2) = 1.
  K1 = 1
  go to 147
  146 TCOS(1) = -2.
  TCOS(2) = 1.
  TCOS(3) = -1.
  K1 = 2
  147 K2 = 1
  K3 = 0
  K4 = 0
  go to 150
  148 DO 149 I=1,MR
     B(I) = Q(I,2)
     B2(I) = Q(I,3)
     B3(I) = Q(I,1)
  149 CONTINUE
  call COSGEN (3,1,0.5,0.0,TCOS)
  TCOS(4) = -1.
  TCOS(5) = 1.
  TCOS(6) = -1.
  TCOS(7) = 1.
  K1 = 3
  K2 = 2
  K3 = 1
  K4 = 1
  150 call TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
  DO 151 I=1,MR
     B(I) = B(I)+B2(I)+B3(I)
  151 CONTINUE
  go to (153,153,152),NP
  152 TCOS(1) = 2.
  call TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  153 DO 154 I=1,MR
     Q(I,2) = B(I)
     B(I) = Q(I,1)+B(I)
  154 CONTINUE
  TCOS(1) = -1.+4.*FNUM
  call TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  DO 155 I=1,MR
     Q(I,1) = B(I)
  155 CONTINUE
  JR = 1
  I2R = 0
  go to 188
!
!     CASE N = 2**P+1
!
  156 CONTINUE
  DO 157 I=1,MR
     B(I) = Q(I,J)+Q(I,1)-Q(I,JM1)+Q(I,NLAST)-Q(I,JM2)
  157 CONTINUE
  go to (158,160,158),NP
  158 DO 159 I=1,MR
     B2(I) = Q(I,1)+Q(I,NLAST)+Q(I,J)-Q(I,JM1)-Q(I,JP1)
     B3(I) = 0.
  159 CONTINUE
  K1 = NLAST-1
  K2 = NLAST+JR-1
  call COSGEN (JR-1,1,0.0,1.0,TCOS(NLAST))
  TCOS(K2) = 2*NP-4
  call COSGEN (JR,1,0.5-FNUM,0.5,TCOS(K2+1))
  K3 = (3-NP)/2
  call S1MERG (TCOS,K1,JR-K3,K2-K3,JR+K3,0)
  K1 = K1-1+K3
  call COSGEN (JR,1,FNUM,0.5,TCOS(K1+1))
  K2 = JR
  K3 = 0
  K4 = 0
  go to 162
  160 DO 161 I=1,MR
     FI = (Q(I,J)-Q(I,JM1)-Q(I,JP1))/2.
     B2(I) = Q(I,1)+FI
     B3(I) = Q(I,NLAST)+FI
  161 CONTINUE
  K1 = NLAST+JR-1
  K2 = K1+JR-1
  call COSGEN (JR-1,1,0.0,1.0,TCOS(K1+1))
  call COSGEN (NLAST,1,0.5,0.0,TCOS(K2+1))
  call S1MERG (TCOS,K1,JR-1,K2,NLAST,0)
  K3 = K1+NLAST-1
  K4 = K3+JR
  call COSGEN (JR,1,0.5,0.5,TCOS(K3+1))
  call COSGEN (JR,1,0.0,0.5,TCOS(K4+1))
  call S1MERG (TCOS,K3,JR,K4,JR,K1)
  K2 = NLAST-1
  K3 = JR
  K4 = JR
  162 call TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
  DO 163 I=1,MR
     B(I) = B(I)+B2(I)+B3(I)
  163 CONTINUE
  if (NP  /=  3) go to 164
  TCOS(1) = 2.
  call TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  164 DO 165 I=1,MR
     Q(I,J) = B(I)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
     B(I) = Q(I,J)+Q(I,1)
  165 CONTINUE
  call COSGEN (JR,1,FNUM,0.5,TCOS)
  call TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
  DO 166 I=1,MR
     Q(I,1) = Q(I,1)-Q(I,JM1)+B(I)
  166 CONTINUE
  go to 188
!
!     CASE OF GENERAL N WITH NR = 3 .
!
  167 CONTINUE
  DO 168 I=1,MR
     B(I) = Q(I,1)-Q(I,JM1)+Q(I,J)
  168 CONTINUE
  if (NROD  /=  0) go to 170
  DO 169 I=1,MR
     II = IP+I
     B(I) = B(I)+P(II)
  169 CONTINUE
  go to 172
  170 DO 171 I=1,MR
     B(I) = B(I)+Q(I,NLAST)-Q(I,JM2)
  171 CONTINUE
  172 CONTINUE
  DO 173 I=1,MR
     T = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
     Q(I,J) = T
     B2(I) = Q(I,NLAST)+T
     B3(I) = Q(I,1)+T
  173 CONTINUE
  K1 = KR+2*JR
  call COSGEN (JR-1,1,0.0,1.0,TCOS(K1+1))
  K2 = K1+JR
  TCOS(K2) = 2*NP-4
  K4 = (NP-1)*(3-NP)
  K3 = K2+1-K4
  call COSGEN (KR+JR+K4,1,K4/2.,1.-K4,TCOS(K3))
  K4 = 1-NP/3
  call S1MERG (TCOS,K1,JR-K4,K2-K4,KR+JR+K4,0)
  if (NP  ==  3) K1 = K1-1
  K2 = KR+JR
  K4 = K1+K2
  call COSGEN (KR,1,FNUM2,0.5,TCOS(K4+1))
  K3 = K4+KR
  call COSGEN (JR,1,FNUM,0.5,TCOS(K3+1))
  call S1MERG (TCOS,K4,KR,K3,JR,K1)
  K4 = K3+JR
  call COSGEN (LR,1,FNUM2,0.5,TCOS(K4+1))
  call S1MERG (TCOS,K3,JR,K4,LR,K1+K2)
  call COSGEN (KR,1,FNUM2,0.5,TCOS(K3+1))
  K3 = KR
  K4 = KR
  call TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
  DO 174 I=1,MR
     B(I) = B(I)+B2(I)+B3(I)
  174 CONTINUE
  if (NP  /=  3) go to 175
  TCOS(1) = 2.
  call TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  175 DO 176 I=1,MR
     Q(I,J) = Q(I,J)+B(I)
     B(I) = Q(I,1)+Q(I,J)
  176 CONTINUE
  call COSGEN (JR,1,FNUM,0.5,TCOS)
  call TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
  if (JR  /=  1) go to 178
  DO 177 I=1,MR
     Q(I,1) = B(I)
  177 CONTINUE
  go to 188
  178 CONTINUE
  DO 179 I=1,MR
     Q(I,1) = Q(I,1)-Q(I,JM1)+B(I)
  179 CONTINUE
  go to 188
  180 CONTINUE
!
!     CASE OF GENERAL N AND NR = 2 .
!
  DO 181 I=1,MR
     II = IP+I
     B3(I) = 0.
     B(I) = Q(I,1)+P(II)
     Q(I,1) = Q(I,1)-Q(I,JM1)
     B2(I) = Q(I,1)+Q(I,NLAST)
  181 CONTINUE
  K1 = KR+JR
  K2 = K1+JR
  call COSGEN (JR-1,1,0.0,1.0,TCOS(K1+1))
  go to (182,183,182),NP
  182 TCOS(K2) = 2*NP-4
  call COSGEN (KR,1,0.0,1.0,TCOS(K2+1))
  go to 184
  183 call COSGEN (KR+1,1,0.5,0.0,TCOS(K2))
  184 K4 = 1-NP/3
  call S1MERG (TCOS,K1,JR-K4,K2-K4,KR+K4,0)
  if (NP  ==  3) K1 = K1-1
  K2 = KR
  call COSGEN (KR,1,FNUM2,0.5,TCOS(K1+1))
  K4 = K1+KR
  call COSGEN (LR,1,FNUM2,0.5,TCOS(K4+1))
  K3 = LR
  K4 = 0
  call TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
  DO 185 I=1,MR
     B(I) = B(I)+B2(I)
  185 CONTINUE
  if (NP  /=  3) go to 186
  TCOS(1) = 2.
  call TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  186 DO 187 I=1,MR
     Q(I,1) = Q(I,1)+B(I)
  187 CONTINUE
  188 CONTINUE
!
!     START BACK SUBSTITUTION.
!
  J = NLAST-JR
  DO 189 I=1,MR
     B(I) = Q(I,NLAST)+Q(I,J)
  189 CONTINUE
  JM2 = NLAST-I2R
  if (JR  /=  1) go to 191
  DO 190 I=1,MR
     Q(I,NLAST) = 0.
  190 CONTINUE
  go to 195
  191 CONTINUE
  if (NROD  /=  0) go to 193
  DO 192 I=1,MR
     II = IP+I
     Q(I,NLAST) = P(II)
  192 CONTINUE
  IP = IP-MR
  go to 195
  193 DO 194 I=1,MR
     Q(I,NLAST) = Q(I,NLAST)-Q(I,JM2)
  194 CONTINUE
  195 CONTINUE
  call COSGEN (KR,1,FNUM2,0.5,TCOS)
  call COSGEN (LR,1,FNUM2,0.5,TCOS(KR+1))
  call TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
  DO 196 I=1,MR
     Q(I,NLAST) = Q(I,NLAST)+B(I)
  196 CONTINUE
  NLASTP = NLAST
  197 CONTINUE
  JSTEP = JR
  JR = I2R
  I2R = I2R/2
  if (JR  ==  0) go to 210
  JSTART = 1+JR
  KR = KR-JR
  if (NLAST+JR  >  N) go to 198
  KR = KR-JR
  NLAST = NLAST+JR
  JSTOP = NLAST-JSTEP
  go to 199
  198 CONTINUE
  JSTOP = NLAST-JR
  199 CONTINUE
  LR = KR-JR
  call COSGEN (JR,1,0.5,0.0,TCOS)
  DO 209 J=JSTART,JSTOP,JSTEP
     JM2 = J-JR
     JP2 = J+JR
     if (J  /=  JR) go to 201
     DO 200 I=1,MR
        B(I) = Q(I,J)+Q(I,JP2)
  200    CONTINUE
     go to 203
  201    CONTINUE
     DO 202 I=1,MR
        B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  202    CONTINUE
  203    CONTINUE
     if (JR  /=  1) go to 205
     DO 204 I=1,MR
        Q(I,J) = 0.
  204    CONTINUE
     go to 207
  205    CONTINUE
     JM1 = J-I2R
     JP1 = J+I2R
     DO 206 I=1,MR
        Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  206    CONTINUE
  207    CONTINUE
     call TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
     DO 208 I=1,MR
        Q(I,J) = Q(I,J)+B(I)
  208    CONTINUE
  209 CONTINUE
  NROD = 1
  if (NLAST+I2R  <=  N) NROD = 0
  if (NLASTP  /=  NLAST) go to 188
  go to 197
  210 CONTINUE
!
!     return STORAGE REQUIREMENTS FOR P VECTORS.
!
  W(1) = IPSTOR
  return
end
