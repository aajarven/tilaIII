subroutine POISN2 (M, N, ISTAG, MIXBND, A, BB, C, Q, IDIMQ, B, B2, &
     B3, W, W2, W3, D, TCOS, P)
!
!! POISN2 is subsidiary to GENBUN.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (POISN2-S, CMPOSN-C)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Subroutine to solve Poisson's equation with Neumann boundary
!     conditions.
!
!     ISTAG = 1 if the last diagonal block is A.
!     ISTAG = 2 if the last diagonal block is A-I.
!     MIXBND = 1 if have Neumann boundary conditions at both boundaries.
!     MIXBND = 2 if have Neumann boundary conditions at bottom and
!     Dirichlet condition at top.  (for this case, must have ISTAG = 1.)
!
!***SEE ALSO  GENBUN
!***ROUTINES CALLED  COSGEN, S1MERG, TRI3, TRIX
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   920130  Modified to use merge routine S1MERG rather than deleted
!           routine MERGE.  (WRB)
!***END PROLOGUE  POISN2
!
  DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,*) , &
                  B(*)       ,B2(*)      ,B3(*)      ,W(*)       , &
                  W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    , &
                  K(4)       ,P(*)
  EQUIVALENCE     (K(1),K1)  ,(K(2),K2)  ,(K(3),K3)  ,(K(4),K4)
!***FIRST EXECUTABLE STATEMENT  POISN2
  FISTAG = 3-ISTAG
  FNUM = 1./ISTAG
  FDEN = 0.5*(ISTAG-1)
  MR = M
  IP = -MR
  IPSTOR = 0
  I2R = 1
  JR = 2
  NR = N
  NLAST = N
  KR = 1
  LR = 0
  go to (101,103),ISTAG
  101 CONTINUE
  DO 102 I=1,MR
     Q(I,N) = .5*Q(I,N)
  102 CONTINUE
  go to (103,104),MIXBND
  103 if (N  <=  3) go to 155
  104 CONTINUE
  JR = 2*I2R
  NROD = 1
  if ((NR/2)*2  ==  NR) NROD = 0
  go to (105,106),MIXBND
  105 JSTART = 1
  go to 107
  106 JSTART = JR
  NROD = 1-NROD
  107 CONTINUE
  JSTOP = NLAST-JR
  if (NROD  ==  0) JSTOP = JSTOP-I2R
  call COSGEN (I2R,1,0.5,0.0,TCOS)
  I2RBY2 = I2R/2
  if (JSTOP  >=  JSTART) go to 108
  J = JR
  go to 116
  108 CONTINUE
!
!     REGULAR REDUCTION.
!
  DO 115 J=JSTART,JSTOP,JR
     JP1 = J+I2RBY2
     JP2 = J+I2R
     JP3 = JP2+I2RBY2
     JM1 = J-I2RBY2
     JM2 = J-I2R
     JM3 = JM2-I2RBY2
     if (J  /=  1) go to 109
     JM1 = JP1
     JM2 = JP2
     JM3 = JP3
  109    CONTINUE
     if (I2R  /=  1) go to 111
     if (J  ==  1) JM2 = JP2
     DO 110 I=1,MR
        B(I) = 2.*Q(I,J)
        Q(I,J) = Q(I,JM2)+Q(I,JP2)
  110    CONTINUE
     go to 113
  111    CONTINUE
     DO 112 I=1,MR
        FI = Q(I,J)
        Q(I,J) = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
        B(I) = FI+Q(I,J)-Q(I,JM3)-Q(I,JP3)
  112    CONTINUE
  113    CONTINUE
     call TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
     DO 114 I=1,MR
        Q(I,J) = Q(I,J)+B(I)
  114    CONTINUE
!
!     END OF REDUCTION FOR REGULAR UNKNOWNS.
!
  115 CONTINUE
!
!     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
!
  J = JSTOP+JR
  116 NLAST = J
  JM1 = J-I2RBY2
  JM2 = J-I2R
  JM3 = JM2-I2RBY2
  if (NROD  ==  0) go to 128
!
!     ODD NUMBER OF UNKNOWNS
!
  if (I2R  /=  1) go to 118
  DO 117 I=1,MR
     B(I) = FISTAG*Q(I,J)
     Q(I,J) = Q(I,JM2)
  117 CONTINUE
  go to 126
  118 DO 119 I=1,MR
     B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  119 CONTINUE
  if (NRODPR  /=  0) go to 121
  DO 120 I=1,MR
     II = IP+I
     Q(I,J) = Q(I,JM2)+P(II)
  120 CONTINUE
  IP = IP-MR
  go to 123
  121 CONTINUE
  DO 122 I=1,MR
     Q(I,J) = Q(I,J)-Q(I,JM1)+Q(I,JM2)
  122 CONTINUE
  123 if (LR  ==  0) go to 124
  call COSGEN (LR,1,0.5,FDEN,TCOS(KR+1))
  go to 126
  124 CONTINUE
  DO 125 I=1,MR
     B(I) = FISTAG*B(I)
  125 CONTINUE
  126 CONTINUE
  call COSGEN (KR,1,0.5,FDEN,TCOS)
  call TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
  DO 127 I=1,MR
     Q(I,J) = Q(I,J)+B(I)
  127 CONTINUE
  KR = KR+I2R
  go to 151
  128 CONTINUE
!
!     EVEN NUMBER OF UNKNOWNS
!
  JP1 = J+I2RBY2
  JP2 = J+I2R
  if (I2R  /=  1) go to 135
  DO 129 I=1,MR
     B(I) = Q(I,J)
  129 CONTINUE
  call TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  IP = 0
  IPSTOR = MR
  go to (133,130),ISTAG
  130 DO 131 I=1,MR
     P(I) = B(I)
     B(I) = B(I)+Q(I,N)
  131 CONTINUE
  TCOS(1) = 1.
  TCOS(2) = 0.
  call TRIX (1,1,MR,A,BB,C,B,TCOS,D,W)
  DO 132 I=1,MR
     Q(I,J) = Q(I,JM2)+P(I)+B(I)
  132 CONTINUE
  go to 150
  133 CONTINUE
  DO 134 I=1,MR
     P(I) = B(I)
     Q(I,J) = Q(I,JM2)+2.*Q(I,JP2)+3.*B(I)
  134 CONTINUE
  go to 150
  135 CONTINUE
  DO 136 I=1,MR
     B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  136 CONTINUE
  if (NRODPR  /=  0) go to 138
  DO 137 I=1,MR
     II = IP+I
     B(I) = B(I)+P(II)
  137 CONTINUE
  go to 140
  138 CONTINUE
  DO 139 I=1,MR
     B(I) = B(I)+Q(I,JP2)-Q(I,JP1)
  139 CONTINUE
  140 CONTINUE
  call TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
  IP = IP+MR
  IPSTOR = MAX(IPSTOR,IP+MR)
  DO 141 I=1,MR
     II = IP+I
     P(II) = B(I)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
     B(I) = P(II)+Q(I,JP2)
  141 CONTINUE
  if (LR  ==  0) go to 142
  call COSGEN (LR,1,0.5,FDEN,TCOS(I2R+1))
  call S1MERG (TCOS,0,I2R,I2R,LR,KR)
  go to 144
  142 DO 143 I=1,I2R
     II = KR+I
     TCOS(II) = TCOS(I)
  143 CONTINUE
  144 call COSGEN (KR,1,0.5,FDEN,TCOS)
  if (LR  /=  0) go to 145
  go to (146,145),ISTAG
  145 CONTINUE
  call TRIX (KR,KR,MR,A,BB,C,B,TCOS,D,W)
  go to 148
  146 CONTINUE
  DO 147 I=1,MR
     B(I) = FISTAG*B(I)
  147 CONTINUE
  148 CONTINUE
  DO 149 I=1,MR
     II = IP+I
     Q(I,J) = Q(I,JM2)+P(II)+B(I)
  149 CONTINUE
  150 CONTINUE
  LR = KR
  KR = KR+JR
  151 CONTINUE
  go to (152,153),MIXBND
  152 NR = (NLAST-1)/JR+1
  if (NR  <=  3) go to 155
  go to 154
  153 NR = NLAST/JR
  if (NR  <=  1) go to 192
  154 I2R = JR
  NRODPR = NROD
  go to 104
  155 CONTINUE
!
!      BEGIN SOLUTION
!
  J = 1+JR
  JM1 = J-I2R
  JP1 = J+I2R
  JM2 = NLAST-I2R
  if (NR  ==  2) go to 184
  if (LR  /=  0) go to 170
  if (N  /=  3) go to 161
!
!     CASE N = 3.
!
  go to (156,168),ISTAG
  156 CONTINUE
  DO 157 I=1,MR
     B(I) = Q(I,2)
  157 CONTINUE
  TCOS(1) = 0.
  call TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  DO 158 I=1,MR
     Q(I,2) = B(I)
     B(I) = 4.*B(I)+Q(I,1)+2.*Q(I,3)
  158 CONTINUE
  TCOS(1) = -2.
  TCOS(2) = 2.
  I1 = 2
  I2 = 0
  call TRIX (I1,I2,MR,A,BB,C,B,TCOS,D,W)
  DO 159 I=1,MR
     Q(I,2) = Q(I,2)+B(I)
     B(I) = Q(I,1)+2.*Q(I,2)
  159 CONTINUE
  TCOS(1) = 0.
  call TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  DO 160 I=1,MR
     Q(I,1) = B(I)
  160 CONTINUE
  JR = 1
  I2R = 0
  go to 194
!
!     CASE N = 2**P+1
!
  161 CONTINUE
  go to (162,170),ISTAG
  162 CONTINUE
  DO 163 I=1,MR
     B(I) = Q(I,J)+.5*Q(I,1)-Q(I,JM1)+Q(I,NLAST)-Q(I,JM2)
  163 CONTINUE
  call COSGEN (JR,1,0.5,0.0,TCOS)
  call TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
  DO 164 I=1,MR
     Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
     B(I) = Q(I,1)+2.*Q(I,NLAST)+4.*Q(I,J)
  164 CONTINUE
  JR2 = 2*JR
  call COSGEN (JR,1,0.0,0.0,TCOS)
  DO 165 I=1,JR
     I1 = JR+I
     I2 = JR+1-I
     TCOS(I1) = -TCOS(I2)
  165 CONTINUE
  call TRIX (JR2,0,MR,A,BB,C,B,TCOS,D,W)
  DO 166 I=1,MR
     Q(I,J) = Q(I,J)+B(I)
     B(I) = Q(I,1)+2.*Q(I,J)
  166 CONTINUE
  call COSGEN (JR,1,0.5,0.0,TCOS)
  call TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
  DO 167 I=1,MR
     Q(I,1) = .5*Q(I,1)-Q(I,JM1)+B(I)
  167 CONTINUE
  go to 194
!
!     CASE OF GENERAL N WITH NR = 3 .
!
  168 DO 169 I=1,MR
     B(I) = Q(I,2)
     Q(I,2) = 0.
     B2(I) = Q(I,3)
     B3(I) = Q(I,1)
  169 CONTINUE
  JR = 1
  I2R = 0
  J = 2
  go to 177
  170 CONTINUE
  DO 171 I=1,MR
     B(I) = .5*Q(I,1)-Q(I,JM1)+Q(I,J)
  171 CONTINUE
  if (NROD  /=  0) go to 173
  DO 172 I=1,MR
     II = IP+I
     B(I) = B(I)+P(II)
  172 CONTINUE
  go to 175
  173 DO 174 I=1,MR
     B(I) = B(I)+Q(I,NLAST)-Q(I,JM2)
  174 CONTINUE
  175 CONTINUE
  DO 176 I=1,MR
     T = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
     Q(I,J) = T
     B2(I) = Q(I,NLAST)+T
     B3(I) = Q(I,1)+2.*T
  176 CONTINUE
  177 CONTINUE
  K1 = KR+2*JR-1
  K2 = KR+JR
  TCOS(K1+1) = -2.
  K4 = K1+3-ISTAG
  call COSGEN (K2+ISTAG-2,1,0.0,FNUM,TCOS(K4))
  K4 = K1+K2+1
  call COSGEN (JR-1,1,0.0,1.0,TCOS(K4))
  call S1MERG (TCOS,K1,K2,K1+K2,JR-1,0)
  K3 = K1+K2+LR
  call COSGEN (JR,1,0.5,0.0,TCOS(K3+1))
  K4 = K3+JR+1
  call COSGEN (KR,1,0.5,FDEN,TCOS(K4))
  call S1MERG (TCOS,K3,JR,K3+JR,KR,K1)
  if (LR  ==  0) go to 178
  call COSGEN (LR,1,0.5,FDEN,TCOS(K4))
  call S1MERG (TCOS,K3,JR,K3+JR,LR,K3-LR)
  call COSGEN (KR,1,0.5,FDEN,TCOS(K4))
  178 K3 = KR
  K4 = KR
  call TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
  DO 179 I=1,MR
     B(I) = B(I)+B2(I)+B3(I)
  179 CONTINUE
  TCOS(1) = 2.
  call TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  DO 180 I=1,MR
     Q(I,J) = Q(I,J)+B(I)
     B(I) = Q(I,1)+2.*Q(I,J)
  180 CONTINUE
  call COSGEN (JR,1,0.5,0.0,TCOS)
  call TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
  if (JR  /=  1) go to 182
  DO 181 I=1,MR
     Q(I,1) = B(I)
  181 CONTINUE
  go to 194
  182 CONTINUE
  DO 183 I=1,MR
     Q(I,1) = .5*Q(I,1)-Q(I,JM1)+B(I)
  183 CONTINUE
  go to 194
  184 CONTINUE
  if (N  /=  2) go to 188
!
!     CASE  N = 2
!
  DO 185 I=1,MR
     B(I) = Q(I,1)
  185 CONTINUE
  TCOS(1) = 0.
  call TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  DO 186 I=1,MR
     Q(I,1) = B(I)
     B(I) = 2.*(Q(I,2)+B(I))*FISTAG
  186 CONTINUE
  TCOS(1) = -FISTAG
  TCOS(2) = 2.
  call TRIX (2,0,MR,A,BB,C,B,TCOS,D,W)
  DO 187 I=1,MR
     Q(I,1) = Q(I,1)+B(I)
  187 CONTINUE
  JR = 1
  I2R = 0
  go to 194
  188 CONTINUE
!
!     CASE OF GENERAL N AND NR = 2 .
!
  DO 189 I=1,MR
     II = IP+I
     B3(I) = 0.
     B(I) = Q(I,1)+2.*P(II)
     Q(I,1) = .5*Q(I,1)-Q(I,JM1)
     B2(I) = 2.*(Q(I,1)+Q(I,NLAST))
  189 CONTINUE
  K1 = KR+JR-1
  TCOS(K1+1) = -2.
  K4 = K1+3-ISTAG
  call COSGEN (KR+ISTAG-2,1,0.0,FNUM,TCOS(K4))
  K4 = K1+KR+1
  call COSGEN (JR-1,1,0.0,1.0,TCOS(K4))
  call S1MERG (TCOS,K1,KR,K1+KR,JR-1,0)
  call COSGEN (KR,1,0.5,FDEN,TCOS(K1+1))
  K2 = KR
  K4 = K1+K2+1
  call COSGEN (LR,1,0.5,FDEN,TCOS(K4))
  K3 = LR
  K4 = 0
  call TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
  DO 190 I=1,MR
     B(I) = B(I)+B2(I)
  190 CONTINUE
  TCOS(1) = 2.
  call TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  DO 191 I=1,MR
     Q(I,1) = Q(I,1)+B(I)
  191 CONTINUE
  go to 194
  192 DO 193 I=1,MR
     B(I) = Q(I,NLAST)
  193 CONTINUE
  go to 196
  194 CONTINUE
!
!     START BACK SUBSTITUTION.
!
  J = NLAST-JR
  DO 195 I=1,MR
     B(I) = Q(I,NLAST)+Q(I,J)
  195 CONTINUE
  196 JM2 = NLAST-I2R
  if (JR  /=  1) go to 198
  DO 197 I=1,MR
     Q(I,NLAST) = 0.
  197 CONTINUE
  go to 202
  198 CONTINUE
  if (NROD  /=  0) go to 200
  DO 199 I=1,MR
     II = IP+I
     Q(I,NLAST) = P(II)
  199 CONTINUE
  IP = IP-MR
  go to 202
  200 DO 201 I=1,MR
     Q(I,NLAST) = Q(I,NLAST)-Q(I,JM2)
  201 CONTINUE
  202 CONTINUE
  call COSGEN (KR,1,0.5,FDEN,TCOS)
  call COSGEN (LR,1,0.5,FDEN,TCOS(KR+1))
  if (LR  /=  0) go to 204
  DO 203 I=1,MR
     B(I) = FISTAG*B(I)
  203 CONTINUE
  204 CONTINUE
  call TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
  DO 205 I=1,MR
     Q(I,NLAST) = Q(I,NLAST)+B(I)
  205 CONTINUE
  NLASTP = NLAST
  206 CONTINUE
  JSTEP = JR
  JR = I2R
  I2R = I2R/2
  if (JR  ==  0) go to 222
  go to (207,208),MIXBND
  207 JSTART = 1+JR
  go to 209
  208 JSTART = JR
  209 CONTINUE
  KR = KR-JR
  if (NLAST+JR  >  N) go to 210
  KR = KR-JR
  NLAST = NLAST+JR
  JSTOP = NLAST-JSTEP
  go to 211
  210 CONTINUE
  JSTOP = NLAST-JR
  211 CONTINUE
  LR = KR-JR
  call COSGEN (JR,1,0.5,0.0,TCOS)
  DO 221 J=JSTART,JSTOP,JSTEP
     JM2 = J-JR
     JP2 = J+JR
     if (J  /=  JR) go to 213
     DO 212 I=1,MR
        B(I) = Q(I,J)+Q(I,JP2)
  212    CONTINUE
     go to 215
  213    CONTINUE
     DO 214 I=1,MR
        B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  214    CONTINUE
  215    CONTINUE
     if (JR  /=  1) go to 217
     DO 216 I=1,MR
        Q(I,J) = 0.
  216    CONTINUE
     go to 219
  217    CONTINUE
     JM1 = J-I2R
     JP1 = J+I2R
     DO 218 I=1,MR
        Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  218    CONTINUE
  219    CONTINUE
     call TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
     DO 220 I=1,MR
        Q(I,J) = Q(I,J)+B(I)
  220    CONTINUE
  221 CONTINUE
  NROD = 1
  if (NLAST+I2R  <=  N) NROD = 0
  if (NLASTP  /=  NLAST) go to 194
  go to 206
  222 CONTINUE
!
!     return STORAGE REQUIREMENTS FOR P VECTORS.
!
  W(1) = IPSTOR
  return
end
