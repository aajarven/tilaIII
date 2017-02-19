subroutine POISD2 (MR, NR, ISTAG, BA, BB, BC, Q, IDIMQ, B, W, D, &
     TCOS, P)
!
!! POISD2 is subsidiary to GENBUN.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (POISD2-S, CMPOSD-C)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Subroutine to solve Poisson's equation for Dirichlet boundary
!     conditions.
!
!     ISTAG = 1 if the last diagonal block is the matrix A.
!     ISTAG = 2 if the last diagonal block is the matrix A+I.
!
!***SEE ALSO  GENBUN
!***ROUTINES CALLED  COSGEN, S1MERG, TRIX
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   920130  Modified to use merge routine S1MERG rather than deleted
!           routine MERGE.  (WRB)
!***END PROLOGUE  POISD2
!
  DIMENSION       Q(IDIMQ,*) ,BA(*)      ,BB(*)      ,BC(*)      , &
                  TCOS(*)    ,B(*)       ,D(*)       ,W(*)       , &
                  P(*)
!***FIRST EXECUTABLE STATEMENT  POISD2
  M = MR
  N = NR
  JSH = 0
  FI = 1./ISTAG
  IP = -M
  IPSTOR = 0
  go to (101,102),ISTAG
  101 KR = 0
  IRREG = 1
  if (N  >  1) go to 106
  TCOS(1) = 0.
  go to 103
  102 KR = 1
  JSTSAV = 1
  IRREG = 2
  if (N  >  1) go to 106
  TCOS(1) = -1.
  103 DO 104 I=1,M
     B(I) = Q(I,1)
  104 CONTINUE
  call TRIX (1,0,M,BA,BB,BC,B,TCOS,D,W)
  DO 105 I=1,M
     Q(I,1) = B(I)
  105 CONTINUE
  go to 183
  106 LR = 0
  DO 107 I=1,M
     P(I) = 0.
  107 CONTINUE
  NUN = N
  JST = 1
  JSP = N
!
!     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2.
!
  108 L = 2*JST
  NODD = 2-2*((NUN+1)/2)+NUN
!
!     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2.
!
  go to (110,109),NODD
  109 JSP = JSP-L
  go to 111
  110 JSP = JSP-JST
  if (IRREG  /=  1) JSP = JSP-L
  111 CONTINUE
!
!     REGULAR REDUCTION
!
  call COSGEN (JST,1,0.5,0.0,TCOS)
  if (L  >  JSP) go to 118
  DO 117 J=L,JSP,L
     JM1 = J-JSH
     JP1 = J+JSH
     JM2 = J-JST
     JP2 = J+JST
     JM3 = JM2-JSH
     JP3 = JP2+JSH
     if (JST  /=  1) go to 113
     DO 112 I=1,M
        B(I) = 2.*Q(I,J)
        Q(I,J) = Q(I,JM2)+Q(I,JP2)
  112    CONTINUE
     go to 115
  113    DO 114 I=1,M
        T = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
        B(I) = T+Q(I,J)-Q(I,JM3)-Q(I,JP3)
        Q(I,J) = T
  114    CONTINUE
  115    CONTINUE
     call TRIX (JST,0,M,BA,BB,BC,B,TCOS,D,W)
     DO 116 I=1,M
        Q(I,J) = Q(I,J)+B(I)
  116    CONTINUE
  117 CONTINUE
!
!     REDUCTION FOR LAST UNKNOWN
!
  118 go to (119,136),NODD
  119 go to (152,120),IRREG
!
!     ODD NUMBER OF UNKNOWNS
!
  120 JSP = JSP+L
  J = JSP
  JM1 = J-JSH
  JP1 = J+JSH
  JM2 = J-JST
  JP2 = J+JST
  JM3 = JM2-JSH
  go to (123,121),ISTAG
  121 CONTINUE
  if (JST  /=  1) go to 123
  DO 122 I=1,M
     B(I) = Q(I,J)
     Q(I,J) = 0.
  122 CONTINUE
  go to 130
  123 go to (124,126),NODDPR
  124 DO 125 I=1,M
     IP1 = IP+I
     B(I) = .5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))+P(IP1)+Q(I,J)
  125 CONTINUE
  go to 128
  126 DO 127 I=1,M
     B(I) = .5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))+Q(I,JP2)-Q(I,JP1)+Q(I,J)
  127 CONTINUE
  128 DO 129 I=1,M
     Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  129 CONTINUE
  130 call TRIX (JST,0,M,BA,BB,BC,B,TCOS,D,W)
  IP = IP+M
  IPSTOR = MAX(IPSTOR,IP+M)
  DO 131 I=1,M
     IP1 = IP+I
     P(IP1) = Q(I,J)+B(I)
     B(I) = Q(I,JP2)+P(IP1)
  131 CONTINUE
  if (LR  /=  0) go to 133
  DO 132 I=1,JST
     KRPI = KR+I
     TCOS(KRPI) = TCOS(I)
  132 CONTINUE
  go to 134
  133 CONTINUE
  call COSGEN (LR,JSTSAV,0.,FI,TCOS(JST+1))
  call S1MERG (TCOS,0,JST,JST,LR,KR)
  134 CONTINUE
  call COSGEN (KR,JSTSAV,0.0,FI,TCOS)
  call TRIX (KR,KR,M,BA,BB,BC,B,TCOS,D,W)
  DO 135 I=1,M
     IP1 = IP+I
     Q(I,J) = Q(I,JM2)+B(I)+P(IP1)
  135 CONTINUE
  LR = KR
  KR = KR+L
  go to 152
!
!     EVEN NUMBER OF UNKNOWNS
!
  136 JSP = JSP+L
  J = JSP
  JM1 = J-JSH
  JP1 = J+JSH
  JM2 = J-JST
  JP2 = J+JST
  JM3 = JM2-JSH
  go to (137,138),IRREG
  137 CONTINUE
  JSTSAV = JST
  IDEG = JST
  KR = L
  go to 139
  138 call COSGEN (KR,JSTSAV,0.0,FI,TCOS)
  call COSGEN (LR,JSTSAV,0.0,FI,TCOS(KR+1))
  IDEG = KR
  KR = KR+JST
  139 if (JST  /=  1) go to 141
  IRREG = 2
  DO 140 I=1,M
     B(I) = Q(I,J)
     Q(I,J) = Q(I,JM2)
  140 CONTINUE
  go to 150
  141 DO 142 I=1,M
     B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  142 CONTINUE
  go to (143,145),IRREG
  143 DO 144 I=1,M
     Q(I,J) = Q(I,JM2)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  144 CONTINUE
  IRREG = 2
  go to 150
  145 CONTINUE
  go to (146,148),NODDPR
  146 DO 147 I=1,M
     IP1 = IP+I
     Q(I,J) = Q(I,JM2)+P(IP1)
  147 CONTINUE
  IP = IP-M
  go to 150
  148 DO 149 I=1,M
     Q(I,J) = Q(I,JM2)+Q(I,J)-Q(I,JM1)
  149 CONTINUE
  150 call TRIX (IDEG,LR,M,BA,BB,BC,B,TCOS,D,W)
  DO 151 I=1,M
     Q(I,J) = Q(I,J)+B(I)
  151 CONTINUE
  152 NUN = NUN/2
  NODDPR = NODD
  JSH = JST
  JST = 2*JST
  if (NUN  >=  2) go to 108
!
!     START SOLUTION.
!
  J = JSP
  DO 153 I=1,M
     B(I) = Q(I,J)
  153 CONTINUE
  go to (154,155),IRREG
  154 CONTINUE
  call COSGEN (JST,1,0.5,0.0,TCOS)
  IDEG = JST
  go to 156
  155 KR = LR+JST
  call COSGEN (KR,JSTSAV,0.0,FI,TCOS)
  call COSGEN (LR,JSTSAV,0.0,FI,TCOS(KR+1))
  IDEG = KR
  156 CONTINUE
  call TRIX (IDEG,LR,M,BA,BB,BC,B,TCOS,D,W)
  JM1 = J-JSH
  JP1 = J+JSH
  go to (157,159),IRREG
  157 DO 158 I=1,M
     Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
  158 CONTINUE
  go to 164
  159 go to (160,162),NODDPR
  160 DO 161 I=1,M
     IP1 = IP+I
     Q(I,J) = P(IP1)+B(I)
  161 CONTINUE
  IP = IP-M
  go to 164
  162 DO 163 I=1,M
     Q(I,J) = Q(I,J)-Q(I,JM1)+B(I)
  163 CONTINUE
  164 CONTINUE
!
!     START BACK SUBSTITUTION.
!
  JST = JST/2
  JSH = JST/2
  NUN = 2*NUN
  if (NUN  >  N) go to 183
  DO 182 J=JST,N,L
     JM1 = J-JSH
     JP1 = J+JSH
     JM2 = J-JST
     JP2 = J+JST
     if (J  >  JST) go to 166
     DO 165 I=1,M
        B(I) = Q(I,J)+Q(I,JP2)
  165    CONTINUE
     go to 170
  166    if (JP2  <=  N) go to 168
     DO 167 I=1,M
        B(I) = Q(I,J)+Q(I,JM2)
  167    CONTINUE
     if (JST  <  JSTSAV) IRREG = 1
     go to (170,171),IRREG
  168    DO 169 I=1,M
        B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  169    CONTINUE
  170    CONTINUE
     call COSGEN (JST,1,0.5,0.0,TCOS)
     IDEG = JST
     JDEG = 0
     go to 172
  171    if (J+L  >  N) LR = LR-JST
     KR = JST+LR
     call COSGEN (KR,JSTSAV,0.0,FI,TCOS)
     call COSGEN (LR,JSTSAV,0.0,FI,TCOS(KR+1))
     IDEG = KR
     JDEG = LR
  172    CONTINUE
     call TRIX (IDEG,JDEG,M,BA,BB,BC,B,TCOS,D,W)
     if (JST  >  1) go to 174
     DO 173 I=1,M
        Q(I,J) = B(I)
  173    CONTINUE
     go to 182
  174    if (JP2  >  N) go to 177
  175    DO 176 I=1,M
        Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
  176    CONTINUE
     go to 182
  177    go to (175,178),IRREG
  178    if (J+JSH  >  N) go to 180
     DO 179 I=1,M
        IP1 = IP+I
        Q(I,J) = B(I)+P(IP1)
  179    CONTINUE
     IP = IP-M
     go to 182
  180    DO 181 I=1,M
        Q(I,J) = B(I)+Q(I,J)-Q(I,JM1)
  181    CONTINUE
  182 CONTINUE
  L = L/2
  go to 164
  183 CONTINUE
!
!     return STORAGE REQUIREMENTS FOR P VECTORS.
!
  W(1) = IPSTOR
  return
end
