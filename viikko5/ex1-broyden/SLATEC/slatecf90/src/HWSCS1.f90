subroutine HWSCS1 (INTL, TS, TF, M, MBDCND, BDTS, BDTF, RS, RF, N, &
     NBDCND, BDRS, BDRF, ELMBDA, F, IDIMF, PERTRB, W, S, AN, BN, CN, &
     R, AM, BM, CM, SINT, BMH)
!
!! HWSCS1 is subsidiary to HWSCSP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (HWSCS1-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  HWSCSP
!***ROUTINES CALLED  BLKTRI
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891009  Removed unreferenced variables.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  HWSCS1
  DIMENSION       F(IDIMF,*) ,BDRS(*)    ,BDRF(*)    ,BDTS(*)    , &
                  BDTF(*)    ,AM(*)      ,BM(*)      ,CM(*)      , &
                  AN(*)      ,BN(*)      ,CN(*)      ,S(*)       , &
                  R(*)       ,SINT(*)    ,BMH(*)     ,W(*)
!***FIRST EXECUTABLE STATEMENT  HWSCS1
  MP1 = M+1
  DTH = (TF-TS)/M
  TDT = DTH+DTH
  HDTH = DTH/2.
  SDTS = 1./(DTH*DTH)
  DO 102 I=1,MP1
     THETA = TS+(I-1)*DTH
     SINT(I) = SIN(THETA)
     if (SINT(I)) 101,102,101
  101    T1 = SDTS/SINT(I)
     AM(I) = T1*SIN(THETA-HDTH)
     CM(I) = T1*SIN(THETA+HDTH)
     BM(I) = -(AM(I)+CM(I))
  102 CONTINUE
  NP1 = N+1
  DR = (RF-RS)/N
  HDR = DR/2.
  TDR = DR+DR
  DR2 = DR*DR
  CZR = 6.*DTH/(DR2*(COS(TS)-COS(TF)))
  DO 103 J=1,NP1
     R(J) = RS+(J-1)*DR
     AN(J) = (R(J)-HDR)**2/DR2
     CN(J) = (R(J)+HDR)**2/DR2
     BN(J) = -(AN(J)+CN(J))
  103 CONTINUE
  MP = 1
  NP = 1
!
! BOUNDARY CONDITION AT PHI=PS
!
  go to (104,104,105,105,106,106,104,105,106),MBDCND
  104 AT = AM(2)
  ITS = 2
  go to 107
  105 AT = AM(1)
  ITS = 1
  CM(1) = CM(1)+AM(1)
  go to 107
  106 ITS = 1
  BM(1) = -4.*SDTS
  CM(1) = -BM(1)
!
! BOUNDARY CONDITION AT PHI=PF
!
  107 go to (108,109,109,108,108,109,110,110,110),MBDCND
  108 CT = CM(M)
  ITF = M
  go to 111
  109 CT = CM(M+1)
  AM(M+1) = AM(M+1)+CM(M+1)
  ITF = M+1
  go to 111
  110 ITF = M+1
  AM(M+1) = 4.*SDTS
  BM(M+1) = -AM(M+1)
  111 WTS = SINT(ITS+1)*AM(ITS+1)/CM(ITS)
  WTF = SINT(ITF-1)*CM(ITF-1)/AM(ITF)
  ITSP = ITS+1
  ITFM = ITF-1
!
! BOUNDARY CONDITION AT R=RS
!
  ICTR = 0
  go to (112,112,113,113,114,114),NBDCND
  112 AR = AN(2)
  JRS = 2
  go to 118
  113 AR = AN(1)
  JRS = 1
  CN(1) = CN(1)+AN(1)
  go to 118
  114 JRS = 2
  ICTR = 1
  S(N) = AN(N)/BN(N)
  DO 115 J=3,N
     L = N-J+2
     S(L) = AN(L)/(BN(L)-CN(L)*S(L+1))
  115 CONTINUE
  S(2) = -S(2)
  DO 116 J=3,N
     S(J) = -S(J)*S(J-1)
  116 CONTINUE
  WTNM = WTS+WTF
  DO 117 I=ITSP,ITFM
     WTNM = WTNM+SINT(I)
  117 CONTINUE
  YPS = CZR*WTNM*(S(2)-1.)
!
! BOUNDARY CONDITION AT R=RF
!
  118 go to (119,120,120,119,119,120),NBDCND
  119 CR = CN(N)
  JRF = N
  go to 121
  120 CR = CN(N+1)
  AN(N+1) = AN(N+1)+CN(N+1)
  JRF = N+1
  121 WRS = AN(JRS+1)*R(JRS)**2/CN(JRS)
  WRF = CN(JRF-1)*R(JRF)**2/AN(JRF)
  WRZ = AN(JRS)/CZR
  JRSP = JRS+1
  JRFM = JRF-1
  MUNK = ITF-ITS+1
  NUNK = JRF-JRS+1
  DO 122 I=ITS,ITF
     BMH(I) = BM(I)
  122 CONTINUE
  ISING = 0
  go to (132,132,123,132,132,123),NBDCND
  123 go to (132,132,124,132,132,124,132,124,124),MBDCND
  124 if (ELMBDA) 132,125,125
  125 ISING = 1
  SUM = WTS*WRS+WTS*WRF+WTF*WRS+WTF*WRF
  if (ICTR) 126,127,126
  126 SUM = SUM+WRZ
  127 DO 129 J=JRSP,JRFM
     R2 = R(J)**2
     DO 128 I=ITSP,ITFM
        SUM = SUM+R2*SINT(I)
  128    CONTINUE
  129 CONTINUE
  DO 130 J=JRSP,JRFM
     SUM = SUM+(WTS+WTF)*R(J)**2
  130 CONTINUE
  DO 131 I=ITSP,ITFM
     SUM = SUM+(WRS+WRF)*SINT(I)
  131 CONTINUE
  HNE = SUM
  132 go to (133,133,133,133,134,134,133,133,134),MBDCND
  133 BM(ITS) = BMH(ITS)+ELMBDA/SINT(ITS)**2
  134 go to (135,135,135,135,135,135,136,136,136),MBDCND
  135 BM(ITF) = BMH(ITF)+ELMBDA/SINT(ITF)**2
  136 DO 137 I=ITSP,ITFM
     BM(I) = BMH(I)+ELMBDA/SINT(I)**2
  137 CONTINUE
  go to (138,138,140,140,142,142,138,140,142),MBDCND
  138 DO 139 J=JRS,JRF
     F(2,J) = F(2,J)-AT*F(1,J)/R(J)**2
  139 CONTINUE
  go to 142
  140 DO 141 J=JRS,JRF
     F(1,J) = F(1,J)+TDT*BDTS(J)*AT/R(J)**2
  141 CONTINUE
  142 go to (143,145,145,143,143,145,147,147,147),MBDCND
  143 DO 144 J=JRS,JRF
     F(M,J) = F(M,J)-CT*F(M+1,J)/R(J)**2
  144 CONTINUE
  go to 147
  145 DO 146 J=JRS,JRF
     F(M+1,J) = F(M+1,J)-TDT*BDTF(J)*CT/R(J)**2
  146 CONTINUE
  147 go to (151,151,153,153,148,148),NBDCND
  148 if (MBDCND-3) 155,149,155
  149 YHLD = F(ITS,1)-CZR/TDT*(SIN(TF)*BDTF(2)-SIN(TS)*BDTS(2))
  DO 150 I=1,MP1
     F(I,1) = YHLD
  150 CONTINUE
  go to 155
  151 RS2 = (RS+DR)**2
  DO 152 I=ITS,ITF
     F(I,2) = F(I,2)-AR*F(I,1)/RS2
  152 CONTINUE
  go to 155
  153 DO 154 I=ITS,ITF
     F(I,1) = F(I,1)+TDR*BDRS(I)*AR/RS**2
  154 CONTINUE
  155 go to (156,158,158,156,156,158),NBDCND
  156 RF2 = (RF-DR)**2
  DO 157 I=ITS,ITF
     F(I,N) = F(I,N)-CR*F(I,N+1)/RF2
  157 CONTINUE
  go to 160
  158 DO 159 I=ITS,ITF
     F(I,N+1) = F(I,N+1)-TDR*BDRF(I)*CR/RF**2
  159 CONTINUE
  160 CONTINUE
  PERTRB = 0.
  if (ISING) 161,170,161
  161 SUM = WTS*WRS*F(ITS,JRS)+WTS*WRF*F(ITS,JRF)+WTF*WRS*F(ITF,JRS)+ &
        WTF*WRF*F(ITF,JRF)
  if (ICTR) 162,163,162
  162 SUM = SUM+WRZ*F(ITS,1)
  163 DO 165 J=JRSP,JRFM
     R2 = R(J)**2
     DO 164 I=ITSP,ITFM
        SUM = SUM+R2*SINT(I)*F(I,J)
  164    CONTINUE
  165 CONTINUE
  DO 166 J=JRSP,JRFM
     SUM = SUM+R(J)**2*(WTS*F(ITS,J)+WTF*F(ITF,J))
  166 CONTINUE
  DO 167 I=ITSP,ITFM
     SUM = SUM+SINT(I)*(WRS*F(I,JRS)+WRF*F(I,JRF))
  167 CONTINUE
  PERTRB = SUM/HNE
  DO 169 J=1,NP1
     DO 168 I=1,MP1
        F(I,J) = F(I,J)-PERTRB
  168    CONTINUE
  169 CONTINUE
  170 DO 172 J=JRS,JRF
     RSQ = R(J)**2
     DO 171 I=ITS,ITF
        F(I,J) = RSQ*F(I,J)
  171    CONTINUE
  172 CONTINUE
  IFLG = INTL
  173 call BLKTRI (IFLG,NP,NUNK,AN(JRS),BN(JRS),CN(JRS),MP,MUNK, &
               AM(ITS),BM(ITS),CM(ITS),IDIMF,F(ITS,JRS),IERROR,W)
  IFLG = IFLG+1
  if (IFLG-1) 174,173,174
  174 if (NBDCND) 177,175,177
  175 DO 176 I=1,MP1
     F(I,JRF+1) = F(I,JRS)
  176 CONTINUE
  177 if (MBDCND) 180,178,180
  178 DO 179 J=1,NP1
     F(ITF+1,J) = F(ITS,J)
  179 CONTINUE
  180 XP = 0.
  if (ICTR) 181,188,181
  181 if (ISING) 186,182,186
  182 SUM = WTS*F(ITS,2)+WTF*F(ITF,2)
  DO 183 I=ITSP,ITFM
     SUM = SUM+SINT(I)*F(I,2)
  183 CONTINUE
  YPH = CZR*SUM
  XP = (F(ITS,1)-YPH)/YPS
  DO 185 J=JRS,JRF
     XPS = XP*S(J)
     DO 184 I=ITS,ITF
        F(I,J) = F(I,J)+XPS
  184    CONTINUE
  185 CONTINUE
  186 DO 187 I=1,MP1
     F(I,1) = XP
  187 CONTINUE
  188 RETURN
end
