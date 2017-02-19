subroutine CBLKT1 (N, AN, BN, CN, M, AM, BM, CM, IDIMY, Y, B, W1, &
     W2, W3, WD, WW, WU, PRDCT, CPRDCT)
!
!! CBLKT1 is subsidiary to CBLKTR.
!
!***LIBRARY   SLATEC
!***TYPE      COMPLEX (BLKTR1-S, CBLKT1-C)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
! CBLKT1 solves the linear system of routine CBLKTR.
!
! B  contains the roots of all the B polynomials.
! W1,W2,W3,WD,WW,WU  are all working arrays.
! PRDCT is either PROCP or PROC depending on whether the boundary
! conditions in the M direction are periodic or not.
! CPRDCT is either CPROCP or CPROC which are called if some of the zeros
! of the B polynomials are complex.
!
!***SEE ALSO  CBLKTR
!***ROUTINES CALLED  INXCA, INXCB, INXCC
!***COMMON BLOCKS    CCBLK
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  CBLKT1
!
  DIMENSION       AN(*)      ,BN(*)      ,CN(*)      ,AM(*)      , &
                  BM(*)      ,CM(*)      ,B(*)       ,W1(*)      , &
                  W2(*)      ,W3(*)      ,WD(*)      ,WW(*)      , &
                  WU(*)      ,Y(IDIMY,*)
  COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        , &
                  NM         ,NCMPLX     ,IK
  COMPLEX         AM         ,BM         ,CM         ,Y          , &
                  W1         ,W2         ,W3         ,WD         , &
                  WW         ,WU
!***FIRST EXECUTABLE STATEMENT  CBLKT1
  KDO = K-1
  DO 109 L=1,KDO
     IR = L-1
     I2 = 2**IR
     I1 = I2/2
     I3 = I2+I1
     I4 = I2+I2
     IRM1 = IR-1
     call INXCB (I2,IR,IM2,NM2)
     call INXCB (I1,IRM1,IM3,NM3)
     call INXCB (I3,IRM1,IM1,NM1)
     call PRDCT (NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,Y(1,I2),W3, &
                 M,AM,BM,CM,WD,WW,WU)
     if = 2**K
     DO 108 I=I4,IF,I4
        if (I-NM) 101,101,108
  101       IPI1 = I+I1
        IPI2 = I+I2
        IPI3 = I+I3
        call INXCC (I,IR,IDXC,NC)
        if (I-IF) 102,108,108
  102       call INXCA (I,IR,IDXA,NA)
        call INXCB (I-I1,IRM1,IM1,NM1)
        call INXCB (IPI2,IR,IP2,NP2)
        call INXCB (IPI1,IRM1,IP1,NP1)
        call INXCB (IPI3,IRM1,IP3,NP3)
        call PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W3,W1,M,AM, &
                    BM,CM,WD,WW,WU)
        if (IPI2-NM) 105,105,103
  103       DO 104 J=1,M
           W3(J) = (0.,0.)
           W2(J) = (0.,0.)
  104       CONTINUE
        go to 106
  105       call PRDCT (NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM, &
                    Y(1,IPI2),W3,M,AM,BM,CM,WD,WW,WU)
        call PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W3,W2,M,AM, &
                    BM,CM,WD,WW,WU)
  106       DO 107 J=1,M
           Y(J,I) = W1(J)+W2(J)+Y(J,I)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
  if (NPP) 132,110,132
!
!     THE PERIODIC CASE IS TREATED USING THE CAPACITANCE MATRIX METHOD
!
  110 if = 2**K
  I = IF/2
  I1 = I/2
  call INXCB (I-I1,K-2,IM1,NM1)
  call INXCB (I+I1,K-2,IP1,NP1)
  call INXCB (I,K-1,IZ,NZ)
  call PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,Y(1,I),W1,M,AM, &
              BM,CM,WD,WW,WU)
  IZR = I
  DO 111 J=1,M
     W2(J) = W1(J)
  111 CONTINUE
  DO 113 LL=2,K
     L = K-LL+1
     IR = L-1
     I2 = 2**IR
     I1 = I2/2
     I = I2
     call INXCC (I,IR,IDXC,NC)
     call INXCB (I,IR,IZ,NZ)
     call INXCB (I-I1,IR-1,IM1,NM1)
     call INXCB (I+I1,IR-1,IP1,NP1)
     call PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W1,W1,M,AM,BM, &
                 CM,WD,WW,WU)
     DO 112 J=1,M
        W1(J) = Y(J,I)+W1(J)
  112    CONTINUE
     call PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,W1,M,AM, &
                 BM,CM,WD,WW,WU)
  113 CONTINUE
  DO 118 LL=2,K
     L = K-LL+1
     IR = L-1
     I2 = 2**IR
     I1 = I2/2
     I4 = I2+I2
     IFD = IF-I2
     DO 117 I=I2,IFD,I4
        if (I-I2-IZR) 117,114,117
  114       if (I-NM) 115,115,118
  115       call INXCA (I,IR,IDXA,NA)
        call INXCB (I,IR,IZ,NZ)
        call INXCB (I-I1,IR-1,IM1,NM1)
        call INXCB (I+I1,IR-1,IP1,NP1)
        call PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W2,W2,M,AM, &
                    BM,CM,WD,WW,WU)
        DO 116 J=1,M
           W2(J) = Y(J,I)+W2(J)
  116       CONTINUE
        call PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W2,W2,M, &
                    AM,BM,CM,WD,WW,WU)
        IZR = I
        if (I-NM) 117,119,117
  117    CONTINUE
  118 CONTINUE
  119 DO 120 J=1,M
     Y(J,NM+1) = Y(J,NM+1)-CN(NM+1)*W1(J)-AN(NM+1)*W2(J)
  120 CONTINUE
  call INXCB (IF/2,K-1,IM1,NM1)
  call INXCB (IF,K-1,IP,NP)
  if (NCMPLX) 121,122,121
  121 call CPRDCT (NM+1,B(IP),NM1,B(IM1),0,DUM,0,DUM,Y(1,NM+1), &
               Y(1,NM+1),M,AM,BM,CM,W1,W3,WW)
  go to 123
  122 call PRDCT (NM+1,B(IP),NM1,B(IM1),0,DUM,0,DUM,Y(1,NM+1), &
              Y(1,NM+1),M,AM,BM,CM,WD,WW,WU)
  123 DO 124 J=1,M
     W1(J) = AN(1)*Y(J,NM+1)
     W2(J) = CN(NM)*Y(J,NM+1)
     Y(J,1) = Y(J,1)-W1(J)
     Y(J,NM) = Y(J,NM)-W2(J)
  124 CONTINUE
  DO 126 L=1,KDO
     IR = L-1
     I2 = 2**IR
     I4 = I2+I2
     I1 = I2/2
     I = I4
     call INXCA (I,IR,IDXA,NA)
     call INXCB (I-I2,IR,IM2,NM2)
     call INXCB (I-I2-I1,IR-1,IM3,NM3)
     call INXCB (I-I1,IR-1,IM1,NM1)
     call PRDCT (NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,W1,W1,M,AM, &
                 BM,CM,WD,WW,WU)
     call PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W1,W1,M,AM,BM, &
                 CM,WD,WW,WU)
     DO 125 J=1,M
        Y(J,I) = Y(J,I)-W1(J)
  125    CONTINUE
  126 CONTINUE
!
  IZR = NM
  DO 131 L=1,KDO
     IR = L-1
     I2 = 2**IR
     I1 = I2/2
     I3 = I2+I1
     I4 = I2+I2
     IRM1 = IR-1
     DO 130 I=I4,IF,I4
        IPI1 = I+I1
        IPI2 = I+I2
        IPI3 = I+I3
        if (IPI2-IZR) 127,128,127
  127       if (I-IZR) 130,131,130
  128       call INXCC (I,IR,IDXC,NC)
        call INXCB (IPI2,IR,IP2,NP2)
        call INXCB (IPI1,IRM1,IP1,NP1)
        call INXCB (IPI3,IRM1,IP3,NP3)
        call PRDCT (NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM,W2,W2,M, &
                    AM,BM,CM,WD,WW,WU)
        call PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W2,W2,M,AM, &
                    BM,CM,WD,WW,WU)
        DO 129 J=1,M
           Y(J,I) = Y(J,I)-W2(J)
  129       CONTINUE
        IZR = I
        go to 131
  130    CONTINUE
  131 CONTINUE
!
! BEGIN BACK SUBSTITUTION PHASE
!
  132 DO 144 LL=1,K
     L = K-LL+1
     IR = L-1
     IRM1 = IR-1
     I2 = 2**IR
     I1 = I2/2
     I4 = I2+I2
     IFD = IF-I2
     DO 143 I=I2,IFD,I4
        if (I-NM) 133,133,143
  133       IMI1 = I-I1
        IMI2 = I-I2
        IPI1 = I+I1
        IPI2 = I+I2
        call INXCA (I,IR,IDXA,NA)
        call INXCC (I,IR,IDXC,NC)
        call INXCB (I,IR,IZ,NZ)
        call INXCB (IMI1,IRM1,IM1,NM1)
        call INXCB (IPI1,IRM1,IP1,NP1)
        if (I-I2) 134,134,136
  134       DO 135 J=1,M
           W1(J) = (0.,0.)
  135       CONTINUE
        go to 137
  136       call PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),Y(1,IMI2), &
                    W1,M,AM,BM,CM,WD,WW,WU)
  137       if (IPI2-NM) 140,140,138
  138       DO 139 J=1,M
           W2(J) = (0.,0.)
  139       CONTINUE
        go to 141
  140       call PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),Y(1,IPI2), &
                    W2,M,AM,BM,CM,WD,WW,WU)
  141       DO 142 J=1,M
           W1(J) = Y(J,I)+W1(J)+W2(J)
  142       CONTINUE
        call PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,Y(1,I), &
                    M,AM,BM,CM,WD,WW,WU)
  143    CONTINUE
  144 CONTINUE
  return
end
