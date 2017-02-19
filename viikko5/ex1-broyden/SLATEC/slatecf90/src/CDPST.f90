subroutine CDPST (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML, &
     MU, N, NDE, NQ, SAVE2, T, USERS, Y, YH, YWT, UROUND, NFE, NJE, &
     A, DFDY, FAC, IER, IPVT, SAVE1, ISWFLG, BND, JSTATE)
!
!! CDPST evaluates the Jacobian matrix of the right hand side ...
!  of the differential equations.
!
!***LIBRARY   SLATEC (SDRIVE)
!***TYPE      COMPLEX (SDPST-S, DDPST-D, CDPST-C)
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***DESCRIPTION
!
!  If MITER is 1, 2, 4, or 5, the matrix
!  P = I - L(0)*H*Jacobian is stored in DFDY and subjected to LU
!  decomposition, with the results also stored in DFDY.
!
!***ROUTINES CALLED  CGBFA, CGEFA, SCNRM2
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   900329  Initial submission to SLATEC.
!***END PROLOGUE  CDPST
  INTEGER I, IFLAG, IMAX, IMPL, INFO, ISWFLG, J, J2, JSTATE, K, &
          MATDIM, MITER, ML, MU, MW, N, NDE, NFE, NJE, NQ
  COMPLEX A(MATDIM,*), CFCTR, DFDY(MATDIM,*), DY, FAC(*), SAVE1(*), &
          SAVE2(*), Y(*), YH(N,*), YJ, YS, YWT(*)
  REAL BL, BND, BP, BR, BU, DFDYMX, DIFF, EL(13,12), FACMAX, FACMIN, &
       FACTOR, H, SCALE, SCNRM2, T, UROUND, ZMAX, ZMIN
  INTEGER IPVT(*)
  LOGICAL IER
  PARAMETER(FACMAX = .5E0, BU = 0.5E0)
!***FIRST EXECUTABLE STATEMENT  CDPST
  NJE = NJE + 1
  IER = .FALSE.
  if (MITER  ==  1 .OR. MITER  ==  2) THEN
    if (MITER  ==  1) THEN
      call JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
      if (N  ==  0) THEN
        JSTATE = 8
        return
      end if
      if (ISWFLG  ==  3) BND = SCNRM2(N*N, DFDY, 1)
      FACTOR = -EL(1,NQ)*H
      DO 110 J = 1,N
        DO 110 I = 1,N
 110          DFDY(I,J) = FACTOR*DFDY(I,J)
    ELSE if (MITER  ==  2) THEN
      BR = UROUND**(.875E0)
      BL = UROUND**(.75E0)
      BP = UROUND**(-.15E0)
      FACMIN = UROUND**(.78E0)
      DO 170 J = 1,N
        if (ABS(Y(J))  >  ABS(YWT(J))) THEN
          YS = Y(J)
        ELSE
          YS = YWT(J)
        end if
 120        DY = FAC(J)*YS
        if (DY  ==  0.E0) THEN
          if (REAL(FAC(J))  <  FACMAX) THEN
            FAC(J) = MIN(100.E0*REAL(FAC(J)), FACMAX)
            go to 120
          ELSE
            DY = YS
          end if
        end if
        DY = (Y(J) + DY) - Y(J)
        YJ = Y(J)
        Y(J) = Y(J) + DY
        call F (N, T, Y, SAVE1)
        if (N  ==  0) THEN
          JSTATE = 6
          return
        end if
        Y(J) = YJ
        CFCTR = -EL(1,NQ)*H/DY
        DO 140 I = 1,N
 140          DFDY(I,J) = (SAVE1(I) - SAVE2(I))*CFCTR
!                                                                 Step 1
        DIFF = ABS(SAVE2(1) - SAVE1(1))
        IMAX = 1
        DO 150 I = 2,N
          if (ABS(SAVE2(I) - SAVE1(I))  >  DIFF) THEN
            IMAX = I
            DIFF = ABS(SAVE2(I) - SAVE1(I))
          end if
 150          CONTINUE
!                                                                 Step 2
        if (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))  >  0.E0) THEN
          SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
!                                                                 Step 3
          if (DIFF  >  BU*SCALE) THEN
            FAC(J) = MAX(FACMIN, REAL(FAC(J))*.5E0)
          ELSE if (BR*SCALE  <=  DIFF .AND. DIFF  <=  BL*SCALE) THEN
            FAC(J) = MIN(REAL(FAC(J))*2.E0, FACMAX)
!                                                                 Step 4
          ELSE if (DIFF  <  BR*SCALE) THEN
            FAC(J) = MIN(BP*REAL(FAC(J)), FACMAX)
          end if
        end if
 170        CONTINUE
      if (ISWFLG  ==  3) BND = SCNRM2(N*N, DFDY, 1)/(-EL(1,NQ)*H)
      NFE = NFE + N
    end if
    if (IMPL  ==  0) THEN
      DO 190 I = 1,N
 190        DFDY(I,I) = DFDY(I,I) + 1.E0
    ELSE if (IMPL  ==  1) THEN
      call FA (N, T, Y, A, MATDIM, ML, MU, NDE)
      if (N  ==  0) THEN
        JSTATE = 9
        return
      end if
      DO 210 J = 1,N
        DO 210 I = 1,N
 210          DFDY(I,J) = DFDY(I,J) + A(I,J)
    ELSE if (IMPL  ==  2) THEN
      call FA (N, T, Y, A, MATDIM, ML, MU, NDE)
      if (N  ==  0) THEN
        JSTATE = 9
        return
      end if
      DO 230 I = 1,NDE
 230        DFDY(I,I) = DFDY(I,I) + A(I,1)
    ELSE if (IMPL  ==  3) THEN
      call FA (N, T, Y, A, MATDIM, ML, MU, NDE)
      if (N  ==  0) THEN
        JSTATE = 9
        return
      end if
      DO 220 J = 1,NDE
        DO 220 I = 1,NDE
 220          DFDY(I,J) = DFDY(I,J) + A(I,J)
    end if
    call CGEFA (DFDY, MATDIM, N, IPVT, INFO)
    if (INFO  /=  0) IER = .TRUE.
  ELSE if (MITER  ==  4 .OR. MITER  ==  5) THEN
    if (MITER  ==  4) THEN
      call JACOBN (N, T, Y, DFDY(ML+1,1), MATDIM, ML, MU)
      if (N  ==  0) THEN
        JSTATE = 8
        return
      end if
      FACTOR = -EL(1,NQ)*H
      MW = ML + MU + 1
      DO 260 J = 1,N
        DO 260 I = MAX(ML+1, MW+1-J), MIN(MW+N-J, MW+ML)
 260          DFDY(I,J) = FACTOR*DFDY(I,J)
    ELSE if (MITER  ==  5) THEN
      BR = UROUND**(.875E0)
      BL = UROUND**(.75E0)
      BP = UROUND**(-.15E0)
      FACMIN = UROUND**(.78E0)
      MW = ML + MU + 1
      J2 = MIN(MW, N)
      DO 340 J = 1,J2
        DO 290 K = J,N,MW
          if (ABS(Y(K))  >  ABS(YWT(K))) THEN
            YS = Y(K)
          ELSE
            YS = YWT(K)
          end if
 280          DY = FAC(K)*YS
          if (DY  ==  0.E0) THEN
            if (REAL(FAC(K))  <  FACMAX) THEN
              FAC(K) = MIN(100.E0*REAL(FAC(K)), FACMAX)
              go to 280
            ELSE
              DY = YS
            end if
          end if
          DY = (Y(K) + DY) - Y(K)
          DFDY(MW,K) = Y(K)
 290          Y(K) = Y(K) + DY
        call F (N, T, Y, SAVE1)
        if (N  ==  0) THEN
          JSTATE = 6
          return
        end if
        DO 330 K = J,N,MW
          DY = Y(K) - DFDY(MW,K)
          Y(K) = DFDY(MW,K)
          CFCTR = -EL(1,NQ)*H/DY
          DO 300 I = MAX(ML+1, MW+1-K), MIN(MW+N-K, MW+ML)
 300            DFDY(I,K) = CFCTR*(SAVE1(I+K-MW) - SAVE2(I+K-MW))
!                                                                 Step 1
          IMAX = MAX(1, K - MU)
          DIFF = ABS(SAVE2(IMAX) - SAVE1(IMAX))
          DO 310 I = MAX(1, K - MU)+1, MIN(K + ML, N)
            if (ABS(SAVE2(I) - SAVE1(I))  >  DIFF) THEN
              IMAX = I
              DIFF = ABS(SAVE2(I) - SAVE1(I))
            end if
 310            CONTINUE
!                                                                 Step 2
          if (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))  > 0.E0) THEN
            SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
!                                                                 Step 3
            if (DIFF  >  BU*SCALE) THEN
              FAC(J) = MAX(FACMIN, REAL(FAC(J))*.5E0)
            ELSE if (BR*SCALE  <= DIFF .AND. DIFF  <= BL*SCALE) THEN
              FAC(J) = MIN(REAL(FAC(J))*2.E0, FACMAX)
!                                                                 Step 4
            ELSE if (DIFF  <  BR*SCALE) THEN
              FAC(K) = MIN(BP*REAL(FAC(K)), FACMAX)
            end if
          end if
 330          CONTINUE
 340        CONTINUE
      NFE = NFE + J2
    end if
    if (ISWFLG  ==  3) THEN
      DFDYMX = 0.E0
      DO 345 J = 1,N
        DO 345 I = MAX(ML+1, MW+1-J), MIN(MW+N-J, MW+ML)
          ZMAX = MAX(ABS(REAL(DFDY(I,J))), ABS(AIMAG(DFDY(I,J))))
          ZMIN = MIN(ABS(REAL(DFDY(I,J))), ABS(AIMAG(DFDY(I,J))))
          if (ZMAX  /=  0.E0) &
          DFDYMX = MAX(DFDYMX, ZMAX*SQRT(1.E0+ (ZMIN/ZMAX)**2))
 345          CONTINUE
      BND = 0.E0
      if (DFDYMX  /=  0.E0) THEN
        DO 350 J = 1,N
          DO 350 I = MAX(ML+1, MW+1-J), MIN(MW+N-J, MW+ML)
            BND = BND + (REAL(DFDY(I,J))/DFDYMX)**2 + &
            (AIMAG(DFDY(I,J))/DFDYMX)**2
 350            CONTINUE
        BND = DFDYMX*SQRT(BND)/(-EL(1,NQ)*H)
      end if
    end if
    if (IMPL  ==  0) THEN
      DO 360 J = 1,N
 360        DFDY(MW,J) = DFDY(MW,J) + 1.E0
    ELSE if (IMPL  ==  1) THEN
      call FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
      if (N  ==  0) THEN
        JSTATE = 9
        return
      end if
      DO 380 J = 1,N
        DO 380 I = MAX(ML+1, MW+1-J), MIN(MW+N-J, MW+ML)
 380          DFDY(I,J) = DFDY(I,J) + A(I,J)
    ELSE if (IMPL  ==  2) THEN
      call FA (N, T, Y, A, MATDIM, ML, MU, NDE)
      if (N  ==  0) THEN
        JSTATE = 9
        return
      end if
      DO 400 J = 1,NDE
 400        DFDY(MW,J) =  DFDY(MW,J) + A(J,1)
    ELSE if (IMPL  ==  3) THEN
      call FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
      if (N  ==  0) THEN
        JSTATE = 9
        return
      end if
      DO 390 J = 1,NDE
        DO 390 I = MAX(ML+1, MW+1-J), MIN(MW+NDE-J, MW+ML)
 390          DFDY(I,J) = DFDY(I,J) + A(I,J)
    end if
    call CGBFA (DFDY, MATDIM, N, ML, MU, IPVT, INFO)
    if (INFO  /=  0) IER = .TRUE.
  ELSE if (MITER  ==  3) THEN
    IFLAG = 1
    call USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL, &
                N, NDE, IFLAG)
    if (IFLAG  ==  -1) THEN
      IER = .TRUE.
      return
    end if
    if (N  ==  0) THEN
      JSTATE = 10
      return
    end if
  end if
  return
end
