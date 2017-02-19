subroutine DDCOR (DFDY, EL, FA, H, IERROR, IMPL, IPVT, MATDIM, &
     MITER, ML, MU, N, NDE, NQ, T, USERS, Y, YH, YWT, EVALFA, SAVE1, &
     SAVE2, A, D, JSTATE)
!
!! DDCOR computes corrections to the Y array for DDRIVE.
!
!***LIBRARY   SLATEC (SDRIVE)
!***TYPE      DOUBLE PRECISION (SDCOR-S, DDCOR-D, CDCOR-C)
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***DESCRIPTION
!
!  In the case of functional iteration, update Y directly from the
!  result of the last call to F.
!  In the case of the chord method, compute the corrector error and
!  solve the linear system with that as right hand side and DFDY as
!  coefficient matrix, using the LU decomposition if MITER is 1, 2, 4,
!  or 5.
!
!***ROUTINES CALLED  DGBSL, DGESL, DNRM2
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   900329  Initial submission to SLATEC.
!***END PROLOGUE  DDCOR
  INTEGER I, IERROR, IFLAG, IMPL, J, JSTATE, MATDIM, MITER, ML, MU, &
          MW, N, NDE, NQ
  DOUBLE PRECISION A(MATDIM,*), D, DFDY(MATDIM,*), EL(13,12), H, &
       SAVE1(*), SAVE2(*), DNRM2, T, Y(*), YH(N,*), YWT(*)
  INTEGER IPVT(*)
  LOGICAL EVALFA
!***FIRST EXECUTABLE STATEMENT  DDCOR
  if (MITER  ==  0) THEN
    if (IERROR  ==  1 .OR. IERROR  ==  5) THEN
      DO 100 I = 1,N
 100        SAVE1(I) = (H*SAVE2(I) - YH(I,2) - SAVE1(I))/YWT(I)
    ELSE
      DO 102 I = 1,N
        SAVE1(I) = (H*SAVE2(I) - YH(I,2) - SAVE1(I))/ &
        MAX(ABS(Y(I)), YWT(I))
 102        CONTINUE
    end if
    D = DNRM2(N, SAVE1, 1)/SQRT(DBLE(N))
    DO 105 I = 1,N
 105      SAVE1(I) = H*SAVE2(I) - YH(I,2)
  ELSE if (MITER  ==  1 .OR. MITER  ==  2) THEN
    if (IMPL  ==  0) THEN
      DO 130 I = 1,N
 130        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
    ELSE if (IMPL  ==  1) THEN
      if (EVALFA) THEN
        call FA (N, T, Y, A, MATDIM, ML, MU, NDE)
        if (N  ==  0) THEN
          JSTATE = 9
          return
        end if
      ELSE
        EVALFA = .TRUE.
      end if
      DO 150 I = 1,N
 150        SAVE2(I) = H*SAVE2(I)
      DO 160 J = 1,N
        DO 160 I = 1,N
 160          SAVE2(I) = SAVE2(I) - A(I,J)*(YH(J,2) + SAVE1(J))
    ELSE if (IMPL  ==  2) THEN
      if (EVALFA) THEN
        call FA (N, T, Y, A, MATDIM, ML, MU, NDE)
        if (N  ==  0) THEN
          JSTATE = 9
          return
        end if
      ELSE
        EVALFA = .TRUE.
      end if
      DO 180 I = 1,N
 180        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
    ELSE if (IMPL  ==  3) THEN
      if (EVALFA) THEN
        call FA (N, T, Y, A, MATDIM, ML, MU, NDE)
        if (N  ==  0) THEN
          JSTATE = 9
          return
        end if
      ELSE
        EVALFA = .TRUE.
      end if
      DO 140 I = 1,N
 140        SAVE2(I) = H*SAVE2(I)
      DO 170 J = 1,NDE
        DO 170 I = 1,NDE
 170          SAVE2(I) = SAVE2(I) - A(I,J)*(YH(J,2) + SAVE1(J))
    end if
    call DGESL (DFDY, MATDIM, N, IPVT, SAVE2, 0)
    if (IERROR  ==  1 .OR. IERROR  ==  5) THEN
      DO 200 I = 1,N
        SAVE1(I) = SAVE1(I) + SAVE2(I)
 200        SAVE2(I) = SAVE2(I)/YWT(I)
    ELSE
      DO 205 I = 1,N
        SAVE1(I) = SAVE1(I) + SAVE2(I)
 205        SAVE2(I) = SAVE2(I)/MAX(ABS(Y(I)), YWT(I))
    end if
    D = DNRM2(N, SAVE2, 1)/SQRT(DBLE(N))
  ELSE if (MITER  ==  4 .OR. MITER  ==  5) THEN
    if (IMPL  ==  0) THEN
      DO 230 I = 1,N
 230        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
    ELSE if (IMPL  ==  1) THEN
      if (EVALFA) THEN
        call FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
        if (N  ==  0) THEN
          JSTATE = 9
          return
        end if
      ELSE
        EVALFA = .TRUE.
      end if
      DO 250 I = 1,N
 250        SAVE2(I) = H*SAVE2(I)
      MW = ML + 1 + MU
      DO 260 J = 1,N
        DO 260 I = MAX(ML+1, MW+1-J), MIN(MW+N-J, MW+ML)
          SAVE2(I+J-MW) = SAVE2(I+J-MW) &
                          - A(I,J)*(YH(J,2) + SAVE1(J))
 260        CONTINUE
    ELSE if (IMPL  ==  2) THEN
      if (EVALFA) THEN
        call FA (N, T, Y, A, MATDIM, ML, MU, NDE)
        if (N  ==  0) THEN
          JSTATE = 9
          return
        end if
      ELSE
        EVALFA = .TRUE.
      end if
      DO 280 I = 1,N
 280        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
    ELSE if (IMPL  ==  3) THEN
      if (EVALFA) THEN
        call FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
        if (N  ==  0) THEN
          JSTATE = 9
          return
        end if
      ELSE
        EVALFA = .TRUE.
      end if
      DO 270 I = 1,N
 270        SAVE2(I) = H*SAVE2(I)
      MW = ML + 1 + MU
      DO 290 J = 1,NDE
        DO 290 I = MAX(ML+1, MW+1-J), MIN(MW+NDE-J, MW+ML)
          SAVE2(I+J-MW) = SAVE2(I+J-MW) &
                          - A(I,J)*(YH(J,2) + SAVE1(J))
 290        CONTINUE
    end if
    call DGBSL (DFDY, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
    if (IERROR  ==  1 .OR. IERROR  ==  5) THEN
      DO 300 I = 1,N
        SAVE1(I) = SAVE1(I) + SAVE2(I)
 300        SAVE2(I) = SAVE2(I)/YWT(I)
    ELSE
      DO 305 I = 1,N
        SAVE1(I) = SAVE1(I) + SAVE2(I)
 305        SAVE2(I) = SAVE2(I)/MAX(ABS(Y(I)), YWT(I))
    end if
    D = DNRM2(N, SAVE2, 1)/SQRT(DBLE(N))
  ELSE if (MITER  ==  3) THEN
    IFLAG = 2
    call USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL, &
                N, NDE, IFLAG)
    if (N  ==  0) THEN
      JSTATE = 10
      return
    end if
    if (IERROR  ==  1 .OR. IERROR  ==  5) THEN
      DO 320 I = 1,N
        SAVE1(I) = SAVE1(I) + SAVE2(I)
 320        SAVE2(I) = SAVE2(I)/YWT(I)
    ELSE
      DO 325 I = 1,N
        SAVE1(I) = SAVE1(I) + SAVE2(I)
 325        SAVE2(I) = SAVE2(I)/MAX(ABS(Y(I)), YWT(I))
    end if
    D = DNRM2(N, SAVE2, 1)/SQRT(DBLE(N))
  end if
  return
end
