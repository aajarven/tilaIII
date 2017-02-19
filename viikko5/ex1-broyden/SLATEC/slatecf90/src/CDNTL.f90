subroutine CDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM, &
     MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T, UROUND, USERS, &
     Y, YWT, H, MNTOLD, MTROLD, NFE, RC, YH, A, CONVRG, EL, FAC, &
     IER, IPVT, NQ, NWAIT, RH, RMAX, SAVE2, TQ, TREND, ISWFLG, &
     JSTATE)
!
!! CDNTL sets parameters on the first call to CDSTP, on an internal ...
!  restart, or when the user has altered MINT, MITER, and/or H.
!
!***LIBRARY   SLATEC (SDRIVE)
!***TYPE      COMPLEX (SDNTL-S, DDNTL-D, CDNTL-C)
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***DESCRIPTION
!
!  On the first call, the order is set to 1 and the initial derivatives
!  are calculated.  RMAX is the maximum ratio by which H can be
!  increased in one step.  It is initially RMINIT to compensate
!  for the small initial H, but then is normally equal to RMNORM.
!  If a failure occurs (in corrector convergence or error test), RMAX
!  is set at RMFAIL for the next increase.
!  If the caller has changed MINT, or if JTASK = 0, CDCST is called
!  to set the coefficients of the method.  If the caller has changed H,
!  YH must be rescaled.  If H or MINT has been changed, NWAIT is
!  reset to NQ + 2 to prevent further increases in H for that many
!  steps.  Also, RC is reset.  RC is the ratio of new to old values of
!  the coefficient L(0)*H.  If the caller has changed MITER, RC is
!  set to 0 to force the partials to be updated, if partials are used.
!
!***ROUTINES CALLED  CDCST, CDSCL, CGBFA, CGBSL, CGEFA, CGESL, SCNRM2
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   900329  Initial submission to SLATEC.
!***END PROLOGUE  CDNTL
  INTEGER I, IFLAG, IMPL, INFO, ISWFLG, JSTATE, JTASK, MATDIM, &
          MAXORD, MINT, MITER, ML, MNTOLD, MTROLD, MU, N, NDE, NFE, &
          NQ, NWAIT
  COMPLEX A(MATDIM,*), FAC(*), SAVE1(*), SAVE2(*), Y(*), YH(N,*), &
          YWT(*)
  REAL EL(13,12), EPS, H, HMAX, HOLD, OLDL0, RC, RH, RMAX, &
       RMINIT, SCNRM2, SUM, T, TQ(3,12), TREND, UROUND
  INTEGER IPVT(*)
  LOGICAL CONVRG, IER
  PARAMETER(RMINIT = 10000.E0)
!***FIRST EXECUTABLE STATEMENT  CDNTL
  IER = .FALSE.
  if (JTASK  >=  0) THEN
    if (JTASK  ==  0) THEN
      call CDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
      RMAX = RMINIT
    end if
    RC = 0.E0
    CONVRG = .FALSE.
    TREND = 1.E0
    NQ = 1
    NWAIT = 3
    call F (N, T, Y, SAVE2)
    if (N  ==  0) THEN
      JSTATE = 6
      return
    end if
    NFE = NFE + 1
    if (IMPL  /=  0) THEN
      if (MITER  ==  3) THEN
        IFLAG = 0
        call USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL, IMPL, N, &
                    NDE, IFLAG)
        if (IFLAG  ==  -1) THEN
          IER = .TRUE.
          return
        end if
        if (N  ==  0) THEN
          JSTATE = 10
          return
        end if
      ELSE if (IMPL  ==  1) THEN
        if (MITER  ==  1 .OR. MITER  ==  2) THEN
          call FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          if (N  ==  0) THEN
            JSTATE = 9
            return
          end if
          call CGEFA (A, MATDIM, N, IPVT, INFO)
          if (INFO  /=  0) THEN
            IER = .TRUE.
            return
          end if
          call CGESL (A, MATDIM, N, IPVT, SAVE2, 0)
        ELSE if (MITER  ==  4 .OR. MITER  ==  5) THEN
          call FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
          if (N  ==  0) THEN
            JSTATE = 9
            return
          end if
          call CGBFA (A, MATDIM, N, ML, MU, IPVT, INFO)
          if (INFO  /=  0) THEN
            IER = .TRUE.
            return
          end if
          call CGBSL (A, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
        end if
      ELSE if (IMPL  ==  2) THEN
        call FA (N, T, Y, A, MATDIM, ML, MU, NDE)
        if (N  ==  0) THEN
          JSTATE = 9
          return
        end if
        DO 150 I = 1,NDE
          if (A(I,1)  ==  0.E0) THEN
            IER = .TRUE.
            return
          ELSE
            SAVE2(I) = SAVE2(I)/A(I,1)
          end if
 150          CONTINUE
        DO 155 I = NDE+1,N
 155          A(I,1) = 0.E0
      ELSE if (IMPL  ==  3) THEN
        if (MITER  ==  1 .OR. MITER  ==  2) THEN
          call FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          if (N  ==  0) THEN
            JSTATE = 9
            return
          end if
          call CGEFA (A, MATDIM, NDE, IPVT, INFO)
          if (INFO  /=  0) THEN
            IER = .TRUE.
            return
          end if
          call CGESL (A, MATDIM, NDE, IPVT, SAVE2, 0)
        ELSE if (MITER  ==  4 .OR. MITER  ==  5) THEN
          call FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
          if (N  ==  0) THEN
            JSTATE = 9
            return
          end if
          call CGBFA (A, MATDIM, NDE, ML, MU, IPVT, INFO)
          if (INFO  /=  0) THEN
            IER = .TRUE.
            return
          end if
          call CGBSL (A, MATDIM, NDE, ML, MU, IPVT, SAVE2, 0)
        end if
      end if
    end if
    DO 170 I = 1,NDE
 170      SAVE1(I) = SAVE2(I)/MAX(1.E0, ABS(YWT(I)))
    SUM = SCNRM2(NDE, SAVE1, 1)/SQRT(REAL(NDE))
    if (SUM  >  EPS/ABS(H)) H = SIGN(EPS/SUM, H)
    DO 180 I = 1,N
 180      YH(I,2) = H*SAVE2(I)
    if (MITER  ==  2 .OR. MITER  ==  5 .OR. ISWFLG  ==  3) THEN
      DO 20 I = 1,N
 20         FAC(I) = SQRT(UROUND)
    end if
  ELSE
    if (MITER  /=  MTROLD) THEN
      MTROLD = MITER
      RC = 0.E0
      CONVRG = .FALSE.
    end if
    if (MINT  /=  MNTOLD) THEN
      MNTOLD = MINT
      OLDL0 = EL(1,NQ)
      call CDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
      RC = RC*EL(1,NQ)/OLDL0
      NWAIT = NQ + 2
    end if
    if (H  /=  HOLD) THEN
      NWAIT = NQ + 2
      RH = H/HOLD
      call CDSCL (HMAX, N, NQ, RMAX,  HOLD, RC, RH, YH)
    end if
  end if
  return
end
