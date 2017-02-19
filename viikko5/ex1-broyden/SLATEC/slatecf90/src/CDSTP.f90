subroutine CDSTP (EPS, F, FA, HMAX, IMPL, IERROR, JACOBN, MATDIM, &
     MAXORD, MINT, MITER, ML, MU, N, NDE, YWT, UROUND, USERS, AVGH, &
     AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD, NFE, NJE, NQUSED, &
     NSTEP, T, Y, YH, A, CONVRG, DFDY, EL, FAC, HOLD, IPVT, JSTATE, &
     JSTEPL, NQ, NWAIT, RC, RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG, &
     MTRSV, MXRDSV)
!
!! CDSTP performs one step of the integration of an initial value problem ...
!  for a system of ordinary differential equations.
!
!***LIBRARY   SLATEC (SDRIVE)
!***TYPE      COMPLEX (SDSTP-S, DDSTP-D, CDSTP-C)
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***DESCRIPTION
!
!  Communication with CDSTP is done with the following variables:
!
!    YH      An N by MAXORD+1 array containing the dependent variables
!              and their scaled derivatives.  MAXORD, the maximum order
!              used, is currently 12 for the Adams methods and 5 for the
!              Gear methods.  YH(I,J+1) contains the J-th derivative of
!              Y(I), scaled by H**J/factorial(J).  Only Y(I),
!              1  <=  I  <=  N, need be set by the calling program on
!              the first entry.  The YH array should not be altered by
!              the calling program.  When referencing YH as a
!              2-dimensional array, use a column length of N, as this is
!              the value used in CDSTP.
!    DFDY    A block of locations used for partial derivatives if MITER
!              is not 0.  If MITER is 1 or 2 its length must be at least
!              N*N.  If MITER is 4 or 5 its length must be at least
!              (2*ML+MU+1)*N.
!    YWT     An array of N locations used in convergence and error tests
!    SAVE1
!    SAVE2   Arrays of length N used for temporary storage.
!    IPVT    An integer array of length N used by the linear system
!              solvers for the storage of row interchange information.
!    A       A block of locations used to store the matrix A, when using
!              the implicit method.  If IMPL is 1, A is a MATDIM by N
!              array.  If MITER is 1 or 2 MATDIM is N, and if MITER is 4
!              or 5 MATDIM is 2*ML+MU+1.  If IMPL is 2 its length is N.
!              If IMPL is 3, A is a MATDIM by NDE array.
!    JTASK   An integer used on input.
!              It has the following values and meanings:
!                  ==  0  Perform the first step.  This value enables
!                         the subroutine to initialize itself.
!                 >  0  Take a new step continuing from the last.
!                         Assumes the last step was successful and
!                         user has not changed any parameters.
!                  <  0  Take a new step with a new value of H and/or
!                         MINT and/or MITER.
!    JSTATE  A completion code with the following meanings:
!                1  The step was successful.
!                2  A solution could not be obtained with H  /=  0.
!                3  A solution was not obtained in MXTRY attempts.
!                4  For IMPL  /=  0, the matrix A is singular.
!              On a return with JSTATE  >  1, the values of T and
!              the YH array are as of the beginning of the last
!              step, and H is the last step size attempted.
!
!***ROUTINES CALLED  CDCOR, CDCST, CDNTL, CDPSC, CDPST, CDSCL, SCNRM2
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   900329  Initial submission to SLATEC.
!***END PROLOGUE  CDSTP
  EXTERNAL F, JACOBN, FA, USERS
  INTEGER I, IERROR, IMPL, IPVT(*), ISWFLG, ITER, J, JSTATE, JSTEPL, &
          JTASK, MATDIM, MAXORD, MINT, MITER, ML, MNTOLD, MTROLD, &
          MTRSV, MU, MXFAIL, MXITER, MXRDSV, MXTRY, N, NDE, NDJSTP, &
          NFAIL, NFE, NJE, NQ, NQUSED, NSTEP, NSV, NTRY, NWAIT
  COMPLEX A(MATDIM,*), DFDY(MATDIM,*), FAC(*), SAVE1(*), SAVE2(*), &
          Y(*), YH(N,*), YWT(*)
  REAL AVGH, AVGORD, BIAS1, BIAS2, BIAS3, BND, CTEST, D, DENOM, D1, &
       EL(13,12), EPS, ERDN, ERUP, ETEST, H, HMAX, HN, HOLD, HS, &
       HUSED, NUMER, RC, RCTEST, RH, RH1, RH2, RH3, RMAX, RMFAIL, &
       RMNORM, SCNRM2, T, TOLD, TQ(3,12), TREND, TRSHLD, UROUND, &
       Y0NRM
  LOGICAL CONVRG, EVALFA, EVALJC, IER, SWITCH
  PARAMETER(BIAS1 = 1.3E0, BIAS2 = 1.2E0, BIAS3 = 1.4E0, MXFAIL = 3, &
            MXITER = 3, MXTRY = 50, RCTEST = .3E0, RMFAIL = 2.E0, &
            RMNORM = 10.E0, TRSHLD = 1.E0)
  PARAMETER (NDJSTP = 10)
  DATA IER /.FALSE./
!***FIRST EXECUTABLE STATEMENT  CDSTP
  NSV = N
  BND = 0.E0
  SWITCH = .FALSE.
  NTRY = 0
  TOLD = T
  NFAIL = 0
  if (JTASK  <=  0) THEN
    call CDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM, &
                MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T, &
                UROUND, USERS, Y, YWT,  H, MNTOLD, MTROLD, NFE, RC, &
                YH,  A, CONVRG, EL, FAC, IER, IPVT, NQ, NWAIT, RH, &
                RMAX, SAVE2, TQ, TREND, ISWFLG, JSTATE)
    if (N  ==  0) go to 440
    if (H  ==  0.E0) go to 400
    if (IER) go to 420
  end if
 100  NTRY = NTRY + 1
  if (NTRY  >  MXTRY) go to 410
  T = T + H
  call CDPSC (1, N, NQ,  YH)
  EVALJC = (((ABS(RC - 1.E0)  >  RCTEST) .OR. &
    (NSTEP  >=  JSTEPL + NDJSTP)) .AND. (MITER  /=  0))
  EVALFA = .NOT. EVALJC
!
 110  ITER = 0
  DO 115 I = 1,N
 115    Y(I) = YH(I,1)
  call F (N, T, Y, SAVE2)
  if (N  ==  0) THEN
    JSTATE = 6
    go to 430
  end if
  NFE = NFE + 1
  if (EVALJC .OR. IER) THEN
    call CDPST (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML, &
                MU, N, NDE, NQ, SAVE2, T, USERS, Y, YH, YWT, UROUND, &
                NFE, NJE,  A, DFDY, FAC, IER, IPVT, SAVE1, ISWFLG, &
                BND, JSTATE)
    if (N  ==  0) go to 430
    if (IER) go to 160
    CONVRG = .FALSE.
    RC = 1.E0
    JSTEPL = NSTEP
  end if
  DO 125 I = 1,N
 125    SAVE1(I) = 0.E0
!                      Up to MXITER corrector iterations are taken.
!                      Convergence is tested by requiring the r.m.s.
!                      norm of changes to be less than EPS.  The sum of
!                      the corrections is accumulated in the vector
!                      SAVE1(I).  It is approximately equal to the L-th
!                      derivative of Y multiplied by
!                      H**L/(factorial(L-1)*EL(L,NQ)), and is thus
!                      proportional to the actual errors to the lowest
!                      power of H present (H**L).  The YH array is not
!                      altered in the correction loop.  The norm of the
!                      iterate difference is stored in D.  If
!                      ITER  >  0, an estimate of the convergence rate
!                      constant is stored in TREND, and this is used in
!                      the convergence test.
!
 130  call CDCOR (DFDY, EL, FA, H, IERROR, IMPL, IPVT, MATDIM, MITER, &
              ML, MU, N, NDE, NQ, T, USERS, Y, YH, YWT,  EVALFA, &
              SAVE1, SAVE2,  A, D, JSTATE)
    if (N  ==  0) go to 430
  if (ISWFLG  ==  3 .AND. MINT  ==  1) THEN
    if (ITER  ==  0) THEN
      NUMER = SCNRM2(N, SAVE1, 1)
      DO 132 I = 1,N
 132        DFDY(1,I) = SAVE1(I)
      Y0NRM = SCNRM2(N, YH, 1)
    ELSE
      DENOM = NUMER
      DO 134 I = 1,N
 134        DFDY(1,I) = SAVE1(I) - DFDY(1,I)
      NUMER = SCNRM2(N, DFDY, MATDIM)
      if (EL(1,NQ)*NUMER  <=  100.E0*UROUND*Y0NRM) THEN
        if (RMAX  ==  RMFAIL) THEN
          SWITCH = .TRUE.
          go to 170
        end if
      end if
      DO 136 I = 1,N
 136        DFDY(1,I) = SAVE1(I)
      if (DENOM  /=  0.E0) &
      BND = MAX(BND, NUMER/(DENOM*ABS(H)*EL(1,NQ)))
    end if
  end if
  if (ITER  >  0) TREND = MAX(.9E0*TREND, D/D1)
  D1 = D
  CTEST = MIN(2.E0*TREND, 1.E0)*D
  if (CTEST  <=  EPS) go to 170
  ITER = ITER + 1
  if (ITER  <  MXITER) THEN
    DO 140 I = 1,N
 140      Y(I) = YH(I,1) + EL(1,NQ)*SAVE1(I)
    call F (N, T, Y, SAVE2)
    if (N  ==  0) THEN
      JSTATE = 6
      go to 430
    end if
    NFE = NFE + 1
    go to 130
  end if
!                     The corrector iteration failed to converge in
!                     MXITER tries.  If partials are involved but are
!                     not up to date, they are reevaluated for the next
!                     try.  Otherwise the YH array is retracted to its
!                     values before prediction, and H is reduced, if
!                     possible.  If not, a no-convergence exit is taken.
  if (CONVRG) THEN
    EVALJC = .TRUE.
    EVALFA = .FALSE.
    go to 110
  end if
 160  T = TOLD
  call CDPSC (-1, N, NQ,  YH)
  NWAIT = NQ + 2
  if (JTASK  /=  0 .AND. JTASK  /=  2) RMAX = RMFAIL
  if (ITER  ==  0) THEN
    RH = .3E0
  ELSE
    RH = .9E0*(EPS/CTEST)**(.2E0)
  end if
  if (RH*H  ==  0.E0) go to 400
  call CDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
  go to 100
!                          The corrector has converged.  CONVRG is set
!                          to .TRUE. if partial derivatives were used,
!                          to indicate that they may need updating on
!                          subsequent steps.  The error test is made.
 170  CONVRG = (MITER  /=  0)
  if (IERROR  ==  1 .OR. IERROR  ==  5) THEN
    DO 180 I = 1,NDE
 180      SAVE2(I) = SAVE1(I)/YWT(I)
  ELSE
    DO 185 I = 1,NDE
 185      SAVE2(I) = SAVE1(I)/MAX(ABS(Y(I)), ABS(YWT(I)))
  end if
  ETEST = SCNRM2(NDE, SAVE2, 1)/(TQ(2,NQ)*SQRT(REAL(NDE)))
!
!                           The error test failed.  NFAIL keeps track of
!                           multiple failures.  Restore T and the YH
!                           array to their previous values, and prepare
!                           to try the step again.  Compute the optimum
!                           step size for this or one lower order.
  if (ETEST  >  EPS) THEN
    T = TOLD
    call CDPSC (-1, N, NQ,  YH)
    NFAIL = NFAIL + 1
    if (NFAIL  <  MXFAIL .OR. NQ  ==  1) THEN
      if (JTASK  /=  0 .AND. JTASK  /=  2) RMAX = RMFAIL
      RH2 = 1.E0/(BIAS2*(ETEST/EPS)**(1.E0/(NQ+1)))
      if (NQ  >  1) THEN
        if (IERROR  ==  1 .OR. IERROR  ==  5) THEN
          DO 190 I = 1,NDE
 190            SAVE2(I) = YH(I,NQ+1)/YWT(I)
        ELSE
          DO 195 I = 1,NDE
 195            SAVE2(I) = YH(I,NQ+1)/MAX(ABS(Y(I)), ABS(YWT(I)))
        end if
        ERDN = SCNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(REAL(NDE)))
        RH1 = 1.E0/MAX(1.E0, BIAS1*(ERDN/EPS)**(1.E0/NQ))
        if (RH2  <  RH1) THEN
          NQ = NQ - 1
          RC = RC*EL(1,NQ)/EL(1,NQ+1)
          RH = RH1
        ELSE
          RH = RH2
        end if
      ELSE
        RH = RH2
      end if
      NWAIT = NQ + 2
      if (RH*H  ==  0.E0) go to 400
      call CDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
      go to 100
    end if
!                Control reaches this section if the error test has
!                failed MXFAIL or more times.  It is assumed that the
!                derivatives that have accumulated in the YH array have
!                errors of the wrong order.  Hence the first derivative
!                is recomputed, the order is set to 1, and the step is
!                retried.
    NFAIL = 0
    JTASK = 2
    DO 215 I = 1,N
 215      Y(I) = YH(I,1)
    call CDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM, &
                MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T, &
                UROUND, USERS, Y, YWT,  H, MNTOLD, MTROLD, NFE, RC, &
                YH,  A, CONVRG, EL, FAC, IER, IPVT, NQ, NWAIT, RH, &
                RMAX, SAVE2, TQ, TREND, ISWFLG, JSTATE)
    RMAX = RMNORM
    if (N  ==  0) go to 440
    if (H  ==  0.E0) go to 400
    if (IER) go to 420
    go to 100
  end if
!                          After a successful step, update the YH array.
  NSTEP = NSTEP + 1
  HUSED = H
  NQUSED = NQ
  AVGH = ((NSTEP-1)*AVGH + H)/NSTEP
  AVGORD = ((NSTEP-1)*AVGORD + NQ)/NSTEP
  DO 230 J = 1,NQ+1
    DO 230 I = 1,N
 230      YH(I,J) = YH(I,J) + EL(J,NQ)*SAVE1(I)
  DO 235 I = 1,N
 235    Y(I) = YH(I,1)
!                                          If ISWFLG is 3, consider
!                                          changing integration methods.
  if (ISWFLG  ==  3) THEN
    if (BND  /=  0.E0) THEN
      if (MINT  ==  1 .AND. NQ  <=  5) THEN
        HN = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.E0/(NQ+1)))
        HN = MIN(HN, 1.E0/(2.E0*EL(1,NQ)*BND))
        HS = ABS(H)/MAX(UROUND, &
        (ETEST/(EPS*EL(NQ+1,1)))**(1.E0/(NQ+1)))
        if (HS  >  1.2E0*HN) THEN
          MINT = 2
          MNTOLD = MINT
          MITER = MTRSV
          MTROLD = MITER
          MAXORD = MIN(MXRDSV, 5)
          RC = 0.E0
          RMAX = RMNORM
          TREND = 1.E0
          call CDCST (MAXORD, MINT, ISWFLG, EL, TQ)
          NWAIT = NQ + 2
        end if
      ELSE if (MINT  ==  2) THEN
        HS = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.E0/(NQ+1)))
        HN = ABS(H)/MAX(UROUND, &
        (ETEST*EL(NQ+1,1)/EPS)**(1.E0/(NQ+1)))
        HN = MIN(HN, 1.E0/(2.E0*EL(1,NQ)*BND))
        if (HN  >=  HS) THEN
          MINT = 1
          MNTOLD = MINT
          MITER = 0
          MTROLD = MITER
          MAXORD = MIN(MXRDSV, 12)
          RMAX = RMNORM
          TREND = 1.E0
          CONVRG = .FALSE.
          call CDCST (MAXORD, MINT, ISWFLG, EL, TQ)
          NWAIT = NQ + 2
        end if
      end if
    end if
  end if
  if (SWITCH) THEN
    MINT = 2
    MNTOLD = MINT
    MITER = MTRSV
    MTROLD = MITER
    MAXORD = MIN(MXRDSV, 5)
    NQ = MIN(NQ, MAXORD)
    RC = 0.E0
    RMAX = RMNORM
    TREND = 1.E0
    call CDCST (MAXORD, MINT, ISWFLG, EL, TQ)
    NWAIT = NQ + 2
  end if
!                           Consider changing H if NWAIT = 1.  Otherwise
!                           decrease NWAIT by 1.  If NWAIT is then 1 and
!                           NQ < MAXORD, then SAVE1 is saved for use in
!                           a possible order increase on the next step.
!
  if (JTASK  ==  0 .OR. JTASK  ==  2) THEN
    RH = 1.E0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.E0/(NQ+1)))
    if (RH > TRSHLD) call CDSCL (HMAX, N, NQ, RMAX, H, RC, RH, YH)
  ELSE if (NWAIT  >  1) THEN
    NWAIT = NWAIT - 1
    if (NWAIT  ==  1 .AND. NQ  <  MAXORD) THEN
      DO 250 I = 1,NDE
 250        YH(I,MAXORD+1) = SAVE1(I)
    end if
!             If a change in H is considered, an increase or decrease in
!             order by one is considered also.  A change in H is made
!             only if it is by a factor of at least TRSHLD.  Factors
!             RH1, RH2, and RH3 are computed, by which H could be
!             multiplied at order NQ - 1, order NQ, or order NQ + 1,
!             respectively.  The largest of these is determined and the
!             new order chosen accordingly.  If the order is to be
!             increased, we compute one additional scaled derivative.
!             If there is a change of order, reset NQ and the
!             coefficients.  In any case H is reset according to RH and
!             the YH array is rescaled.
  ELSE
    if (NQ  ==  1) THEN
      RH1 = 0.E0
    ELSE
      if (IERROR  ==  1 .OR. IERROR  ==  5) THEN
        DO 270 I = 1,NDE
 270          SAVE2(I) = YH(I,NQ+1)/YWT(I)
      ELSE
        DO 275 I = 1,NDE
 275          SAVE2(I) = YH(I,NQ+1)/MAX(ABS(Y(I)), ABS(YWT(I)))
      end if
      ERDN = SCNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(REAL(NDE)))
      RH1 = 1.E0/MAX(UROUND, BIAS1*(ERDN/EPS)**(1.E0/NQ))
    end if
    RH2 = 1.E0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.E0/(NQ+1)))
    if (NQ  ==  MAXORD) THEN
      RH3 = 0.E0
    ELSE
      if (IERROR  ==  1 .OR. IERROR  ==  5) THEN
        DO 290 I = 1,NDE
 290          SAVE2(I) = (SAVE1(I) - YH(I,MAXORD+1))/YWT(I)
      ELSE
        DO 295 I = 1,NDE
          SAVE2(I) = (SAVE1(I) - YH(I,MAXORD+1))/ &
          MAX(ABS(Y(I)), ABS(YWT(I)))
 295          CONTINUE
      end if
      ERUP = SCNRM2(NDE, SAVE2, 1)/(TQ(3,NQ)*SQRT(REAL(NDE)))
      RH3 = 1.E0/MAX(UROUND, BIAS3*(ERUP/EPS)**(1.E0/(NQ+2)))
    end if
    if (RH1  >  RH2 .AND. RH1  >=  RH3) THEN
      RH = RH1
      if (RH  <=  TRSHLD) go to 380
      NQ = NQ - 1
      RC = RC*EL(1,NQ)/EL(1,NQ+1)
    ELSE if (RH2  >=  RH1 .AND. RH2  >=  RH3) THEN
      RH = RH2
      if (RH  <=  TRSHLD) go to 380
    ELSE
      RH = RH3
      if (RH  <=  TRSHLD) go to 380
      DO 360 I = 1,N
 360        YH(I,NQ+2) = SAVE1(I)*EL(NQ+1,NQ)/(NQ+1)
      NQ = NQ + 1
      RC = RC*EL(1,NQ)/EL(1,NQ-1)
    end if
    if (ISWFLG  ==  3 .AND. MINT  ==  1) THEN
      if (BND /= 0.E0) RH = MIN(RH, 1.E0/(2.E0*EL(1,NQ)*BND*ABS(H)))
    end if
    call CDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
    RMAX = RMNORM
 380    NWAIT = NQ + 2
  end if
!               All returns are made through this section.  H is saved
!               in HOLD to allow the caller to change H on the next step
  JSTATE = 1
  HOLD = H
  return
!
 400  JSTATE = 2
  HOLD = H
  DO 405 I = 1,N
 405    Y(I) = YH(I,1)
  return
!
 410  JSTATE = 3
  HOLD = H
  return
!
 420  JSTATE = 4
  HOLD = H
  return
!
 430  T = TOLD
  call CDPSC (-1, NSV, NQ,  YH)
  DO 435 I = 1,NSV
 435    Y(I) = YH(I,1)
 440  HOLD = H
  return
end
