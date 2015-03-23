      PROGRAM HWIGPR
C---Debug and test program for new versions of Herwig
      INCLUDE 'HERWIG65.INC'
      EXTERNAL HWUDAT
      LOGICAL RSUD
      INTEGER JMIN,JMAX,J,N,JPRO,JPROC(35),NRS(2)
      DATA JPROC/1500,105,200,250,300,400,500,550,1399,1450,1612,1705,
     $     1706,1800,1999,2000,2100,2150,2200,2399,2400,2599,2699,2799,
     $     2800,2810,2820,2915,5000,5500,8000,9000,9010,9130,9599/
      RSUD=.TRUE.
c      RSUD=.FALSE.
      JMIN=1
c      PRINT *,' NUMBER OF RUNS?'
c      READ *, JMAX      
      JMAX=20
c      PRINT *,' RANDOM NUMBER SEEDS?'
c      READ *, NRS
      NRS(1)=54321
      NRS(2)=98765
      OPEN(UNIT=56,FILE='hermass.txt')
      OPEN(UNIT=57,FILE='hermass.top')
      MAXEV=10000
      IPROC=13010
      JPRO=MOD(IPROC,10000)
      PBEAM1=7000.
      PBEAM2=7000.
      PART1='P'
      PART2='P'
C---INITIALISE OTHER COMMON BLOCKS
      CALL HWIGIN
C---USER CAN RESET PARAMETERS AT
C   THIS POINT, OTHERWISE VALUES
C   SET IN HWIGIN WILL BE USED.
      EFFMIN=1D-5
      LWSUD=0
CJEM      LWSUD=77
      LRSUD=77
CJEM      LRSUD=0
      LWDEC=0
CJEM      LWDEC=88
      LRDEC=88
CJEM      LRDEC=0
      DO 999 J=JMIN,JMAX
         NRN(1)=NRS(1)
         NRN(2)=NRS(2)
         REWIND(LRDEC)
         IF (RSUD) THEN
            RSUD=.FALSE.
         ELSE
C---LRSUD.lt.0 suppresses reading and computing sud
            LRSUD=-1
         ENDIF
         MAXPR=0
         PRVTX=.FALSE.
         MAXER=100
         TLOUT=15.
         VPCUT=1.
         ZMXISR=ZERO
         PRESPL=.FALSE.
         TMNISR=1D-2
C---READ IN SUSY INPUT FILE, IN THIS CASE SNOWMASS POINT 1a
         OPEN(UNIT=LRSUSY,FORM='FORMATTED',STATUS='UNKNOWN',
     &        FILE='susyhit_softsusy_ISAWIG-test.out')
C     &        FILE='sps_pt1a.1200.in')
         CALL HWISSP
         CLOSE(UNIT=LRSUSY)
C---COMPUTE PARAMETER-DEPENDENT CONSTANTS
         CALL HWUINC
C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
         CALL HWUSTA('PI0     ')
C---USER'S INITIAL CALCULATIONS
         CALL HWABEG
C---INITIALISE ELEMENTARY PROCESS
         CALL HWEINI
C---LOOP OVER EVENTS
         DO 100 N=1,MAXEV
C---INITIALISE EVENT
            CALL HWUINE
C---GENERATE HARD SUBPROCESS
            CALL HWEPRO
C---GENERATE PARTON CASCADES
            CALL HWBGEN
C---DO HEAVY OBJECT DECAYS
            CALL HWDHOB
C---DO CLUSTER HADRONIZATION
c     CALL HWCFOR
C---DO CLUSTER DECAY
c     CALL HWCDEC
C---DO UNSTABLE PARTICLE DECAYS
c     CALL HWDHAD
C---DO HEAVY FLAVOUR DECAYS
c     CALL HWDHVY
C---ADD SOFT UNDERLYING EVENT IF NEEDED
c     CALL HWMEVT
C---FINALISE EVENT
            CALL HWUFNE
            NRS(1)=NRN(1)
            NRS(2)=NRN(2)
C---USER'S EVENT ANALYSIS
            CALL HWANAL
C---CHECK TO SEE IF ENOUGH GOOD EVENTS GENERATED
            IF (IERROR.EQ.191) GOTO 200
 100     CONTINUE
C---TERMINATE ELEMENTARY PROCESS
 200     CALL HWEFIN
C---USER'S TERMINAL CALCULATIONS
         CALL HWAEND
 999  CONTINUE
      END
CDECK  ID>, HWABEG. 
*CMZ :-        -26/04/91  10.18.56  by  Bryan Webber
*-- Author :    Bryan Webber
C----------------------------------------------------------------------
      SUBROUTINE HWABEG
C     USER'S ROUTINE FOR INITIALIZATION
C----------------------------------------------------------------------
      DOUBLE PRECISION MDA,MDC,MDM,MLSP
      INTEGER MDMAX,MDEVT,MDCOR
      PARAMETER (MDMAX=25)
      COMMON/MDCOM/MDA(8,8,2,MDMAX),MDC(8,2,MDMAX),MDM(8,MDMAX),
     & MLSP,MDEVT,MDCOR
      MDEVT=0
      END
CDECK  ID>, HWAEND.
*CMZ :-        -26/04/91  10.18.56  by  Bryan Webber
*-- Author :    Bryan Webber
C----------------------------------------------------------------------
      SUBROUTINE HWAEND
C     USER'S ROUTINE FOR TERMINAL CALCULATIONS, HISTOGRAM OUTPUT, ETC
C----------------------------------------------------------------------
      INTEGER I
      EXTERNAL MDFCN
C--FIT MASSES USING MINUIT
      OPEN(UNIT=55,FILE='hermass.in')
      CALL MINTIO(55,56,7)
      CALL MINUIT(MDFCN,0)
      REWIND(UNIT=55)
      END
CDECK  ID>, HWANAL. 
*CMZ :-        -26/04/91  11.11.54  by  Bryan Webber
*-- Author :    Bryan Webber
C----------------------------------------------------------------------
      SUBROUTINE HWANAL
C     USER'S ROUTINE TO ANALYSE DATA FROM EVENT
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INTEGER IHEP,I,IDP,IP1,IP2,ISQ,IC,IDCH(8)
      DOUBLE PRECISION PCH(5,8),PTM(2),MCH(8),PST(5),HWRGAU,RESLN,RESCA
      PARAMETER (RESLN=0.05D0)
      IF (IERROR.NE.0) RETURN
      ISQ=0
      DO IHEP=1,NHEP
         IF (ISTHEP(IHEP).EQ.155) THEN
            IDP=IDHEP(IHEP)
C            IF (ABS(MOD(IDP,1000000)).LT.6) THEN
C--FOUND SQUARK (NOT STOP)
            IF (ABS(IDP).GT.1000000.AND.ABS(IDP).LT.1000005) THEN
C--FOUND LEFT SQUARK (NOT ST OR SB) N.B. RIGHT SQUARKS HAVE SMALL B.R.
               ISQ=ISQ+1
               IF (ISQ.GT.2) THEN
                  PRINT *,' FOUND 3RD SQUARK: EVENT,IHEP =',NEVHEP,IHEP
                  CALL HWUEPR
                  RETURN
               ENDIF
               IC=4*ISQ-3
               IP1=JDAHEP(1,IHEP)
               IP2=JDAHEP(2,IHEP)
               IP2=JDAHEP(1,IP2)
               MCH(IC)=PHEP(5,IHEP)
               IDCH(IC)=IDHW(IP2)
               CALL HWVEQU(5,PHEP(1,IP2),PCH(1,IC))
               IC=IC+1
               IP1=JDAHEP(1,IP1)
               IP1=JDAHEP(1,IP1)
               MCH(IC)=PHEP(5,IP1)
               IP1=JDAHEP(1,IP1)
               IP2=JDAHEP(2,IP1)
C--NEAR LEPTON
               IDCH(IC)=IDHW(IP2)
               CALL HWVEQU(5,PHEP(1,IP2),PCH(1,IC))
               IC=IC+1
               IP1=JDAHEP(1,IP1)
               IP1=JDAHEP(1,IP1)
               IP2=JDAHEP(2,IP1)
               MCH(IC)=PHEP(5,IP1)
C--FAR LEPTON
               IDCH(IC)=IDHW(IP2)
               CALL HWVEQU(5,PHEP(1,IP2),PCH(1,IC))
               IC=IC+1
C--LSP
               IP2=JDAHEP(1,IP1)
               MCH(IC)=PHEP(5,IP2)
               IDCH(IC)=IDHW(IP2)
               CALL HWVEQU(5,PHEP(1,IP2),PCH(1,IC))
            ENDIF
         ENDIF
      ENDDO
      IF (ISQ.NE.2) RETURN
C--KEEP ONLY UNLIKE-FLAVOUR LEPTONS
      IF (IDCH(2).EQ.IDCH(6).OR.IDCH(2).EQ.IDCH(7)) RETURN
C--NEGLECT JET MASSES: RESCALE ENERGY
C      PCH(4,1)=SQRT(PCH(1,1)**2+PCH(2,1)**2+PCH(3,1)**2)
C      PCH(4,5)=SQRT(PCH(1,5)**2+PCH(2,5)**2+PCH(3,5)**2)
C--OR RESCALE MOMENTUM
C      RESCA= PCH(4,1)/SQRT(PCH(1,1)**2+PCH(2,1)**2+PCH(3,1)**2)
C      DO I=1,3
C         PCH(I,1)=PCH(I,1)*RESCA
C      ENDDO
C      RESCA= PCH(4,5)/SQRT(PCH(1,5)**2+PCH(2,5)**2+PCH(3,5)**2)
C      DO I=1,3
C         PCH(I,5)=PCH(I,5)*RESCA
C      ENDDO
C      PCH(5,1)=0D0
C      PCH(5,5)=0D0
C--RESCALE MOMENTA WITH RMS SMEARING RESLN
      IF (RESLN.GT.0D0) THEN
         DO IC=1,8
            RESCA=HWRGAU(IC,ONE,RESLN)
            DO I=1,3
               PCH(I,IC)=PCH(I,IC)*RESCA
            ENDDO
            PCH(4,IC)=SQRT(PCH(1,IC)**2+PCH(2,IC)**2+PCH(3,IC)**2
     &                    +PCH(5,IC)**2)
         ENDDO
      ENDIF
      DO I=1,2
         PTM(I)=PCH(I,4)+PCH(I,8)
      ENDDO
      CALL MDSET(PTM,PCH,MCH,1)
      IF (IERROR.NE.0) RETURN
C      IF (NEVHEP.LE.50) THEN
C         PRINT *
C         PRINT *,' EVENT ',NEVHEP
C         DO IC=1,8
C            PRINT '(I3,X,F9.2,X,A8,X,5F9.2)',IC,MCH(IC),
C     &           RNAME(IDCH(IC)),(PCH(I,IC),I=1,5)
C         ENDDO
C      ENDIF
C--CHAINS WITH 1 AND 5 TRANSPOSED
      CALL HWVEQU(5,PCH,PST)
      CALL HWVEQU(5,PCH(1,5),PCH)
      CALL HWVEQU(5,PST,PCH(1,5))
      CALL MDSET(PTM,PCH,MCH,2)
      END
C----------------------------------------------------------------------
      SUBROUTINE MDFCN(NVAR,GRAD,CHISQ,MVAR,IFLAG,FUTIL)
C--FUNCTION CALLED BY MINUIT TO CALCULATE CHISQ
      IMPLICIT NONE
      DOUBLE PRECISION GRAD,CHISQ,MVAR,FUTIL,MDFIT,MLO,MHI,MPT
      INTEGER NVAR,IFLAG
      DIMENSION GRAD(NVAR),MVAR(NVAR)
      DOUBLE PRECISION MDA,MDC,MDM,MLSP
      INTEGER MDMAX,MDEVT,MDCOR
      PARAMETER (MDMAX=25)
      COMMON/MDCOM/MDA(8,8,2,MDMAX),MDC(8,2,MDMAX),MDM(8,MDMAX),
     & MLSP,MDEVT,MDCOR
      IF (IFLAG.EQ.1) THEN
         PRINT *
         PRINT *,' NUMBER OF GOOD EVENTS = ',MDEVT
         IF (NVAR.EQ.3) THEN
            PRINT *
            PRINT *,' INPUT LSP MASS?'
            READ *,MLSP
            IF (MLSP.LT.-1D3) STOP
         ENDIF
      ENDIF
      CHISQ=MDFIT(MVAR,NVAR,IFLAG)
      IF (IFLAG.EQ.3) THEN
C--PRINT OUT BEST FIT
         IF (NVAR.EQ.3) THEN
            PRINT *, MVAR,MLSP,CHISQ
C--Write  Msq  Mn2  Msl  Mlsp  Chisq  Nev  Ngood
            WRITE(57,'(5F10.2,2I5)') MVAR,MLSP,CHISQ,MDEVT,MDCOR
         ELSE
            PRINT *, MVAR,CHISQ
            WRITE(57,'(5F10.2,2I5)') MVAR,CHISQ,MDEVT,MDCOR
C--Making a plot of chisq vs Mvar(4)
c            WRITE(57,'(''('',5F10.2,2I5)') MVAR,CHISQ,MDEVT,MDCOR
c            MLO=80D0
c            MHI=110D0
c            DO MPT=MLO,MHI,0.25D0
c               MVAR(4)=MPT
c               CHISQ=MDFIT(MVAR,NVAR,IFLAG)
c               WRITE(57,'(2E12.4)') MPT,CHISQ
c            ENDDO
c            WRITE(57,*) ' JOIN 1' 
         ENDIF
      ENDIF
      END
C----------------------------------------------------------------------
      FUNCTION MDFIT(MVAR,NVAR,IFLAG)
C--MASS DETERMINATION FOR PAIR OF 3-STEP CHAINS
C     Z->1+Y, Y->2+X, X->3+4,    Z'->5+Y', Y'->6+X', X'->7+8
C  WHERE 4 AND 8 ARE INVISIBLE (N & N')
C  INPUT MCH = (MZ,MY,MX,MN,MZ',MY',MN') HYPOTHESIS
C        NCH = NUMBER OF MASSES TO VARY
C OUTPUT MASDET = CHISQ FOR THIS HYPOTHESIS
C----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION MDFIT,MVAR,MCH(8),MMCH(8),MMTR(8),PINV(5,2),
     & MMINV,CHISQ,CHILO,CHIEVT,CHISUM
      INTEGER NVAR,IFLAG,I,MDLIM,JCH,ICH,JPRM,IPRM,NPRM(4,2)
      DIMENSION MVAR(NVAR)
      DOUBLE PRECISION MDA,MDC,MDM,MLSP
      INTEGER MDMAX,MDEVT,MDCOR
      PARAMETER (MDMAX=25)
      COMMON/MDCOM/MDA(8,8,2,MDMAX),MDC(8,2,MDMAX),MDM(8,MDMAX),
     & MLSP,MDEVT,MDCOR
      IF (NVAR.EQ.8) THEN
         DO I=1,8
            MCH(I)=MVAR(I)
         ENDDO
      ELSEIF (NVAR.EQ.4) THEN
         DO I=1,4
            MCH(I)=MVAR(I)
            MCH(I+4)=MVAR(I)
         ENDDO
      ELSEIF (NVAR.EQ.3) THEN
         DO I=1,3
            MCH(I)=MVAR(I)
            MCH(I+4)=MVAR(I)
         ENDDO
         MCH(4)=MLSP
         MCH(8)=MLSP
      ELSE
         PRINT *,' MDFIT: BAD NVAR = ',NVAR
         STOP
      ENDIF
      DO I=1,8
         MMCH(I)=MCH(I)**2
      ENDDO
      CHISQ=0D0
      IF (MCH(1).LT.0D0) CHISQ=CHISQ+MMCH(1)**2
      IF (MCH(2).LT.0D0) CHISQ=CHISQ+MMCH(2)**2
      IF (MCH(3).LT.0D0) CHISQ=CHISQ+MMCH(3)**2
      IF (MCH(4).LT.0D0) CHISQ=CHISQ+MMCH(4)**2
      IF (MCH(5).LT.0D0) CHISQ=CHISQ+MMCH(5)**2
      IF (MCH(6).LT.0D0) CHISQ=CHISQ+MMCH(6)**2
      IF (MCH(7).LT.0D0) CHISQ=CHISQ+MMCH(7)**2
      IF (MCH(8).LT.0D0) CHISQ=CHISQ+MMCH(8)**2
      IF (MCH(1).LT.MCH(2)) CHISQ=CHISQ+(MMCH(1)-MMCH(2))**2
      IF (MCH(2).LT.MCH(3)) CHISQ=CHISQ+(MMCH(2)-MMCH(3))**2
      IF (MCH(3).LT.MCH(4)) CHISQ=CHISQ+(MMCH(3)-MMCH(4))**2
      IF (MCH(5).LT.MCH(6)) CHISQ=CHISQ+(MMCH(5)-MMCH(6))**2
      IF (MCH(6).LT.MCH(7)) CHISQ=CHISQ+(MMCH(6)-MMCH(7))**2
      IF (MCH(7).LT.MCH(8)) CHISQ=CHISQ+(MMCH(7)-MMCH(8))**2
      CHILO=CHISQ
      CHISUM=0D0
      MDLIM=MDEVT
      MDEVT=0
      DO JCH=1,2
         DO JPRM=1,4
            NPRM(JPRM,JCH)=0
         ENDDO
      ENDDO
      DO WHILE (MDEVT.LT.MDLIM)
         MDEVT=MDEVT+1
C         IF (IFLAG.EQ.3) THEN
C            PRINT *
C            PRINT *,' EVENT  ',MDEVT,'  FIT & TRUTH'
C            PRINT '(8F8.2)',MCH,(MDM(I,MDEVT),I=1,8)
C            DO I=1,8
C               MMTR(I)=MDM(I,MDEVT)**2
C            ENDDO
C            CALL MDGET(MMTR,1,1,PINV)
C            PRINT '(5F8.2)',PINV
C         ENDIF
         CHIEVT=1D50
         DO JCH=1,2
            DO JPRM=1,4
               CALL MDGET(MMCH,JCH,JPRM,PINV)
               MMINV=SIGN(PINV(5,1)**2,PINV(5,1))
               CHISQ=CHILO+(MMCH(4)-MMINV)**2
               MMINV=SIGN(PINV(5,2)**2,PINV(5,2))
               CHISQ=CHISQ+(MMCH(8)-MMINV)**2
               IF (CHISQ.LT.CHIEVT) THEN
                  CHIEVT=CHISQ
                  IPRM=JPRM
                  ICH=JCH
               ENDIF
C               IF (IFLAG.EQ.3) THEN
C                  PRINT *,' CHAIN,PERM   ', JCH,JPRM,'  CHISQ = ',CHISQ
C                  PRINT '(5F8.2)',PINV
C               ENDIF
            ENDDO
         ENDDO
         IF (IFLAG.EQ.3) THEN
C            PRINT *,' SELECT ', ICH,IPRM,'  CHISQ = ',CHIEVT
            NPRM(IPRM,ICH)=NPRM(IPRM,ICH)+1
         ENDIF
         CHISUM=CHISUM+CHIEVT
      ENDDO
      IF (IFLAG.EQ.3) THEN
c         PRINT *,' NPRM = ',NPRM
         MDCOR=NPRM(1,1)
      ENDIF
      MDFIT=CHISUM*1D-8
      END
C----------------------------------------------------------------------
      SUBROUTINE MDGET(MMCH,JCH,JPRM,PINV)
C--MASS DETERMINATION FOR PAIR OF 3-STEP CHAINS
C     Z->1+Y, Y->2+X, X->3+4,    Z'->5+Y', Y'->6+X', X'->7+8
C  WHERE 4 AND 8 ARE INVISIBLE (N & N')
C  INPUT MMCH = (MZ^2,MY^2,MX^2,MN^2,MZ'^2,MY'^2,MN'^2) HYPOTHESIS
C        JPRM = PERMUTATION (WITHIN CHAINS)
C           A = INVERSE A-MATRIX FROM MDSET
C           C = C-VECTOR FROM MDSET
C OUTPUT PINV(I,J) = INVISIBLE 5-MOMENTA (J=1,2)
C----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION MMCH(8),PINV(5,2),A(8,8),C(8),S(8)
      INTEGER I,J,JJ,K,KK,JCH,JPRM,PERM(8,4)
      DOUBLE PRECISION MDA,MDC,MDM,MLSP
      INTEGER MDMAX,MDEVT,MDCOR
      PARAMETER (MDMAX=25)
      COMMON/MDCOM/MDA(8,8,2,MDMAX),MDC(8,2,MDMAX),MDM(8,MDMAX),
     & MLSP,MDEVT,MDCOR
      DATA PERM/1,2,3,4,5,6,7,8, 1,3,2,4,5,6,7,8,
     &          1,2,3,4,5,7,6,8, 1,3,2,4,5,7,6,8/
      DO J=1,8
         JJ=PERM(J,JPRM)
         CALL HWVEQU(8,MDA(1,JJ,JCH,MDEVT),A(1,J))
      ENDDO
      CALL HWVEQU(8,MDC(1,JCH,MDEVT),C)
      DO I=1,7
         S(I)=C(I)-MMCH(I)+MMCH(I+1)
      ENDDO
      S(4)=C(4)
      S(8)=C(8)
C--SOLVE FOR THE INVISIBLE MOMENTA
      KK=0
      DO K=1,2
         DO I=1,4
            PINV(I,K)=0D0
            DO J=1,8
               PINV(I,K)=PINV(I,K)+A(I+KK,J)*S(J)
            ENDDO
         ENDDO
         CALL HWUMAS(PINV(1,K))
         KK=4
      ENDDO
      END
C----------------------------------------------------------------------
      SUBROUTINE MDSET(PTM,PCH,MCH,JCH)
C--MASS DETERMINATION FOR PAIR OF 3-STEP CHAINS
C     Z->1+Y, Y->2+X, X->3+4,    Z'->5+Y', Y'->6+X', X'->7+8
C  WHERE 4 AND 8 ARE INVISIBLE (N & N')
C        PTM(I) = MISSING PTX,PTY (I=1,2)
C        PCH(I,J) = VISIBLE 5-MOMENTA (J=1,2,3,5,6,7)
C OUTPUT INVERSE MATRIX A AND VECTOR C FOR SOLVING FOR INVISIBLES
C----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION PTM(2),PCH(5,8),MCH(8),A(8,8),C(8),HWULDO
      INTEGER JCH,I,J,IR(8),IFAIL
      DOUBLE PRECISION MDA,MDC,MDM,MLSP
      INTEGER MDMAX,MDEVT,MDCOR
      PARAMETER (MDMAX=25)
      COMMON/MDCOM/MDA(8,8,2,MDMAX),MDC(8,2,MDMAX),MDM(8,MDMAX),
     & MLSP,MDEVT,MDCOR
      IF (MDEVT.GE.MDMAX) THEN
         MDEVT=MDMAX
         CALL HWWARN('MDSET',191)
         RETURN
      ENDIF
C--FILL MATRIX A
      DO I=1,3
         DO J=1,4
            A(I,J)=2D0*PCH(J,I)
            A(I,J+4)=0D0
         ENDDO
         A(I,4)=-A(I,4)
      ENDDO
      DO I=5,7
         DO J=1,4
            A(I,J)=0D0
            A(I,J+4)=2D0*PCH(J,I)
         ENDDO
         A(I,8)=-A(I,8)
      ENDDO
      DO J=1,8
         A(4,J)=0D0
         A(8,J)=0D0
      ENDDO
      A(4,1)=1D0
      A(4,5)=1D0
      A(8,2)=1D0
      A(8,6)=1D0
C--CERNLIB inversion routine. NB it destroys input matrix
      CALL DINV(8,A,8,IR,IFAIL)
      IF (IFAIL.NE.0) THEN
         PRINT *
         PRINT *,' MDSET: INVERSION FAILED, IFAIL = ',IFAIL
         IF (JCH.NE.1) MDEVT=MDEVT-1
         RETURN
      ENDIF
C--FILL THE COMPONENTS OF C
      C(1)=PCH(5,1)**2+2D0*HWULDO(PCH(1,1),PCH(1,2))
     &     +2D0*HWULDO(PCH(1,1),PCH(1,3))
      C(2)=PCH(5,2)**2+2D0*HWULDO(PCH(1,2),PCH(1,3))
      C(3)=PCH(5,3)**2
      C(4)=PTM(1)
      C(5)=PCH(5,5)**2+2D0*HWULDO(PCH(1,5),PCH(1,6))
     &     +2D0*HWULDO(PCH(1,5),PCH(1,7))
      C(6)=PCH(5,6)**2+2D0*HWULDO(PCH(1,6),PCH(1,7))
      C(7)=PCH(5,7)**2
      C(8)=PTM(2)
      IF (JCH.EQ.1) MDEVT=MDEVT+1
      CALL HWVEQU(64,A,MDA(1,1,JCH,MDEVT))
      CALL HWVEQU(8,C,MDC(1,JCH,MDEVT))
      CALL HWVEQU(8,MCH,MDM(1,MDEVT))
      END
*
* $Id: dinv.F,v 1.1.1.1 1996/02/15 17:48:49 mclareni Exp $
*
* $Log: dinv.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:49  mclareni
* Kernlib
*
*
C#include "kernnum/pilot.h"
      SUBROUTINE DINV(N,A,IDIM,R,IFAIL)
      INTEGER R(N)
      REAL T1,T2,T3
      DOUBLE PRECISION A(IDIM,N),DET,TEMP,S,
     $                 C11,C12,C13,C21,C22,C23,C31,C32,C33
      CHARACTER*6 NAME
      DATA NAME/'DINV'/,KPRNT/0/
C
C     ******************************************************************
C
C     REPLACES A BY ITS INVERSE.
C
C     (PARAMETERS AS FOR DEQINV.)
C
C     CALLS ... DFACT, DFINV, F010PR, ABEND.
C
C     ******************************************************************
C
C  TEST FOR PARAMETER ERRORS.
C
      IF((N.LT.1).OR.(N.GT.IDIM)) GO TO 7
C
C  TEST FOR N.LE.3.
C
      IF(N.GT.3) GO TO 6
      IFAIL=0
      IF(N.LT.3) GO TO 4
C
C  N=3 CASE.
C
C     COMPUTE COFACTORS.
      C11=A(2,2)*A(3,3)-A(2,3)*A(3,2)
      C12=A(2,3)*A(3,1)-A(2,1)*A(3,3)
      C13=A(2,1)*A(3,2)-A(2,2)*A(3,1)
      C21=A(3,2)*A(1,3)-A(3,3)*A(1,2)
      C22=A(3,3)*A(1,1)-A(3,1)*A(1,3)
      C23=A(3,1)*A(1,2)-A(3,2)*A(1,1)
      C31=A(1,2)*A(2,3)-A(1,3)*A(2,2)
      C32=A(1,3)*A(2,1)-A(1,1)*A(2,3)
      C33=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      T1=ABS(SNGL(A(1,1)))
      T2=ABS(SNGL(A(2,1)))
      T3=ABS(SNGL(A(3,1)))
C
C     (SET TEMP=PIVOT AND DET=PIVOT*DET.)
      IF(T1.GE.T2) GO TO 1
         IF(T3.GE.T2) GO TO 2
C        (PIVOT IS A21)
            TEMP=A(2,1)
            DET=C13*C32-C12*C33
            GO TO 3
    1 IF(T3.GE.T1) GO TO 2
C     (PIVOT IS A11)
         TEMP=A(1,1)
         DET=C22*C33-C23*C32
         GO TO 3
C     (PIVOT IS A31)
    2    TEMP=A(3,1)
         DET=C23*C12-C22*C13
C
C     SET ELEMENTS OF INVERSE IN A.
    3 IF(DET.EQ.0D0) GO TO 8
      S=TEMP/DET
      A(1,1)=S*C11
      A(1,2)=S*C21
      A(1,3)=S*C31
      A(2,1)=S*C12
      A(2,2)=S*C22
      A(2,3)=S*C32
      A(3,1)=S*C13
      A(3,2)=S*C23
      A(3,3)=S*C33
      RETURN
C
    4 IF(N.LT.2) GO TO 5
C
C  N=2 CASE BY CRAMERS RULE.
C
      DET=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      IF(DET.EQ.0D0) GO TO 8
      S=1D0/DET
      C11   =S*A(2,2)
      A(1,2)=-S*A(1,2)
      A(2,1)=-S*A(2,1)
      A(2,2)=S*A(1,1)
      A(1,1)=C11
      RETURN
C
C  N=1 CASE.
C
    5 IF(A(1,1).EQ.0D0) GO TO 8
      A(1,1)=1D0/A(1,1)
      RETURN
C
C  N.GT.3 CASES.  FACTORIZE MATRIX AND INVERT.
C
    6 CALL DFACT(N,A,IDIM,R,IFAIL,DET,JFAIL)
      IF(IFAIL.NE.0) RETURN
      CALL DFINV(N,A,IDIM,R)
      RETURN
C
C  ERROR EXITS.
C
    7 IFAIL=+1
c      CALL F010PR(NAME,N,IDIM,K,KPRNT)
      RETURN
C
    8 IFAIL=-1
      RETURN
C
      END
*
* $Id: dfact.F,v 1.1.1.1 1996/02/15 17:49:03 mclareni Exp $
*
* $Log: dfact.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*
C#include "kernnum/pilot.h"
          SUBROUTINE          DFACT(N,A,IDIM,IR,IFAIL,DET,JFAIL)
          INTEGER             IR(*),    IPAIRF
          DOUBLE PRECISION    A(IDIM,*),DET,      ZERO,     ONE,X,Y,TF
          REAL                G1,       G2
          REAL                PIVOTF,   P,        Q,        SIZEF,  T
          DOUBLE PRECISION    S11, S12, DOTF
          CHARACTER*6         HNAME
          IPAIRF(J,K)  =  J*2**12 + K
          PIVOTF(X)    =  ABS(SNGL(X))
          SIZEF(X)     =  ABS(SNGL(X))
          DOTF(X,Y,S11)  =  X * Y + S11
          DATA      G1, G2              /  1.E-37,  1.E37  /
C#endif
C#if defined(CERNLIB_NUME38)
c          DATA      G1, G2              /  1.E-19,  1.E19  /
C#endif
C#if defined(CERNLIB_NUME999)
c          DATA      ?????  G1, G2 NOT DEFINED  ?????
C#endif
          DATA      HNAME               /  ' DFACT'  /
          DATA      ZERO, ONE           /  0.D0, 1.D0  /
          DATA      NORMAL, IMPOSS      /  0, -1  /
          DATA      JRANGE, JOVER, JUNDER  /  0, +1, -1  /
C#include "fact.inc"
*
* $Id: fact.inc,v 1.1.1.1 1996/02/15 17:49:03 mclareni Exp $
*
* $Log: fact.inc,v $
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*
*
* fact.inc
*
          IF(IDIM .GE. N  .AND.  N .GT. 0)  GOTO 110
C             CALL TMPRNT(HNAME,N,IDIM,0)
          IFAIL=101
             RETURN
 110      IFAIL  =  NORMAL
          JFAIL  =  JRANGE
          NXCH   =  0
          DET    =  ONE
          DO 144    J  =  1, N
 120         K  =  J
             P  =  PIVOTF(A(J,J))
             IF(J .EQ. N)  GOTO 122
             JP1  =  J+1
             DO 121    I  =  JP1, N
                Q  =  PIVOTF(A(I,J))
                IF(Q .LE. P)  GOTO 121
                   K  =  I
                   P  =  Q
 121            CONTINUE
             IF(K .NE. J)  GOTO 123
 122         IF(P .GT. 0.)  GOTO 130
                DET    =  ZERO
                IFAIL  =  IMPOSS
                JFAIL  =  JRANGE
                RETURN
 123         DO 124    L  =  1, N
                TF      =  A(J,L)
                A(J,L)  =  A(K,L)
                A(K,L)  =  TF
 124            CONTINUE
             NXCH      =  NXCH + 1
             IR(NXCH)  =  IPAIRF(J,K)
 130         DET     =  DET * A(J,J)
             A(J,J)  =  ONE / A(J,J)
             T  =  SIZEF(DET)
             IF(T .LT. G1)  THEN
                DET    =  ZERO
                IF(JFAIL .EQ. JRANGE)  JFAIL  =  JUNDER
             ELSEIF(T .GT. G2)  THEN
                DET    =  ONE
                IF(JFAIL .EQ. JRANGE)  JFAIL  =  JOVER
             ENDIF
             IF(J .EQ. N)  GOTO 144
             JM1  =  J-1
             JP1  =  J+1
             DO 143   K  =  JP1, N
                S11  =  -A(J,K)
                S12  =  -A(K,J+1)
                IF(J .EQ. 1)  GOTO 142
                DO 141  I  =  1, JM1
                   S11  =  DOTF(A(I,K),A(J,I),S11)
                   S12  =  DOTF(A(I,J+1),A(K,I),S12)
 141               CONTINUE
 142            A(J,K)    =  -S11 * A(J,J)
                A(K,J+1)  =  -DOTF(A(J,J+1),A(K,J),S12)
 143            CONTINUE
 144         CONTINUE
 150      IF(MOD(NXCH,2) .NE. 0)  DET  =  -DET
          IF(JFAIL .NE. JRANGE)   DET  =  ZERO
          IR(N)  =  NXCH
          RETURN
          END
*
* $Id: dfinv.F,v 1.1.1.1 1996/02/15 17:49:03 mclareni Exp $
*
* $Log: dfinv.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*
C#include "kernnum/pilot.h"
          SUBROUTINE          DFINV(N,A,IDIM,IR)
          INTEGER             IR(*)
          DOUBLE PRECISION    A(IDIM,*),ZERO,     X, Y, TI
          DOUBLE PRECISION    S31, S32, S33, S34, DOTF
          CHARACTER*6         HNAME
          DATA      HNAME               /  ' DFINV'  /
          DOTF(X,Y,S31)  =  X*Y + S31
          DATA      ZERO      /  0.D0  /
C#include "finv.inc"
*
* $Id: finv.inc,v 1.1.1.1 1996/02/15 17:49:03 mclareni Exp $
*
* $Log: finv.inc,v $
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*
*
* finv.inc
*
          IF(IDIM .GE. N  .AND.  N .GT. 0)  GOTO 310
C             CALL TMPRNT(HNAME,N,IDIM,0)
          IFAIL=102
             RETURN
 310      IF(N .EQ. 1)  RETURN
          A(2,1)  =  -A(2,2) * DOTF(A(1,1),A(2,1),ZERO)
          A(1,2)  =  -A(1,2)
          IF(N .EQ. 2)  GOTO 330
          DO 314    I  =  3, N
             IM2  =  I-2
             DO 312 J  =  1, IM2
                S31  =  ZERO
                S32  =  A(J,I)
                DO 311  K  =  J, IM2
                   S31  =  DOTF(A(K,J),A(I,K),S31)
                   S32  =  DOTF(A(J,K+1),A(K+1,I),S32)
 311               CONTINUE
                A(I,J)  =  -A(I,I) * DOTF(A(I-1,J),A(I,I-1),S31)
                A(J,I)  =  -S32
 312            CONTINUE
             A(I,I-1)  =  -A(I,I) * DOTF(A(I-1,I-1),A(I,I-1),ZERO)
             A(I-1,I)  =  -A(I-1,I)
 314         CONTINUE
 330      NM1  =  N-1
          DO 335   I  =  1, NM1
             NMI  =  N-I
             DO 332   J  =  1, I
                S33  =  A(I,J)
                DO 331   K  =  1, NMI
                   S33  =  DOTF(A(I+K,J),A(I,I+K),S33)
 331               CONTINUE
                A(I,J)  =  S33
 332            CONTINUE
             DO 334   J  =  1, NMI
                S34  =  ZERO
                DO 333   K  =  J, NMI
                   S34  =  DOTF(A(I+K,I+J),A(I,I+K),S34)
 333               CONTINUE
                A(I,I+J)  =  S34
 334            CONTINUE
 335         CONTINUE
          NXCH  =  IR(N)
          IF(NXCH .EQ. 0)  RETURN
            DO 342 M  =  1, NXCH
             K   =  NXCH - M+1
             IJ  =  IR(K)
             I   =  IJ / 4096
             J   =  MOD(IJ,4096)
             DO 341  K  =  1, N
                TI      =  A(K,I)
                A(K,I)  =  A(K,J)
                A(K,J)  =  TI
 341            CONTINUE
 342         CONTINUE
          RETURN
          END
C----------------------------------------------------------------------
      SUBROUTINE UPEVNT
      END
C----------------------------------------------------------------------
      SUBROUTINE UPINIT
      END
C----------------------------------------------------------------------
      SUBROUTINE PDFSET(PARM,VAL)
C----------------------------------------------------------------------
C     DUMMY SUBROUTINE: DELETE AND SET MODPDF(I)
C     IN MAIN PROGRAM IF YOU USE PDFLIB CERN-LIBRARY
C     PACKAGE FOR NUCLEON STRUCTURE FUNCTIONS
C----------------------------------------------------------------------
      DOUBLE PRECISION VAL(20)
      CHARACTER*20 PARM(20)
      WRITE (6,10)
   10 FORMAT(/10X,'PDFSET CALLED BUT NOT LINKED')
      STOP
      END
C-----------------------------------------------------------------------
      SUBROUTINE STRUCTM(X,QSCA,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
C-----------------------------------------------------------------------
C     DUMMY SUBROUTINE: DELETE IF YOU USE PDFLIB CERN-LIBRARY
C     PACKAGE FOR NUCLEON STRUCTURE FUNCTIONS
C-----------------------------------------------------------------------
      DOUBLE PRECISION X,QSCA,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU
      WRITE (6,10)
  10  FORMAT(/10X,'STRUCTM CALLED BUT NOT LINKED')
      STOP
      END
