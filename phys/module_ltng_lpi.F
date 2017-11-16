MODULE MODULE_LTNG_LPI
!YAIR, Y., B. LYNN, C. PRICE, V. KOTRONI, K. LAGOUVARDOS, E. MORIN,
!A. MAGNAI, AND M. DEL CARMEN LLASAT (2010), PREDICTING THE POTENTIAL FOR
!LIGHTNING ACTIVITY IN MEDITERRANEAN STORMS BASED ON THE WEATHER
!RESEARCH AND FORECASTING (WRF) MODEL DYNAMIC AND MICROPHYSICAL
!FIELDS, J. GEOPHYS. RES., 115, D04205, DOI:10.1029/2008JD010868.
! HOWEVER, WE DON'T CHECK FOR COLLAPSING CELL (SO AS NOT TO REQUIRE USE OF HALO).
! THIS MEANS THAT LPI IS ALSO CALCULATED IN CELLS THAT ARE NO LONGER (ON AVERAGE) GROWING
! FOR A "COMPLETE" LIGHTNING FORECAST SCHEME, PLEASE SEE:
!HTTP://JOURNALS.AMETSOC.ORG/DOI/ABS/10.1175/WAF-D-11-00144.1
!(PREDICTING CLOUD-TO-GROUND AND INTRACLOUD LIGHTNING IN WEATHER FORECAST MODELS)

CONTAINS
!===================================================================
!
  SUBROUTINE CALCLPI(QV,QC, QR, QI, QS, QG, QH                            &
                 ,W,Z,DZ8W,PI_PHY,TH_PHY,P_PHY,RHO_PHY                    &
                 ,LPI&
                 ,IDS,IDE, JDS,JDE, KDS,KDE                        &
                 ,IMS,IME, JMS,JME, KMS,KME                        &
                 ,ITS,ITE, JTS,JTE, KTS,KTE                        &
                                                                   )
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!
!
  INTEGER,      INTENT(IN   )    ::   IDS,IDE, JDS,JDE, KDS,KDE , &
                                      IMS,IME, JMS,JME, KMS,KME , &
                                      ITS,ITE, JTS,JTE, KTS,KTE
  REAL, DIMENSION( IMS:IME , KMS:KME , JMS:JME ),                 &
        INTENT(IN) ::                                          &
                                                              QV, &
                                                              QC, &
                                                              QI, &
                                                              QR, &
                                                              QS, &
                                                              QG,QH

      REAL, DIMENSION( IMS:IME, KMS:KME, JMS:JME ),                  &
         INTENT(IN ) ::  W, Z
      REAL, INTENT(IN),     DIMENSION(IMS:IME, KMS:KME, JMS:JME)::      &
     &                      DZ8W,PI_PHY,P_PHY,RHO_PHY
      REAL, INTENT(IN),  DIMENSION(IMS:IME, KMS:KME, JMS:JME)::      &
     &                      TH_PHY
      REAL, INTENT(INOUT),  DIMENSION(IMS:IME,JMS:JME)::      &
     &                      LPI




      REAL, DIMENSION(KMS:KME)::    TEMPK,RH
      REAL, DIMENSION(KMS:KME):: QV1D,P1D,RHO1D,QTI1D
      REAL, DIMENSION(KMS:KME):: TEMP,QC1D,QL1D,QI1D,QS1D,QG1D,LPI1D
      REAL, DIMENSION(0:KME):: W1D,HEIGHT
      REAL, DIMENSION(KMS:KME):: E1D,HEIGHT_T,W1D_T
      REAL Z_FULL,QRS,TETEN,RELHUM,LOC,TD_850,TD_700,PC_DWPT
      INTEGER LEVEL
      REAL :: DT_LPI,T_BASE,T_TOP
      INTEGER I_COLLAPSE
      LOGICAL LOOK_T
      INTEGER I_START,I_END,J_START,J_END


  INTEGER ::               I,J,K
!-------------------------------------------------------------------
      DO J = JTS,JTE
      DO I = ITS,ITE
        Z_FULL=0.
        HEIGHT(0)=Z_FULL
        W1D(0)=W(I,1,J)
      DO K = KTS,KTE-1
          IF (K.LT.KTE-1)THEN
           W1D(K)=W(I,K+1,J)
          ELSE
           W1D(K)=0.
          END IF
          TEMP(K) = TH_PHY(I,K,J)*PI_PHY(I,K,J)-273.16
          TEMPK(K) = TH_PHY(I,K,J)*PI_PHY(I,K,J)
          QV1D(K)=QV(I,K,J)
          P1D(K)=P_PHY(I,K,J)
          RHO1D(K)=RHO_PHY(I,K,J)
          Z_FULL=Z_FULL+DZ8W(I,K,J)
          HEIGHT(K)=Z_FULL
          QC1D(K)=QC(I,K,J)
          QL1D(K)=QC(I,K,J)+QR(I,K,J)
          QI1D(K)=QI(I,K,J)
          QTI1D(K)=QI(I,K,J)+QS(I,K,J)+QG(I,K,J)+QH(I,K,J)
          QS1D(K)=QS(I,K,J)
!         QG1D(K)=QG(I,K,J)+QH(I,K,J)
! HAIL DOESN'T USUALLY CHARGE
          QG1D(K)=QG(I,K,J)
! FOR CONSERVATIVE ADVECTION MULTIPLY BY RHO1D AND DIVIDE BY IT BELOW
      ENDDO
      DO K = KTS,KTE-1
       HEIGHT_T(K)=0.5*(HEIGHT(K-1)+HEIGHT(K))
       W1D_T(K)=0.5*(W1D(K-1)+W1D(K))
      END DO
      T_BASE=-0
      T_TOP=-20
      CALL CALC_LPI(QL1D,QI1D,QS1D,QG1D,W1D,TEMP,HEIGHT,LPI(I,J),T_BASE,T_TOP,KME,KTE)
      END DO
      END DO
      RETURN
      END SUBROUTINE CALCLPI
      SUBROUTINE &
     &  CALC_LPI(QL3D,QI3D,QS3D,QG3D,W3D,T3D,HEIGHT,LPI,T_BASE,T_TOP,NK,NKE)
      IMPLICIT NONE
      INTEGER NK,NKE
      REAL T_BASE,T_TOP
      REAL QL3D(NK)
      REAL QG3D(NK)
      REAL QI3D(NK)
      REAL QS3D(NK)
      REAL W3D(0:NK)
      REAL T3D(NK)
      REAL HEIGHT(0:NK)
      REAL LPI
      REAL DEL_Z(NK)
      REAL W_AVE(NK)
      INTEGER IC,JC,ICNT,I,J,K,I_COLLAPSE
      REAL I_DIST,J_DIST,DEL_Z_TOT
      REAL TOP, BOT
      REAL NUM,DEN,AVE_Z
      REAL NUM_S,DEN_S
      REAL NUM_I,DEN_I
      REAL Q_ISG
      ICNT=0
      DO K=1,NKE
        TOP=HEIGHT(K)
        BOT=HEIGHT(K-1)
        DEL_Z(K)=TOP-BOT
        W_AVE(K)=0.5*(W3D(K)+W3D(K-1))
      END DO
!
!     CHECK FOR COLLAPSING CELL
! HERE, WE DON'T CHECK, SINCE IT REQUIRES A HALO.
      AVE_Z=0
      DEL_Z_TOT=0
      LPI=0
      DO K=1,NKE-1
       IF (T3D(K).LE.T_BASE.AND.T3D(K).GT.T_TOP)THEN ! SET TEMP RANGE
        
        DEN_I = QI3D(K)+QG3D(K)     
        DEN_S = QS3D(K)+QG3D(K)
        IF (QS3D(K).EQ.0.OR.QG3D(K).EQ.0.)THEN !CHECKS FOR ZEROES
         DEN_S=10000.
         NUM_S = 0.
        ELSE
         NUM_S = SQRT(QS3D(K)*QG3D(K))   
        END IF
        IF (QI3D(K).EQ.0.OR.QG3D(K).EQ.0.)THEN  ! CHECKS FOR ZEROES
         DEN_I=10000.
         NUM_I = 0.
        ELSE
         NUM_I = SQRT(QI3D(K)*QG3D(K))
        END IF
        Q_ISG = QG3D(K)*(NUM_I/DEN_I+NUM_S/DEN_S)  ! ICE "FRACT"-CONTENT

        IF (QL3D(K).EQ.0.OR.Q_ISG.EQ.0)THEN
          NUM=0
          DEN=10000.
        ELSE
         NUM = SQRT(QL3D(K)*Q_ISG)
         DEN = QL3D(K)+Q_ISG
        END IF
        DEL_Z_TOT=DEL_Z_TOT+DEL_Z(K)
        IF (NUM.GT.0)THEN
         AVE_Z=AVE_Z+DEL_Z(K)*(2.*NUM/DEN)*W_AVE(K)**2 ! LIGHTNING POTENTIAL INDEX J/UNIT-MASS
        END IF
       END IF
      END DO
!
      IF (DEL_Z_TOT.EQ.0)DEL_Z_TOT=100000
      LPI=AVE_Z/DEL_Z_TOT
       
!
      RETURN
      END SUBROUTINE CALC_LPI
  END MODULE MODULE_LTNG_LPI
