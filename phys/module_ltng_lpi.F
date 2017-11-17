MODULE MODULE_LTNG_LPI
!REFERENCE

      IMPLICIT NONE
      REAL, PARAMETER, PRIVATE:: T_0 = 273.15
      REAL, PARAMETER, PRIVATE:: T_BASE= 0
      REAL, PARAMETER, PRIVATE:: T_TOP=-20
CONTAINS
!===================================================================
!
      SUBROUTINE CALCLPI(                                         &
                         QV,QC,QR,QI,QS,QG,QH                     &
                        ,W,Z,DZ8W,PI_PHY,TH_PHY,P_PHY,RHO_PHY     &
                        ,LPI                                      &
                        ,IDS,IDE, JDS,JDE, KDS,KDE                &
                        ,IMS,IME, JMS,JME, KMS,KME                &
                        ,ITS,ITE, JTS,JTE, KTS,KTE                &
                                                                  )
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
!
!
      INTEGER, INTENT(IN) ::          IDS,IDE, JDS,JDE, KDS,KDE , &
                                      IMS,IME, JMS,JME, KMS,KME , &
                                      ITS,ITE, JTS,JTE, KTS,KTE
      REAL, DIMENSION( IMS:IME , KMS:KME , JMS:JME ),             &
            INTENT(IN) ::                                         &
                                                             QV , &
                                                             QC , &
                                                             QR , &
                                                             QI , &
                                                             QS , &
                                                             QG , &
                                                             QH

      REAL, DIMENSION( IMS:IME, KMS:KME, JMS:JME ),               &
            INTENT(IN ) ::  W,Z,DZ8W,PI_PHY,TH_PHY,P_PHY,RHO_PHY

      REAL, INTENT(INOUT), DIMENSION(IMS:IME,JMS:JME)::      LPI

! Define Local Data
      REAL, DIMENSION(KMS:KME):: TEMPC
      REAL, DIMENSION(KMS:KME):: QL1D,QI1D,QS1D,QG1D
      REAL, DIMENSION(0:KME):: W1D,HEIGHT
      REAL Z_FULL


      INTEGER ::               I,J,K
!-------------------------------------------------------------------
      DO J = JTS,JTE
      DO I = ITS,ITE
        Z_FULL=0.
        HEIGHT(0)=0.
        W1D(0)=W(I,1,J)
        DO K = KTS,KTE-1
          IF (K.LT.KTE-1)THEN
             W1D(K)=W(I,K+1,J)
          ELSE
             W1D(K)=0.
          END IF
          TEMPC(K) = TH_PHY(I,K,J)*PI_PHY(I,K,J)-T_0
          Z_FULL=Z_FULL+DZ8W(I,K,J)
          HEIGHT(K)=Z_FULL
          QL1D(K)=QC(I,K,J)+QR(I,K,J)
          QI1D(K)=QI(I,K,J)
          QS1D(K)=QS(I,K,J)
!         QG1D(K)=QG(I,K,J)+QH(I,K,J)
! HAIL DOESN'T USUALLY CHARGE
          QG1D(K)=QG(I,K,J)
        ENDDO
        CALL CALC_LPI(QL1D,QI1D,QS1D,QG1D,W1D,TEMPC,HEIGHT,LPI(I,J),T_BASE,T_TOP,KME,KTE)
      END DO
      END DO
      RETURN
      END SUBROUTINE CALCLPI



      SUBROUTINE &
     &  CALC_LPI(QL,QI,QS,QG,W,T,HEIGHT,LPI,T_BASE,T_TOP,NK,NKE)
      IMPLICIT NONE
      INTEGER NK,NKE
      REAL T_BASE,T_TOP
      REAL QL(NK)
      REAL QG(NK)
      REAL QI(NK)
      REAL QS(NK)
      REAL W(0:NK)
      REAL T(NK) ! TEMPERATURE IN CELCIUS
      REAL HEIGHT(0:NK)
      REAL LPI
      REAL DEL_Z(NK)
      REAL W_AVE(NK)
      INTEGER K
      REAL DEL_Z_TOT
      REAL TOP, BOT
      REAL NUM,DEN,AVE_Z
      REAL NUM_S,DEN_S
      REAL NUM_I,DEN_I
      REAL Q_ISG

      DO K=1,NKE
          TOP=HEIGHT(K)
          BOT=HEIGHT(K-1)
          DEL_Z(K)=TOP-BOT
          W_AVE(K)=0.5*(W(K)+W(K-1))
      END DO
!
      AVE_Z=0
      DEL_Z_TOT=0
      LPI=0
      DO K=1,NKE-1
          IF (T(K).LE.T_BASE.AND.T(K).GT.T_TOP)THEN ! SET TEMP RANGE
              Q_ISG = 0
              NUM_I = SQRT(QI(K)*QG(K))
              DEN_I = QI(K)+QG(K)
              IF (DEN_I.GT.0) THEN
                  Q_ISG = Q_ISG + QG(K) * (NUM_I/DEN_I)
              END IF

              NUM_S = SQRT(QS(K)*QG(K))   
              DEN_S = QS(K)+QG(K)
              IF (DEN_S.GT.0) THEN
                  Q_ISG = Q_ISG + QG(K) * (NUM_S/DEN_S)  ! ICE "FRACT"-CONTENT
              END IF

              NUM   = SQRT(QL(K)*Q_ISG)
              DEN   = QL(K)+Q_ISG
              IF (DEN.GT.0)THEN
                  AVE_Z = AVE_Z + DEL_Z(K)*(2.*NUM/DEN)*W_AVE(K)**2 ! LIGHTNING POTENTIAL INDEX J/UNIT-MASS
              END IF
              DEL_Z_TOT = DEL_Z_TOT + DEL_Z(K)
          END IF
      END DO
!
      IF (DEL_Z_TOT.GT.0) THEN
          LPI=AVE_Z/DEL_Z_TOT
      END IF
!
      RETURN
      END SUBROUTINE CALC_LPI
  END MODULE MODULE_LTNG_LPI
