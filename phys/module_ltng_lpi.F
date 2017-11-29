!WRF:MODEL_LAYER:PHYSICS
!
MODULE MODULE_LTNG_LPI

CONTAINS
!===================================================================
!
      SUBROUTINE CALCLPI(                                         &
                         QV,QC,QR,QI,QS,QG,QH                     &
                        ,W,Z,DZ8W,PI_PHY,TH_PHY,P_PHY,RHO_PHY     &
                        ,LPI                                      &
                        ,IDS,IDE, JDS,JDE, KDS,KDE                & ! domain
                        ,IMS,IME, JMS,JME, KMS,KME                & ! memory
                        ,ITS,ITE, JTS,JTE, KTS,KTE                & ! tile
                                                                  )
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      REAL, PARAMETER :: T_0 = 273.15
!-------------------------------------------------------------------

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

      REAL, INTENT(OUT), DIMENSION(IMS:IME,JMS:JME)::       LPI

! DEFINE LOCAL DATA
      REAL, DIMENSION(KTS:KTE):: TEMPC
      REAL, DIMENSION(KTS:KTE):: QL1D,QI1D,QS1D,QG1D,W1D,DZ1D

      INTEGER ::               I,J,K

      REAL, PARAMETER :: T_BASE= 0
      REAL, PARAMETER :: T_TOP=-20

      CHARACTER*256 :: MESSAGE

!-------------------------------------------------------------------
      DO J = JTS,JTE
          DO I = ITS,ITE
            DO K = KTS,KTE
              W1D(K)=W(I,K,J)
              DZ1D(K)=DZ8W(I,K,J)
              TEMPC(K) = TH_PHY(I,K,J)*PI_PHY(I,K,J)-T_0
              QL1D(K)=QC(I,K,J)+QR(I,K,J)
              QI1D(K)=QI(I,K,J)
              QS1D(K)=QS(I,K,J)
              ! IGNORE QH. HAIL DOESN'T USUALLY CHARGE
              QG1D(K)=QG(I,K,J)
            ENDDO
            LPI(I,J) = CALC_LPI(QL1D,QI1D,QS1D,QG1D,W1D,TEMPC,DZ1D,T_BASE,T_TOP,KTS,KTE)
          END DO
      END DO
      RETURN
      END SUBROUTINE CALCLPI



      FUNCTION CALC_LPI(QL,QI,QS,QG,W,T,DZ,T_BASE,T_TOP,KTS,KTE)
      IMPLICIT NONE
      INTEGER, INTENT(IN   ) :: KTS,KTE
      REAL,    INTENT(IN   ) :: T_BASE,T_TOP
      REAL, DIMENSION(KTS:KTE), INTENT(IN) :: QL,QI,QS,QG,W
      REAL, DIMENSION(KTS:KTE), INTENT(IN) :: T     ! TEMPERATURE IN CELSIUS
      REAL, DIMENSION(KTS:KTE), INTENT(IN) :: DZ    ! THICKNESS OF LAYER
      REAL :: CALC_LPI


!...DEFINE LOCAL VARIABLES

      INTEGER K
      REAL DEL_Z_TOT
      REAL NUM,DEN,AVE_Z
      REAL NUM_S,DEN_S
      REAL NUM_I,DEN_I
      REAL Q_ISG

      AVE_Z=0
      DEL_Z_TOT=0
      CALC_LPI=0
      DO K=KTS,KTE
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
                  AVE_Z = AVE_Z + DZ(K)*(2.*NUM/DEN)*W(K)**2 ! LIGHTNING POTENTIAL INDEX J/UNIT-MASS
              END IF
              DEL_Z_TOT = DEL_Z_TOT + DZ(K)
          END IF
      END DO
!
      IF (DEL_Z_TOT.GT.0) THEN
          CALC_LPI=AVE_Z/DEL_Z_TOT
      END IF
!
      RETURN
      END FUNCTION CALC_LPI
  END MODULE MODULE_LTNG_LPI
