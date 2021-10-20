

FUNCTION getDebrisThickness(Model, Node) RESULT(debrisTh)

  USE DefUtils
  
  IMPLICIT NONE
  ! the external variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: debrisTh
  !----------------------------------------------------------------------------
  ! internal variables
  !----------------------------------------------------------------------------
  INTEGER :: TimeStep, LastTimeStep=0, NTotSurfNodes, IterN, SurfNodeN, i, j, k, t, Refinement, ii, cc, &
               NSubgridNodes, NSubSegments, NSubTimeSteps=10, NGlacierNodesPre=1, NGlacierNodesPost, &
               GlacierStatus, TerminusNode, TerminusSubNode
  REAL(KIND=dp) :: TimeStepSize, DebrisEisRatio=0.01, MinHeight=4.9_dp, MaxSlope=0.70, YieldStress=8.0e4_dp, &
                     LSubSegments, SubTimeStepSize, DeltaX, DeltaY, DeltaH, DeltaDebris, DeltaDSource, DeltaVelo1, &
                     AdvCorrection, AdvTail=0.50, TotDebrisBeforeSource, TotDebrisAfterSource, &
                     TotDebrisAfterAdvection, TotDebrisAfterTailCut, TotDebrisAfterCorrection, TotDebrisAfterRemoving, &
                     TotDebrisAfterRemovingOutside, DisplacedDebris, DebrisLoss, FracToDisplace, Low, High, Tol, &
                     SlopeDisplacementM, SlopeDisplacementC, MaxSlopeDisplacement, MinSlopeDisplacement
  LOGICAL :: FirstTime, FirstTimeEver=.true., FirstOutside, CutTail, RandomDisplace
  LOGICAL :: DebrisSlideNode=.false., DebrisSlideFull=.false., DebrisSlideYieldStress=.false., DebrisSlidePartial=.true.
  CHARACTER(LEN=MAX_NAME_LEN) :: MonitoringFile='DebrisMonitoring.dat'
  REAL(KIND=dp), ALLOCATABLE :: DebrisPost(:), Debris(:), DebrisSource(:), MassBalance(:), &
                              Coord1(:), Coord2(:), Velo1(:), Velo2(:), Height(:), SurfSlope(:), SubgridX(:), &
                              SubgridD(:), SubgridDSource(:), SubgridDPost(:), SubgridVelo1(:), SubgridY(:), &
                              SubgridH(:), SubgridSurfSlope(:), DrivingStress(:)
  REAL, PARAMETER :: NPI=3.1415927, Rho_d=1900.0, G=9.81
  !INTEGER, POINTER :: NodeIndexes(:)
  REAL(KIND=dp), POINTER :: DebrisValues(:), MassBalanceValues(:), Coord1Values(:), & 
                          Coord2Values(:), Velo1Values(:), Velo2Values(:), HeightValues(:)
  TYPE(Variable_t), POINTER :: DebrisVariable, MassBalanceVariable, Coord1Variable, &
                             Coord2Variable, Velo1Variable, Velo2Variable, HeightVariable  
    
  SAVE FirstTime, FirstTimeEver, FirstOutside, CutTail, LastTimeStep, NTotSurfNodes, IterN, SurfNodeN, &
         MinHeight, DebrisPost, SurfSlope, Refinement, NSubgridNodes, NSubSegments, SubgridX, SubgridD, &
         SubgridDSource, SubgridDPost, SubgridVelo1, SubgridY, SubgridH, SubgridSurfSlope, DrivingStress, &
         NGlacierNodesPre, NGlacierNodesPost, TerminusNode, TerminusSubNode
  SAVE Debris, MassBalance, DebrisSource, Coord1, Coord2, Velo1, Velo2, Height
      
  TimeStepSize = GetTimeStepSize()
  TimeStep = GetTimestep()
  If (TimeStep .ne. LastTimeStep) then
    FirstTime = .true.
    LastTimeStep = TimeStep
    CutTail = .false.
    FirstOutside = .true.
  End if

  Coord1Variable => VariableGet( Model % Variables, 'Coord1Surf' )
  Coord1Values  => Coord1Variable % Values
  Coord2Variable => VariableGet( Model % Variables, 'Coord2Surf' )
  Coord2Values  => Coord2Variable % Values
  Velo1Variable => VariableGet( Model % Variables, 'Velo1Surf' )
  Velo1Values  => Velo1Variable % Values
  Velo2Variable => VariableGet( Model % Variables, 'Velo2Surf' )
  Velo2Values  => Velo2Variable % Values
  HeightVariable => VariableGet( Model % Variables, 'HeightSurf' )
  HeightValues  => HeightVariable % Values
  DebrisVariable => VariableGet( Model % Variables, 'Debris' )
  DebrisValues  => DebrisVariable % Values
  MassBalanceVariable => VariableGet( Model % Variables, 'MassBalance' )
  MassBalanceValues  => MassBalanceVariable % Values
  
  NTotSurfNodes = SIZE(DebrisValues)
  
  If (FirstTimeEver) then
    ALLOCATE(DebrisPost(NTotSurfNodes),DebrisSource(NTotSurfNodes))
    ALLOCATE(Debris(NTotSurfNodes),MassBalance(NTotSurfNodes))
    ALLOCATE(Coord1(NTotSurfNodes),Coord2(NTotSurfNodes))
    ALLOCATE(Velo1(NTotSurfNodes),Velo2(NTotSurfNodes))
    ALLOCATE(Height(NTotSurfNodes),SurfSlope(NTotSurfNodes-1))    
  End if

  Coord1(1:NTotSurfNodes) = Coord1Values(1:NTotSurfNodes)
  Coord2(1:NTotSurfNodes) = Coord2Values(1:NTotSurfNodes)
  Velo1(1:NTotSurfNodes) = Velo1Values(1:NTotSurfNodes)
  Velo2(1:NTotSurfNodes) = Velo2Values(1:NTotSurfNodes)
  Height(1:NTotSurfNodes) = HeightValues(1:NTotSurfNodes)
  Debris(1:NTotSurfNodes) = DebrisValues(1:NTotSurfNodes)
  MassBalance(1:NTotSurfNodes) = MassBalanceValues(1:NTotSurfNodes)
       
  !! Attention: this part is only valid for the special mesh used in 
  !! this debris experiment!!
  If (FirstTimeEver) then
  
    Refinement = 7   ! Please choose an odd number
   
    NSubgridNodes = (NTotSurfNodes-1)*Refinement + NTotSurfNodes
    NSubSegments = Refinement+1
    
    ALLOCATE(SubgridX(NSubgridNodes),SubgridY(NSubgridNodes),SubgridH(NSubgridNodes))
    ALLOCATE(SubgridD(NSubgridNodes),SubgridDSource(NSubgridNodes),SubgridDPost(NSubgridNodes))
    ALLOCATE(SubgridVelo1(NSubgridNodes))
    ALLOCATE(SubgridSurfSlope(NSubgridNodes),DrivingStress(NSubgridNodes))
        
    Do i = 1,NTotSurfNodes-1
      ii = (i-1)*NSubSegments+1
      SubgridX(ii) = Coord1(i)
      SubgridD(ii) = Debris(i)
      DeltaX = (Coord1(i+1)-Coord1(i))/NSubSegments
      DeltaDebris = (Debris(i+1)-Debris(i))/NSubSegments
      Do j = 1,NSubSegments-1
        SubgridX(ii+j) = Coord1(i)+DeltaX*j
        SubgridD(ii+j) = Debris(i)+DeltaDebris*j
      End do
    End do
    SubgridX(NSubgridNodes) = Coord1(NTotSurfNodes)
    SubgridD(NSubgridNodes) = Debris(NTotSurfNodes)    
    
    OPEN(1, file = MonitoringFile, status = 'new') 
    WRITE(1,*) 'BeforeSource  ', 'AfterSource  ', 'AfterAdvection  ', 'AfterTailCut  ', &
               'AfterCorrection  ', 'AfterRemoving  ', 'AfterRemovingOutside'
    CLOSE(1)
    
  End if
  
  
  If (FirstTime) then

    ! Inquire if advancing, retreating or stable glacier
    NGlacierNodesPost=1
    Do while (Height(NGlacierNodesPost)>MinHeight)
      NGlacierNodesPost=NGlacierNodesPost+1
    End do
    NGlacierNodesPost = NGlacierNodesPost-1
    
    If(NGlacierNodesPost>NGlacierNodesPre) then
      GlacierStatus = 1  !'advancing'
    Else if(NGlacierNodesPost<NGlacierNodesPre) then
      GlacierStatus = 2  !'retreating'
    Else if(NGlacierNodesPost==NGlacierNodesPre) then
      GlacierStatus = 3  !'stable'
    End if
    
    NGlacierNodesPre = NGlacierNodesPost
    TerminusNode = NGlacierNodesPost+1
  
    ! Calculate debris source due to melt
    DebrisSource = 0.0_dp
    Do i = 1,NTotSurfNodes
      If(MassBalance(i) .lt. 0.0_dp) then
        DebrisSource(i) = -1.0_dp*MassBalance(i)*TimeStepSize*DebrisEisRatio 
      Else
        DebrisSource(i) = 0.0_dp
      End if  
    End do
    
    
    ! Subsampling of Coord2, DebrisSource, Velo1
    Do i = 1,NTotSurfNodes-1
      ii = (i-1)*NSubSegments+1
      SubgridY(ii) = Coord2(i)
      SubgridH(ii) = Height(i)
      SubgridDSource(ii) = DebrisSource(i)
      SubgridVelo1(ii) = Velo1(i)
      DeltaY = (Coord2(i+1)-Coord2(i))/NSubSegments
      DeltaH = (Height(i+1)-Height(i))/NSubSegments
      DeltaDSource = (DebrisSource(i+1)-DebrisSource(i))/NSubSegments
      DeltaVelo1 = (Velo1(i+1)-Velo1(i))/NSubSegments
      Do j = 1,NSubSegments-1
        SubgridY(ii+j) = Coord2(i)+DeltaY*j
        SubgridH(ii+j) = Height(i)+DeltaH*j
        SubgridDSource(ii+j) = DebrisSource(i)+DeltaDSource*j
        SubgridVelo1(ii+j) = Velo1(i)+DeltaVelo1*j
      End do
    End do
    SubgridY(NSubgridNodes) = Coord2(NTotSurfNodes)
    SubgridH(NSubgridNodes) = Height(NTotSurfNodes)
    SubgridDSource(NSubgridNodes) = DebrisSource(NTotSurfNodes)  
    SubgridVelo1(NSubgridNodes) = Velo1(NTotSurfNodes)  
    
    ! Find terminus in the subgrid frame
    TerminusSubNode=1
    Do while (SubgridH(TerminusSubNode)>MinHeight)
      TerminusSubNode=TerminusSubNode+1
    End do
    
    TotDebrisBeforeSource = SUM(SubgridD)
    
    ! Add debris source due to melt (only upstream of terminus)
    Do i=1,NSubgridNodes
      If((i .le. TerminusSubNode) .and. (SubgridH(i) .gt. MinHeight)) then
        SubgridDPost(i) = SubgridD(i) + SubgridDSource(i)
      Else
        SubgridDPost(i) = SubgridD(i)
      End if  
    End do    
    SubgridD(1:NSubgridNodes) = SubgridDPost(1:NSubgridNodes)    
    TotDebrisAfterSource = SUM(SubgridD)
    
    !!!!!!
    ! DEBRIS ADVECTION        
    !!!!!!
    
    SubTimeStepSize = TimeStepSize/NSubTimeSteps
    SubgridDPost(1) = 0.0_dp
    Do t=1,NSubTimeSteps
      Do i = 2,NSubgridNodes
        If(SubgridH(i) .gt. MinHeight) then
          SubgridDPost(i) = SubgridD(i) + &
            (SubgridVelo1(i)*SubTimeStepSize)/(ABS(SubgridX(i)-SubgridX(i-1)))*(SubgridD(i)-SubgridD(i-1)) + &
            (SubgridD(i)*SubTimeStepSize)/(ABS(SubgridX(i)-SubgridX(i-1)))*(SubgridVelo1(i)-SubgridVelo1(i-1))
        Else 
          SubgridDPost(i) = SubgridD(i)
        End if    
      End do
      SubgridD(1:NSubgridNodes) = SubgridDPost(1:NSubgridNodes)
    End do
    
    TotDebrisAfterAdvection = SUM(SubgridD)
    
    ! Cut advection Tail
    Do i=1,NTotSurfNodes-1
      ii = (i-1)*NSubSegments+1
      If(SubgridD(ii) .gt. 0.0_dp) then
        If(.not. CutTail) then
          If( ((SubgridD(ii+Refinement+1)/SubgridD(ii)) .lt. AdvTail) .and. (SubgridD(ii+Refinement+1) .lt. 0.1) ) then
            CutTail = .true.
          End if
        Else  
          If(SubgridD(ii) .lt. 0.1) then
            SubgridDPost(ii:ii+Refinement) = 0.0_dp
          End if
        End if        
      End if
    End do
    SubgridD(1:NSubgridNodes) = SubgridDPost(1:NSubgridNodes)    
    TotDebrisAfterTailCut = SUM(SubgridD)
    
    ! Compensate advection losses
    AdvCorrection = TotDebrisAfterSource/TotDebrisAfterTailCut
  
    SubgridDPost(1:NSubgridNodes) = SubgridD(1:NSubgridNodes)*AdvCorrection
    SubgridD(1:NSubgridNodes) = SubgridDPost(1:NSubgridNodes)

    TotDebrisAfterCorrection = SUM(SubgridD)
    
    
    ! Debris removal due to slide on steep slopes (from subnode to the next subnode)
    SubgridSurfSlope = 0.0_dp     
    If(DebrisSlideNode) then
      Do i = 10,NSubgridNodes-1          !Skip first nodes with non-realistic surface shape
        SubgridSurfSlope(i) = (SubgridY(i+1)-SubgridY(i))/(SubgridX(i+1)-SubgridX(i))
        If(SubgridSurfSlope(i) .gt. MaxSlope) then
          SubgridDPost(i+1) = SubgridDPost(i+1) + SubgridDPost(i) 
          SubgridDPost(i) = 0.0_dp
        End if
      End do
    End if

    ! Debris removal due to slide on steep slopes (complete removal downstream of subnode)
    SubgridSurfSlope = 0.0_dp     
    If(DebrisSlideFull) then  
      Do i = 10,NSubgridNodes-1         !Skip first nodes with non-realistic surface shape
        SubgridSurfSlope(i) = (SubgridY(i+1)-SubgridY(i))/(SubgridX(i+1)-SubgridX(i))
        If(SubgridSurfSlope(i) .gt. MaxSlope) then
          SubgridDPost(i:NSubgridNodes) = 0.0_dp                    
          Exit
        End if
      End do
    End if    
    
    ! Debris removal due to high driving stress (complete removal downstream of subnode)
    SubgridSurfSlope = 0.0_dp
    DrivingStress = 0.0_dp
    If(DebrisSlideYieldStress) then  
      Do i = 10,NSubgridNodes-1         !Skip first nodes with non-realistic surface shape
        SubgridSurfSlope(i) = (SubgridY(i+1)-SubgridY(i))/(SubgridX(i+1)-SubgridX(i))
        DrivingStress(i) = Rho_d*G*SubgridDPost(i)*ABS(SQRT(SubgridSurfSlope(i)**2/(SubgridSurfSlope(i)**2+1)))
        If((SUM(DrivingStress(i-Refinement+1:i))/Refinement) .gt. YieldStress) then
          SubgridDPost(i-Refinement+1:NSubgridNodes) = 0.0_dp
          Exit
        End if       
      End do
    End if    
        
    ! Partial redistribution of debris due to slope and complete debris removal due to high driving stress
    RandomDisplace = .true.        !  .true. -> random displacement (within a deterministic range)
                                   !  .false. -> deterministic displacement
    SubgridSurfSlope = 0.0_dp
    DrivingStress = 0.0_dp
    FracToDisplace = 0.0_dp
    High = 0.0_dp
    Low = 0.0_dp

    DebrisLoss = 0.0_dp
    SlopeDisplacementM = 1.0_dp          ! parameter M
    SlopeDisplacementC = 0.20_dp         ! parameter C
    MaxSlopeDisplacement = 1.0_dp
    MinSlopeDisplacement = 0.0_dp
    Tol = 0.10_dp                        ! parameter \delta
    
    If(DebrisSlidePartial) then  
      Do i = 10,NSubgridNodes-1         !Skip first nodes with non-realistic surface shape
        
        SubgridSurfSlope(i) = (SubgridY(i+1)-SubgridY(i))/(SubgridX(i+1)-SubgridX(i))
        
        If(SubgridSurfSlope(i) .gt. MaxSlope) then
          
          If(RandomDisplace) then 
            High = SubgridSurfSlope(i)*SlopeDisplacementM + SlopeDisplacementC + Tol
            Low = SubgridSurfSlope(i)*SlopeDisplacementM + SlopeDisplacementC - Tol
            FracToDisplace = RAND()*(High-Low)+Low
          Else
            FracToDisplace = SubgridSurfSlope(i)*SlopeDisplacementM + SlopeDisplacementC
          End if
          
          If(FracToDisplace .gt. MaxSlopeDisplacement) then
            FracToDisplace = MaxSlopeDisplacement
          Else if(FracToDisplace .le. MinSlopeDisplacement) then
            FracToDisplace = MinSlopeDisplacement          
          End if          
                    
          If(SubgridD(i) .gt. 0.0_dp) then
            SubgridDPost(i+1) = SubgridDPost(i+1) + SubgridDPost(i)*FracToDisplace*(1.0_dp-DebrisLoss)
            SubgridDPost(i) = SubgridDPost(i)*(1.0_dp-FracToDisplace)   
          End if  
          
        End if        

        DrivingStress(i) = Rho_d*G*SubgridDPost(i)*ABS(SQRT(SubgridSurfSlope(i)**2/(SubgridSurfSlope(i)**2+1)))
        If((SUM(DrivingStress(i-Refinement+1:i))/Refinement) .gt. YieldStress) then
          SubgridDPost(i-Refinement+1:NSubgridNodes) = 0.0_dp
          Exit
        End if

        FracToDisplace = 0.0_dp 
        High = 0.0_dp
        Low = 0.0_dp
      End do     
    End if     
    SubgridD(1:NSubgridNodes) = SubgridDPost(1:NSubgridNodes)
    TotDebrisAfterRemoving = SUM(SubgridD)

    
    ! Eliminate debris if outside the glacier
    If(GlacierStatus == 2) then
      Do i=1,NTotSurfNodes-1
        ii = (i-1)*NSubSegments+1
        If(SubgridH(ii) .lt. MinHeight) then
          SubgridDPost(ii+Refinement+1:NSubgridNodes) = 0.0_dp
          Exit 
        End if      
      End do   
    End if  

    SubgridD(1:NSubgridNodes) = SubgridDPost(1:NSubgridNodes)
    TotDebrisAfterRemovingOutside = SUM(SubgridD)
    
    OPEN(1, file = MonitoringFile, status = 'old', position="append")  
    WRITE(1,*) TotDebrisBeforeSource, TotDebrisAfterSource, TotDebrisAfterAdvection, &
                 TotDebrisAfterTailCut, TotDebrisAfterCorrection, TotDebrisAfterRemoving, &
                 TotDebrisAfterRemovingOutside
    CLOSE(1) 

    ! Return to the FE mesh
    Do i = 1,NTotSurfNodes
      DebrisPost(i) = SubgridDPost((i-1)*(Refinement+1)+1)
    End do
    
    IterN = 1
    SurfNodeN = 1

  End if
  

  If(MODULO(IterN,2) .ne. 0) then
    debrisTh = DebrisPost(SurfNodeN)
    SurfNodeN = SurfNodeN+1
  End if
    
  IterN = IterN+1
  
  If(IterN .eq. ((NTotSurfNodes-1)*2)+1) then
    debrisTh = DebrisPost(SurfNodeN)
  End if
  
  FirstTime=.false.
  FirstTimeEver=.false.
    
  RETURN

END FUNCTION getDebrisThickness



FUNCTION initDebris(Model, Node, InputArray) RESULT(debrisIni)

  USE DefUtils
  
  IMPLICIT NONE
  ! the external variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: InputArray(1)
  REAL(KIND=dp) :: debrisIni
  !----------------------------------------------------------------------------
  ! internal variables
  !----------------------------------------------------------------------------
  REAL(KIND=dp) :: coord1

  coord1 = InputArray(1)
    
  IF ((coord1>2000.0) .and. (coord1<2100.0)) THEN
     debrisIni = 5.0_dp
  ELSE
     debrisIni = 0.0_dp
  END IF
  
  RETURN

END FUNCTION initDebris



FUNCTION getMassBalance(Model, Node, InputArray) RESULT(MassBalance)

  USE DefUtils
  
  IMPLICIT NONE
  ! the external variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: InputArray(3)
  REAL(KIND=dp) :: MassBalance
  !----------------------------------------------------------------------------
  ! internal variables
  !----------------------------------------------------------------------------
  REAL(KIND=dp) :: Coord2, Debris, MassBalanceCmDayDebris, MassBalanceMYearDebris
  REAL(KIND=dp) :: LW, SW, HumidityG, AirTemp, WindSpeed
  REAL(KIND=dp) :: AccG, AccMax, ReferenceElevation, InitialDebris, AblationDays
  REAL(KIND=dp) :: In, Q, K, Qm, Qh, Tm, Rhoaa, Um, Xm, Xr, Alphad, Alphai, Ustar, Ca, &
                     Lm, Lv, Tf, Eps, Rhoi, Sigma, Kstar, Gamma, PhiD, Humidity
    
  Coord2 = InputArray(1)
  Debris = InputArray(2)
        
  ! altitude gradients of the crucial parameters (radiation from Marty et al., TaAClimat; 2002)
  LW = 2.9          ! W/m^2 /100m                       2.9
  SW = 1.3          ! W/m^2 /100m                       1.3
  HumidityG = 0     ! % /100m         rough estimate
  AirTemp = 0.7     ! C /100m
  WindSpeed = 0.02  ! m/s /100m       rough estimate    0.2
  
  ! accumulation follows a linear increase above the ELA up to a plateau
  AccG = 0.1                    ! m w.e. /100m
  AccMax = 1                    ! m w.e.
  ReferenceElevation = 2200     ! m
  InitialDebris = 0             ! m
  AblationDays = 120            !
  
  In = 100                 ! Wm^-2        incoming long wave
  Q = 500                  ! Wm^-2        incoming short wave
  K = 0.585                ! Wm^-1K^-1    thermal conductivity          0.585
  Qm = 0.0012              ! kg m^-3      measured humiditiy level
  Qh = 0.006               ! kg m^-3      saturated humidity level
  Tm = 2                   ! C            air temperature
  Rhoaa = 1.22             ! kgm^-3       air densitiy
  Um = 1.5                 ! ms^-1        measured wind speed
  Xm = 1.5                 ! ms^-1        measurement height
  Xr = 0.01                ! ms^-1        surface roughness             0.01
  Alphad = 0.07            !              debris albedo                 0.07
  Alphai = 0.4             !              ice ablbedo
  Ustar = 0.16             ! ms^-1        friction velocity             0.16
  Ca = 1000                ! jkg^-1K^-1   specific heat capacity of air
  Lm = 3.34E+05            ! jkg^-1K^-1   latent heat of ice melt
  Lv = 2.50E+06            ! jkg^-1K^-1   latent heat of evaporation
  Tf = 273                 ! K            water freeezing temperature
  Eps = 0.95               !              thermal emissivity
  Rhoi = 900               ! kgm^-3       ice density
  Sigma = 5.67E-08         ! Wm^-2K^-4    Stefan Boltzmann constant
  Kstar = 0.4              !              von kármán constant
  Gamma = 180              ! m^-1         wind speed attenuation        234
  PhiD = 0.01              !              debris packing fraction       0.01
  Humidity = 0.2           !              relative humidity

  MassBalanceCmDayDebris = (((In-(Coord2-ReferenceElevation)*LW/100)-(Eps*Sigma*(Tf*Tf*Tf*Tf))+ &
                  (Q+(Coord2-ReferenceElevation)*SW/100)*(1-Alphad)+ &
                  (Rhoaa*Ca*Ustar*Ustar)/((Um-(Coord2-ReferenceElevation)* &
                  WindSpeed/100)-Ustar*(2-(EXP(Gamma*Xr))))*(Tm-(Coord2- &
                  ReferenceElevation)*AirTemp/100))/((1-PhiD)*Rhoi*Lm)/(1+ &
                  ((Rhoaa*Ca*Ustar*Ustar)/((Um-(Coord2-ReferenceElevation)* &
                  WindSpeed/100)-Ustar*(2-(EXP(Gamma*Xr))))+4*Eps*Sigma*(Tf*Tf*Tf))/ &
                  K*Debris)-(Lv*Ustar*Ustar*(Qh-(Qh*(Humidity-(Coord2- &
                  ReferenceElevation)*HumidityG/100)))*(EXP(-Gamma*Xr)))/((1-PhiD)* &
                  Rhoi*Lm*Ustar)/((((Um-(Coord2-ReferenceElevation)*WindSpeed/100) &
                  -2*Ustar)*EXP(-Gamma*Xr))/Ustar+EXP(Gamma*Debris)))*100*24*60*60      
  
  MassBalanceMYearDebris = -MassBalanceCmDayDebris/100*AblationDays
  MassBalance = MassBalanceMYearDebris
  
  RETURN

END FUNCTION getMassBalance



FUNCTION initMassBalance(Model, Node, InputArray) RESULT(massbalanceIni)

  USE DefUtils
  
  IMPLICIT NONE
  ! the external variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: InputArray
  REAL(KIND=dp) :: massbalanceIni
  !----------------------------------------------------------------------------
  ! internal variables
  !----------------------------------------------------------------------------
  REAL(KIND=dp) :: coord1, zs, MassBalanceCmDay, MassBalanceMYear
  REAL(KIND=dp) :: LW, SW, HumidityG, AirTemp, WindSpeed                   
  REAL(KIND=dp) :: AccG, AccMax, ReferenceElevation, InitialDebris, AblationDays
  REAL(KIND=dp) :: In, Q, K, Qm, Qh, Tm, Rhoaa, Um, Xm, Xr, Alphad, Alphai, Ustar, Ca, &
                     Lm, Lv, Tf, Eps, Rhoi, Sigma, Kstar, Gamma, PhiD, Humidity

  zs = InputArray
  
  ! altitude gradients of the crucial parameters (radiation from Marty et al., TaAClimat; 2002)
  LW = 2.9          ! W/m^2 /100m                       2.9
  SW = 1.3          ! W/m^2 /100m                       1.3
  HumidityG = 0     ! % /100m         rough estimate
  AirTemp = 0.7     ! C /100m
  WindSpeed = 0.02  ! m/s /100m       rough estimate    0.2
  
  ! accumulation follows a linear increase above the ELA up to a plateau
  AccG = 0.1                    ! m w.e. /100m
  AccMax = 1                    ! m w.e.
  ReferenceElevation = 2200  !-700     ! m
  InitialDebris = 0             ! m
  AblationDays = 120            !
  
  In = 100                 ! Wm^-2        incoming long wave
  Q = 500                  ! Wm^-2        incoming short wave
  K = 0.585                ! Wm^-1K^-1    thermal conductivity          0.585
  Qm = 0.0012              ! kg m^-3      measured humiditiy level
  Qh = 0.006               ! kg m^-3      saturated humidity level
  Tm = 2                   ! C            air temperature
  Rhoaa = 1.22             ! kgm^-3       air densitiy
  Um = 1.5                 ! ms^-1        measured wind speed
  Xm = 1.5                 ! ms^-1        measurement height
  Xr = 0.01                ! ms^-1        surface roughness             0.01
  Alphad = 0.07            !              debris albedo                 0.07
  Alphai = 0.4             !              ice ablbedo
  Ustar = 0.16             ! ms^-1        friction velocity             0.16
  Ca = 1000                ! jkg^-1K^-1   specific heat capacity of air
  Lm = 3.34E+05            ! jkg^-1K^-1   latent heat of ice melt
  Lv = 2.50E+06            ! jkg^-1K^-1   latent heat of evaporation
  Tf = 273                 ! K            water freeezing temperature
  Eps = 0.95               !              thermal emissivity
  Rhoi = 900               ! kgm^-3       ice density
  Sigma = 5.67E-08         ! Wm^-2K^-4    Stefan Boltzmann constant
  Kstar = 0.4              !              von kármán constant
  Gamma = 180              ! m^-1         wind speed attenuation        234
  PhiD = 0.01              !              debris packing fraction       0.01
  Humidity = 0.2           !              relative humidity
  
  MassBalanceCmDay = (((In-(zs-ReferenceElevation)*LW/100)-(Eps*Sigma*(Tf*Tf*Tf*Tf))+ &
                        (Q+(zs-ReferenceElevation)*SW/100)*(1-Alphad)+ &
                        (Rhoaa*Ca*Ustar*Ustar)/((Um-(zs-ReferenceElevation)* &
                        WindSpeed/100)-Ustar*(2-(EXP(Gamma*Xr))))*(Tm-(zs- &
                        ReferenceElevation)*AirTemp/100))/((1-PhiD)*Rhoi*Lm)/(1+ &
                        ((Rhoaa*Ca*Ustar*Ustar)/((Um-(zs-ReferenceElevation)* &
                        WindSpeed/100)-Ustar*(2-(EXP(Gamma*Xr))))+4*Eps*Sigma*(Tf*Tf*Tf))/ &
                        K*InitialDebris)-(Lv*Ustar*Ustar*(Qh-(Qh*(Humidity-(zs- &
                        ReferenceElevation)*HumidityG/100)))*(EXP(-Gamma*Xr)))/((1-PhiD)* &
                        Rhoi*Lm*Ustar)/((((Um-(zs-ReferenceElevation)*WindSpeed/100) &
                        -2*Ustar)*EXP(-Gamma*Xr))/Ustar+EXP(Gamma*InitialDebris)))*100*24*60*60
  
  MassBalanceMYear = -MassBalanceCmDay/100*AblationDays

  massbalanceIni = MassBalanceMYear
  
  RETURN

END FUNCTION initMassBalance



FUNCTION getMinZs(Model, Node, InputArray) RESULT(MinZs)

  USE DefUtils
  
  IMPLICIT NONE
  ! the external variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: InputArray
  REAL(KIND=dp) :: MinZs
  !----------------------------------------------------------------------------
  ! internal variables
  !----------------------------------------------------------------------------
  REAL(KIND=dp) :: RefZs

  RefZs = InputArray
  MinZs = RefZs-0.1
  !MinZs = RefZs-4.0
 
  RETURN

END FUNCTION getMinZs



FUNCTION getMaxZs(Model, Node, InputArray) RESULT(MaxZs)

  USE DefUtils
  
  IMPLICIT NONE
  ! the external variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: InputArray
  REAL(KIND=dp) :: MaxZs
  !----------------------------------------------------------------------------
  ! internal variables
  !----------------------------------------------------------------------------
  REAL(KIND=dp) :: RefZs

  RefZs = InputArray
  MaxZs = RefZs+600.0
 
  RETURN

END FUNCTION getMaxZs


