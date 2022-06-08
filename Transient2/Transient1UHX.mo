model Transient1UHX


  SMD_MSR_Modelica.Nuclear.mPKE mPKE( 
    Lam = 2.400E-04,
    lam = {1.240E-02,3.05E-02,1.11E-01,3.01E-01,1.140E+00,3.014E+00},
    beta = {0.000223,0.001457,0.001307,0.002628,0.000766,0.00023},
    tauCoreNom = 8.46,
    tauLoopNom = 16.73) annotation(
    Placement(visible = true, transformation(origin = {-1384, 288}, extent = {{-314, -314}, {314, 314}}, rotation = 0)));
  
  SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(
    P = 8,
    DHYG1 = 9.635981959409105E-01*2.464783802008740E-03,
    DHYG2 = 3.560674858154914E-02*2.464783802008740E-03,
    DHYG3 = 7.950554775404400E-04*2.464783802008740E-03,
    TotalFuelMass = 4.093499709644400e+03,
    DHlamG1 = 0.09453,
    DHlamG2 = 0.004420,
    DHlamG3 = 8.6098E-5) annotation(
    Placement(visible = true, transformation(origin = {-1268, -284}, extent = {{-296, -296}, {296, 296}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.PrimaryHeatExchanger PHX(
    m_PN1 = 87.0959,
    m_PN2 = 87.0959,
    m_PN3 = 87.0959,
    m_PN4 = 87.0959,
    m_TN1 = 101.1662,
    m_TN2 = 101.1662,
    m_SN1 = 24.9200,
    m_SN2 = 24.9200,
    m_SN3 = 24.9200,
    m_SN4 = 24.9200,
    cP_pFluid = 1.9665E-3,
    cP_Tube = 5.778E-04,
    cP_sFluid = 2.39E-3,
    m_dot_pFluidNom = 7.5708E-02*2.14647E+03,
    m_dot_sFluidNom = 5.36265E-02*1.922E3,
    hApnNom = 0.1620,
    hAsnNom = 0.0765,
    tauPHXtoPipe = 3.85,
    tauPHXtoUHX = 4.71) annotation(
    Placement(visible = true, transformation(origin = {1677, 329}, extent = {{-427, -427}, {427, 427}}, rotation = 0)));
  
  SMD_MSR_Modelica.Nuclear.FuelChannel FuelChannel(
    P = 8,
    TotalFuelMass = 4.093499709644400e+03,
    hAnom = 0.0180,
    m_FN1 = 687.3959,
    m_FN2 = 687.3959,
    m_GN = 3.6342e+03,
    cP_fuel = 1.9665E-3,
    cP_graphite = 1.773E-3,
    mdot_fuelNom = 7.5708E-02*2.14647E+03,
    kFN1 = 0.465,
    kFN2 = 0.465,
    kG = 0.07,
    kHT_FN1 = 0.5,
    kHT_FN2 = 0.5,
    tauCoreToDHRS = 1.385) annotation(
    Placement(visible = true, transformation(origin = {-452, -162}, extent = {{-574, -574}, {574, 574}}, rotation = 0)));
    
  SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
    m_DHRS = 1.625049507600000e+02,
    Cp_fluid = 1.9665E-3,
    m_dotNom = 7.5708E-02*2.14647E+03,
    m_Pi_C_PHX = 6.126436643651999e+02,
    DHRS_timeConstant = 10,
    DHRS_MaxPowerRm = 0.04*8,
    DHRS_PowerBleed = 0,
    DHRS_time = 2000,
    tauDHRStoPHX = 1.385) annotation(
    Placement(visible = true, transformation(origin = {634, 392}, extent = {{-298, -298}, {298, 298}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.Pipe pipe(
    m_pi= 1.625049507600000e+02,
    Cp_fluidPi = 1.9665E-3,
    m_dotPiNom = 7.5708E-02*2.14647E+03,
    m_Pi_PHX_C = 1.408917923089200e+03,
    tauPiToCore = 3.835) annotation(
    Placement(visible = true, transformation(origin = {677, -743}, extent = {{481, -481}, {-481, 481}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(
    freeConvectionFF = 0.05,
    primaryPumpK = 0.2,
    tripPrimaryPump = 2000000) annotation(
    Placement(visible = true, transformation(origin = {-900, -870}, extent = {{-96, -96}, {96, 96}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(
    freeConvectionFF = 0.05,
    secondaryPumpK = 0.2, 
    tripSecondaryPump = 2000000) annotation(
    Placement(visible = true, transformation(origin = {1613, -889}, extent = {{-69, -69}, {69, 69}}, rotation = 0)));
  
  SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(
    FuelTempSetPointNode1 = 644.5000, 
    FuelTempSetPointNode2 = 6.57E+02, 
    GrapTempSetPoint = 6.756111E+02, 
    a_F = -8.71E-05, 
    a_G = -6.66E-05, 
    omega = 1, 
    sin_mag = 0, 
    step_mag = 0,
    stepInsertionTime = 2000000,
    sinInsertionTime = 2000000) annotation(
    Placement(visible = true, transformation(origin = {-762, 736}, extent = {{264, -264}, {-264, 264}}, rotation = 0)));
      
  SMD_MSR_Modelica.HeatTransport.UHX uhx(
    m_PN1UHX = 3.924401289158446e+02,
    UHXP = 8,
    cP_pFluidUHX = 2.39E-3,
    m_dot_pFluidUHXnom = 103.0701,
    tauUHXtoPHX = 8.24, 
    tripUHX = 2000) annotation(
    Placement(visible = true, transformation(origin = {2233, -477}, extent = {{-567, 567}, {567, -567}}, rotation = 0)));


equation
  connect(mPKE.n_population, decayHeat.nPop) annotation(
    Line(points = {{-1384, 169}, {-1363, 169}, {-1363, -219}}, color = {20, 36, 248}));
  connect(mPKE.n_population, FuelChannel.nPop) annotation(
    Line(points = {{-1384, 169}, {-1040.5, 169}, {-1040.5, 171}, {-911, 171}}, color = {20, 36, 248}));
  connect(primaryPump.primaryFlowFrac, mPKE.fuelFlowFrac) annotation(
    Line(points = {{-900, -928}, {-1078, -928}, {-1078, 414}, {-1277, 414}}, color = {245, 121, 0}));
  connect(primaryPump.primaryFlowFrac, FuelChannel.fuelFlowFraction) annotation(
    Line(points = {{-900, -928}, {-682, -928}, {-682, -621}}, color = {245, 121, 0}));
  connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
    Line(points = {{-900, -928}, {455, -928}, {455, 332}}, color = {245, 121, 0}));
  connect(reactivityFeedback.feedback, mPKE.feedback) annotation(
    Line(points = {{-852, 784}, {-852, 786.5}, {-1384, 786.5}, {-1384, 414}}, color = {78, 154, 6}));
  connect(FuelChannel.grapNode, reactivityFeedback.grapNode) annotation(
    Line(points = {{-164, -70}, {-42, -70}, {-42, 656}, {-656, 656}}));
  connect(FuelChannel.fuelNode2, reactivityFeedback.fuelNode2) annotation(
    Line(points = {{-164, -334}, {90, -334}, {90, 736}, {-656, 736}}));
  connect(FuelChannel.fuelNode1, reactivityFeedback.fuelNode1) annotation(
    Line(points = {{-164, -576}, {172, -576}, {172, 826}, {-656, 826}}));
  connect(FuelChannel.temp_Out, dhrs.DHRS_tempIN) annotation(
    Line(points = {{-164, 142}, {-164, 458}, {456, 458}}));
  connect(dhrs.DHRS_TempOUT, PHX.T_in_pFluid) annotation(
    Line(points = {{766, 458}, {942, 458}, {942, 457}, {1241, 457}}));
  connect(PHX.T_out_pFluid, pipe.PiTemp_IN) annotation(
    Line(points = {{2368, 458}, {2922, 458}, {2922, -618}, {908, -618}}));
  connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
    Line(points = {{-900, -928}, {458, -928}, {458, -1146}, {908, -1146}, {908, -801}}, color = {245, 121, 0}));
  connect(decayHeat.decayHeat_Out, FuelChannel.P_decay) annotation(
    Line(points = {{-1162, -372}, {-912, -372}, {-912, -196}}, color = {0, 225, 255}));
  connect(decayHeat.decayHeat_Out, dhrs.DHRS_DecayHeat) annotation(
    Line(points = {{-1162, -372}, {-1160, -372}, {-1160, -684}, {-1622, -684}, {-1622, 948}, {456, 948}, {456, 570}}, color = {0, 225, 255}));
  connect(decayHeat.decayHeat_Out, PHX.P_decay) annotation(
    Line(points = {{-1162, -372}, {-1164, -372}, {-1164, -684}, {-1622, -684}, {-1622, 948}, {1796, 948}, {1796, 542}}, color = {0, 225, 255}));
  connect(decayHeat.decayHeat_Out, pipe.PiDecay_Heat) annotation(
    Line(points = {{-1162, -372}, {-1164, -372}, {-1164, -1070}, {526, -1070}, {526, -292}, {908, -292}, {908, -454}}, color = {0, 225, 255}));
  connect(pipe.PiTemp_Out, FuelChannel.temp_In) annotation(
    Line(points = {{542, -618}, {252, -618}, {252, -774}, {-1034, -774}, {-1034, -564}, {-912, -564}}));
  connect(primaryPump.primaryFlowFrac, PHX.primaryFF) annotation(
    Line(points = {{-900, -928}, {1088, -928}, {1088, 792}, {1378, 792}, {1378, 542}}, color = {245, 121, 0}));
  connect(secondaryPump.secondaryFlowFrac, PHX.secondaryFF) annotation(
    Line(points = {{1614, -930}, {1626, -930}, {1626, -54}, {2232, -54}, {2232, 116}}, color = {245, 121, 0}));
  connect(PHX.T_out_sFluid, uhx.UHXtemp_In) annotation(
    Line(points = {{1250, 200}, {1134, 200}, {1134, -466}, {1836, -466}}));
  connect(PHX.T_in_sFluid, uhx.UHXtemp_Out) annotation(
    Line(points = {{2368, 192}, {2662, 192}, {2662, -466}, {2380, -466}}));
  connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
    Line(points = {{1614, -930}, {1836, -930}, {1836, -704}}, color = {245, 121, 0}));
protected


annotation(
    Diagram(coordinateSystem(extent = {{-1720, 1020}, {2940, -1240}})),
    version = "",
    uses);

end Transient1UHX;