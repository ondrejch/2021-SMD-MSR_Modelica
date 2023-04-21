package Transient4
  model MSREoneR
    SMD_MSR_Modelica.Nuclear.mPKE mPKE(Lam = 2.400E-04, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, lam = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, omega = 10, rho_0 = 0.002465140767843, sinInsertionTime = 20000000, sin_mag = 0, stepInsertionTime = 2000, step_mag = -1400E-5, tauCoreNom = 8.46, tauLoopNom = 16.73) annotation(
      Placement(visible = true, transformation(origin = {-1384, 288}, extent = {{-314, -314}, {314, 314}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(P = 8, DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, TotalFuelMass = 4.093499709644400e+03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
      Placement(visible = true, transformation(origin = {-1268, -284}, extent = {{-296, -296}, {296, 296}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryHeatExchanger PHX(m_PN1 = 87.0959, m_PN2 = 87.0959, m_PN3 = 87.0959, m_PN4 = 87.0959, m_TN1 = 101.1662, m_TN2 = 101.1662, m_SN1 = 24.9200, m_SN2 = 24.9200, m_SN3 = 24.9200, m_SN4 = 24.9200, cP_pFluid = 1.9665E-3, cP_Tube = 5.778E-04, cP_sFluid = 2.39E-3, m_dot_pFluidNom = 7.5708E-02*2.14647E+03, m_dot_sFluidNom = 5.36265E-02*1.922E3, hApnNom = 0.1620, hAsnNom = 0.0765, tauPHXtoPipe = 3.85, tauPHXtoUHX = 4.71) annotation(
      Placement(visible = true, transformation(origin = {1677, 329}, extent = {{-427, -427}, {427, 427}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel FuelChannel(P = 8, TF1_0 = 644.5000, TF2_0 = 6.57E+02, TG_0 = 6.756111E+02, cP_fuel = 1.9665E-3, cP_graphite = 1.773E-3, hAnom = 0.0180, kFN1 = 0.465, kFN2 = 0.465, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 687.3959, m_FN2 = 687.3959, m_GN = 3.6342e+03, mdot_fuelNom = 7.5708E-02*2.14647E+03, regionFlowFrac = 1, tauCoreToDHRS = 1.385) annotation(
      Placement(visible = true, transformation(origin = {-452, -162}, extent = {{-574, -574}, {574, 574}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Cp_fluid = 1.9665E-3, DHRS_MaxPowerRm = 0.08*8, DHRS_PowerBleed = 0, DHRS_time = 2000, DHRS_timeConstant = 10, m_DHRS = 1.625049507600000e+02, m_Pi_C_PHX = 6.126436643651999e+02, m_dotNom = 7.5708E-02*2.14647E+03, tauDHRStoPHX = 1.385) annotation(
      Placement(visible = true, transformation(origin = {634, 392}, extent = {{-298, -298}, {298, 298}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(m_pi = 1.625049507600000e+02, Cp_fluidPi = 1.9665E-3, m_dotPiNom = 7.5708E-02*2.14647E+03, m_Pi_PHX_C = 1.408917923089200e+03, tauPiToCore = 3.835) annotation(
      Placement(visible = true, transformation(origin = {677, -743}, extent = {{481, -481}, {-481, 481}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.05, primaryPumpK = 0.2, tripPrimaryPump = 2000) annotation(
      Placement(visible = true, transformation(origin = {-626, -1274}, extent = {{-386, -386}, {386, 386}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.05, secondaryPumpK = 0.2, tripSecondaryPump = 2000) annotation(
      Placement(visible = true, transformation(origin = {1716, -1240}, extent = {{-362, -362}, {362, 362}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 644.5000, FuelTempSetPointNode2 = 6.57E+02, GrapTempSetPoint = 6.756111E+02, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -8.71E-05, a_G = -6.66E-05) annotation(
      Placement(visible = true, transformation(origin = {-762, 736}, extent = {{264, -264}, {-264, 264}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(SetDemand = 0, UHXP = 8, cP_pFluidUHX = 2.39E-3, m_PN1UHX = 3.924401289158446e+02, m_dot_pFluidUHXnom = 103.0701, tauUHXtoPHX = 8.24, tripUHX = 2000) annotation(
      Placement(visible = true, transformation(origin = {2233, -477}, extent = {{-567, 567}, {567, -567}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 8, TotalFuelMass = 4.093499709644400e+03) annotation(
      Placement(visible = true, transformation(origin = {-2139, 175}, extent = {{329, -329}, {-329, 329}}, rotation = 0)));
  equation
    connect(mPKE.n_population, decayHeat.nPop) annotation(
      Line(points = {{-1384, 169}, {-1363, 169}, {-1363, -219}}, color = {20, 36, 248}));
    connect(mPKE.n_population, FuelChannel.nPop) annotation(
      Line(points = {{-1384, 169}, {-1040.5, 169}, {-1040.5, 171}, {-911, 171}}, color = {20, 36, 248}));
    connect(primaryPump.primaryFlowFrac, mPKE.fuelFlowFrac) annotation(
      Line(points = {{-626, -1506}, {-1078, -1506}, {-1078, 414}, {-1277, 414}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, FuelChannel.fuelFlowFraction) annotation(
      Line(points = {{-626, -1506}, {-682, -1506}, {-682, -621}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
      Line(points = {{-626, -1506}, {455, -1506}, {455, 332}}, color = {245, 121, 0}));
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
      Line(points = {{-626, -1506}, {458, -1506}, {458, -1146}, {908, -1146}, {908, -801}}, color = {245, 121, 0}));
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
      Line(points = {{-626, -1506}, {1088, -1506}, {1088, 792}, {1378, 792}, {1378, 542}}, color = {245, 121, 0}));
    connect(secondaryPump.secondaryFlowFrac, PHX.secondaryFF) annotation(
      Line(points = {{1716, -1457}, {1626, -1457}, {1626, -54}, {2232, -54}, {2232, 116}}, color = {245, 121, 0}));
    connect(PHX.T_out_sFluid, uhx.UHXtemp_In) annotation(
      Line(points = {{1250, 200}, {1134, 200}, {1134, -466}, {1836, -466}}));
    connect(PHX.T_in_sFluid, uhx.UHXtemp_Out) annotation(
      Line(points = {{2368, 192}, {2662, 192}, {2662, -466}, {2380, -466}}));
    connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
      Line(points = {{1716, -1457}, {1836, -1457}, {1836, -704}}, color = {245, 121, 0}));
    connect(mPKE.n_population, powerBlock.nPop) annotation(
      Line(points = {{-1384, 168}, {-1634, 168}, {-1634, -124}, {-1968, -124}, {-1968, 44}}, color = {20, 36, 248}));
    connect(decayHeat.decayHeat_Out, powerBlock.P_decay) annotation(
      Line(points = {{-1162, -372}, {-2244, -372}, {-2244, 44}}, color = {0, 225, 255}));
  protected
    annotation(
      Diagram(coordinateSystem(extent = {{-2480, 1020}, {2940, -1680}})),
      version = "",
      uses,
      __OpenModelica_simulationFlags(lv = "LOG_STATS", maxStepSize = "0.01", s = "dassl"));
  end MSREoneR;

  model nineE_MSRE
    parameter Real CpCoefFuel = 1;
    parameter Real HTCcoefCore = 1;
    parameter Real HTCcoefPHX = 1;
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.05, primaryPumpK = 0.2, tripPrimaryPump = 2000) annotation(
      Placement(visible = true, transformation(origin = {-626, -1274}, extent = {{-386, -386}, {386, 386}}, rotation = 0)));
    MSREnineR.nineRcore nineRcore(IF1 = {0.02168, 0.02197, 0.07897, 0.08249, 0.02254, 0.08255, 0.08623, 0.02745, 0.06936}, IF2 = {0.02678, 0.06519, 0.08438, 0.04124, 0.06801, 0.08823, 0.04290, 0.05529, 0.03473}, IG = {0.04443, 0.08835, 0.16671, 0.12077, 0.09181, 0.17429, 0.12612, 0.08408, 0.10343}, Lam = 2.400E-04, P = 8, TF1_0 = {641.0644, 636.1926, 662.2593, 692.4323, 634.6308, 650.7339, 669.3933, 633.4354, 638.9545}, TF2_0 = {651.9893, 647.9801, 677.5167, 699.8892, 641.9048, 660.1706, 673.9817, 635.8835, 640.4922}, TG_0 = {664.0462, 653.7083, 692.6790, 714.5485, 648.7869, 674.7077, 686.8285, 640.8657, 647.7290}, Tmix_0 = 6.5727E+02, TotalFuelMass = 4.093499709644400e+03, aF = -8.71E-05, aG = -6.66E-05, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, cP_fuel = 0.0020*CpCoefFuel, cP_grap = 0.0018, flowFracRegions = {0.061409895150439, 0.138549617100260, 0.234231239311780, 0.565809248437521}, hA = {0.0007056, 0.0021672, 0.00162, 0.0021132, 0.0035586, 0.002745, 0.003573, 0.009801, 0.009648}*HTCcoefCore, kFN1 = {0.01493, 0.02736, 0.04504, 0.05126, 0.03601, 0.06014, 0.06845, 0.06179, 0.09333}, kFN2 = {0.01721, 0.04550, 0.04656, 0.04261, 0.06069, 0.06218, 0.05664, 0.07707, 0.07311}, kHT1 = {0.000946, 0.001685, 0.003029, 0.003447, 0.002216, 0.004044, 0.004603, 0.003920, 0.006277}, kHT2 = {0.001081, 0.003060, 0.003131, 0.002395, 0.004081, 0.004182, 0.003184, 0.005183, 0.004305}, lam = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, mDot = 162.3880, mFN1 = {13.8215, 46.8650, 25.6293, 32.0366, 79.2677, 43.2952, 54.1876, 217.8490, 147.8261}, mFN2 = {14.4622, 31.9451, 25.6293, 62.4256, 54.1876, 43.2952, 105.4462, 126.6819, 248.0549}, mGN = {71.0660, 214.6193, 163.0457, 208.7310, 363.0457, 275.9391, 353.0964, 975.8376, 956.4467}, mMix = 324.7760, omega = 1, regionCoastDown = 0.02, regionFreeConvFF = 0.01, regionTripTime = {200000, 200000, 200000, 200000}, rho_0 = 0.002465140767843, sinReact = 0, sinTime = 20000000, stepReact = -1400E-5, stepTime = 2000, tauCore = 8.46, tauCoreDHRS = 1.3850, tauLoop = 16.73) annotation(
      Placement(visible = true, transformation(origin = {-1361, -215}, extent = {{-399, -399}, {399, 399}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.05, secondaryPumpK = 0.2, tripSecondaryPump = 2000) annotation(
      Placement(visible = true, transformation(origin = {1204, -1400}, extent = {{-362, -362}, {362, 362}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(SetDemand = 0, UHXP = 8, cP_pFluidUHX = 2.39E-3, m_PN1UHX = 3.924401289158446e+02, m_dot_pFluidUHXnom = 100.5793, tauUHXtoPHX = 8.24, tripUHX = 2000) annotation(
      Placement(visible = true, transformation(origin = {2441, -231}, extent = {{-567, 567}, {567, -567}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryHeatExchanger PHX(T_PN1_initial = 6.478750000000000e+02, cP_Tube = 5.778E-04, cP_pFluid = 0.0020*CpCoefFuel, cP_sFluid = 2.39E-3, hApnNom = 0.1620*HTCcoefPHX, hAsnNom = 0.0765, m_PN1 = 87.0959, m_PN2 = 87.0959, m_PN3 = 87.0959, m_PN4 = 87.0959, m_SN1 = 24.9200, m_SN2 = 24.9200, m_SN3 = 24.9200, m_SN4 = 24.9200, m_TN1 = 101.1662, m_TN2 = 101.1662, m_dot_pFluidNom = 162.3880, m_dot_sFluidNom = 100.5793, tauPHXtoPipe = 3.8350, tauPHXtoUHX = 4.71) annotation(
      Placement(visible = true, transformation(origin = {1691, 681}, extent = {{-427, -427}, {427, 427}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Cp_fluid = 0.0020*CpCoefFuel, DHRS_MaxPowerRm = 0.08*8, DHRS_PowerBleed = 0, DHRS_time = 2000, DHRS_timeConstant = 0.02, m_DHRS = 162.3880, m_Pi_C_PHX = 162.3880*3.7700, m_dotNom = 162.3880, tauDHRStoPHX = 1.3850) annotation(
      Placement(visible = true, transformation(origin = {-183, 825}, extent = {{509, -509}, {-509, 509}}, rotation = 180)));
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Cp_fluidPi = 0.0020*CpCoefFuel, m_Pi_PHX_C = 162.3880*8.6700, m_dotPiNom = 162.3880, m_pi = 162.3880, tauPiToCore = 3.8350) annotation(
      Placement(visible = true, transformation(origin = {-86, -892}, extent = {{528, -528}, {-528, 528}}, rotation = 0)));
  equation
    connect(primaryPump.primaryFlowFrac, nineRcore.flowFraction_In) annotation(
      Line(points = {{-626, -1506}, {-1632, -1506}, {-1632, -518}}, color = {245, 121, 0}));
    connect(secondaryPump.secondaryFlowFrac, PHX.secondaryFF) annotation(
      Line(points = {{1204, -1618}, {1498, -1618}, {1498, 314}, {2246, 314}, {2246, 468}}, color = {245, 121, 0}));
    connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
      Line(points = {{1204, -1618}, {1764, -1618}, {1764, -458}, {2044, -458}}, color = {245, 121, 0}));
    connect(uhx.UHXtemp_In, PHX.T_out_sFluid) annotation(
      Line(points = {{2044, -220}, {1264, -220}, {1264, 552}}));
    connect(PHX.T_in_sFluid, uhx.UHXtemp_Out) annotation(
      Line(points = {{2382, 544}, {2796, 544}, {2796, -220}, {2588, -220}}));
    connect(primaryPump.primaryFlowFrac, PHX.primaryFF) annotation(
      Line(points = {{-626, -1506}, {510, -1506}, {510, 894}, {1392, 894}}, color = {245, 121, 0}));
    connect(nineRcore.decayHeat_Out, PHX.P_decay) annotation(
      Line(points = {{-1620, 76}, {-1582, 76}, {-1582, 1104}, {1810, 1104}, {1810, 894}}, color = {0, 225, 255}));
    connect(dhrs.DHRS_TempOUT, PHX.T_in_pFluid) annotation(
      Line(points = {{40, 714}, {1256, 714}, {1256, 810}}));
    connect(nineRcore.decayHeat_Out, dhrs.DHRS_DecayHeat) annotation(
      Line(points = {{-1620, 76}, {-1594, 76}, {-1594, 520}, {-488, 520}}, color = {0, 225, 255}));
    connect(nineRcore.temp_Out, dhrs.DHRS_tempIN) annotation(
      Line(points = {{-1102, 76}, {-1026, 76}, {-1026, 714}, {-488, 714}}));
    connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
      Line(points = {{-626, -1506}, {-722, -1506}, {-722, 926}, {-488, 926}}, color = {245, 121, 0}));
    connect(nineRcore.decayHeat_Out, pipe.PiDecay_Heat) annotation(
      Line(points = {{-1620, 76}, {-1414, 76}, {-1414, 142}, {168, 142}, {168, -575}, {167, -575}}, color = {0, 225, 255}));
    connect(PHX.T_out_pFluid, pipe.PiTemp_IN) annotation(
      Line(points = {{2382, 810}, {3028, 810}, {3028, -820}, {167, -820}, {167, -755}}));
    connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
      Line(points = {{-626, -1506}, {168, -1506}, {168, -956}}, color = {245, 121, 0}));
    connect(pipe.PiTemp_Out, nineRcore.temp_In) annotation(
      Line(points = {{-234, -754}, {-1090, -754}, {-1090, -510}}));
    annotation(
      Diagram(coordinateSystem(extent = {{-1480, 1160}, {3040, -1780}})));
  end nineE_MSRE;

  model nineE_MSREnoDH
    parameter Real CpCoefFuel = 1;
    parameter Real HTCcoefCore = 1;
    parameter Real HTCcoefPHX = 1;
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.05, primaryPumpK = 0.2, tripPrimaryPump = 2000) annotation(
      Placement(visible = true, transformation(origin = {-626, -1274}, extent = {{-386, -386}, {386, 386}}, rotation = 0)));
    MSREnineR.nineRcoreNoDH nineRcoreNoDH(IF1 = {0.02168, 0.02197, 0.07897, 0.08249, 0.02254, 0.08255, 0.08623, 0.02745, 0.06936}, IF2 = {0.02678, 0.06519, 0.08438, 0.04124, 0.06801, 0.08823, 0.04290, 0.05529, 0.03473}, IG = {0.04443, 0.08835, 0.16671, 0.12077, 0.09181, 0.17429, 0.12612, 0.08408, 0.10343}, Lam = 2.400E-04, P = 8, TF1_0 = {641.0644, 636.1926, 662.2593, 692.4323, 634.6308, 650.7339, 669.3933, 633.4354, 638.9545}, TF2_0 = {651.9893, 647.9801, 677.5167, 699.8892, 641.9048, 660.1706, 673.9817, 635.8835, 640.4922}, TG_0 = {664.0462, 653.7083, 692.6790, 714.5485, 648.7869, 674.7077, 686.8285, 640.8657, 647.7290}, Tmix_0 = 663.5300, TotalFuelMass = 4.093499709644400e+03, aF = -8.71E-05, aG = -6.66E-05, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, cP_fuel = 0.0020*CpCoefFuel, cP_grap = 0.0018, flowFracRegions = {0.061409895150439, 0.138549617100260, 0.234231239311780, 0.565809248437521}, hA = {0.0007056, 0.0021672, 0.00162, 0.0021132, 0.0035586, 0.002745, 0.003573, 0.009801, 0.009648}, kFN1 = {0.01493, 0.02736, 0.04504, 0.05126, 0.03601, 0.06014, 0.06845, 0.06179, 0.09333}, kFN2 = {0.01721, 0.04550, 0.04656, 0.04261, 0.06069, 0.06218, 0.05664, 0.07707, 0.07311}, kHT1 = {0.000946, 0.001685, 0.003029, 0.003447, 0.002216, 0.004044, 0.004603, 0.003920, 0.006277}, kHT2 = {0.001081, 0.003060, 0.003131, 0.002395, 0.004081, 0.004182, 0.003184, 0.005183, 0.004305}, lam = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, mDot = 162.3880, mFN1 = {13.8215, 46.8650, 25.6293, 32.0366, 79.2677, 43.2952, 54.1876, 217.8490, 147.8261}, mFN2 = {14.4622, 31.9451, 25.6293, 62.4256, 54.1876, 43.2952, 105.4462, 126.6819, 248.0549}, mGN = {71.0660, 214.6193, 163.0457, 208.7310, 363.0457, 275.9391, 353.0964, 975.8376, 956.4467}, mMix = 324.7760, omega = 1, regionCoastDown = 0.02, regionFreeConvFF = 0.01, regionTripTime = {200000, 200000, 200000, 200000}, rho_0 = 0.002465140767843, sinReact = 0, sinTime = 20000000, stepReact = -1400E-5, stepTime = 2000, tauCore = 8.46, tauCoreDHRS = 1.3850, tauLoop = 16.73) annotation(
      Placement(visible = true, transformation(origin = {-1355, -215}, extent = {{-399, -399}, {399, 399}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.05, secondaryPumpK = 0.2, tripSecondaryPump = 2000) annotation(
      Placement(visible = true, transformation(origin = {1204, -1400}, extent = {{-362, -362}, {362, 362}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(SetDemand = 0, UHXP = 8, cP_pFluidUHX = 2.39E-3, m_PN1UHX = 3.924401289158446e+02, m_dot_pFluidUHXnom = 100.5793, tauUHXtoPHX = 8.24, tripUHX = 2000) annotation(
      Placement(visible = true, transformation(origin = {2441, -231}, extent = {{-567, 567}, {567, -567}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryHeatExchanger PHX(T_PN1_initial = 6.478750000000000e+02, cP_Tube = 5.778E-04, cP_pFluid = 0.0020*CpCoefFuel, cP_sFluid = 2.39E-3, hApnNom = 0.1620*HTCcoefPHX, hAsnNom = 0.0765, m_PN1 = 87.0959, m_PN2 = 87.0959, m_PN3 = 87.0959, m_PN4 = 87.0959, m_SN1 = 24.9200, m_SN2 = 24.9200, m_SN3 = 24.9200, m_SN4 = 24.9200, m_TN1 = 101.1662, m_TN2 = 101.1662, m_dot_pFluidNom = 162.3880, m_dot_sFluidNom = 100.5793, tauPHXtoPipe = 3.8350, tauPHXtoUHX = 4.71) annotation(
      Placement(visible = true, transformation(origin = {1691, 681}, extent = {{-427, -427}, {427, 427}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Cp_fluid = 0.0020*CpCoefFuel, DHRS_MaxPowerRm = 0.08*8, DHRS_PowerBleed = 0, DHRS_time = 2000, DHRS_timeConstant = 0.02, m_DHRS = 162.3880, m_Pi_C_PHX = 162.3880*3.7700, m_dotNom = 162.3880, tauDHRStoPHX = 1.3850) annotation(
      Placement(visible = true, transformation(origin = {-183, 825}, extent = {{509, -509}, {-509, 509}}, rotation = 180)));
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Cp_fluidPi = 0.0020*CpCoefFuel, m_Pi_PHX_C = 162.3880*8.6700, m_dotPiNom = 162.3880, m_pi = 162.3880, tauPiToCore = 3.8350) annotation(
      Placement(visible = true, transformation(origin = {-86, -892}, extent = {{528, -528}, {-528, 528}}, rotation = 0)));
  equation
    connect(secondaryPump.secondaryFlowFrac, PHX.secondaryFF) annotation(
      Line(points = {{1204, -1618}, {1498, -1618}, {1498, 314}, {2246, 314}, {2246, 468}}, color = {245, 121, 0}));
    connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
      Line(points = {{1204, -1618}, {1764, -1618}, {1764, -458}, {2044, -458}}, color = {245, 121, 0}));
    connect(uhx.UHXtemp_In, PHX.T_out_sFluid) annotation(
      Line(points = {{2044, -220}, {1264, -220}, {1264, 552}}));
    connect(PHX.T_in_sFluid, uhx.UHXtemp_Out) annotation(
      Line(points = {{2382, 544}, {2796, 544}, {2796, -220}, {2588, -220}}));
    connect(primaryPump.primaryFlowFrac, PHX.primaryFF) annotation(
      Line(points = {{-626, -1506}, {510, -1506}, {510, 894}, {1392, 894}}, color = {245, 121, 0}));
    connect(dhrs.DHRS_TempOUT, PHX.T_in_pFluid) annotation(
      Line(points = {{40, 714}, {1256, 714}, {1256, 810}}));
    connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
      Line(points = {{-626, -1506}, {-722, -1506}, {-722, 926}, {-488, 926}}, color = {245, 121, 0}));
    connect(PHX.T_out_pFluid, pipe.PiTemp_IN) annotation(
      Line(points = {{2382, 810}, {3028, 810}, {3028, -820}, {167, -820}, {167, -755}}));
    connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
      Line(points = {{-626, -1506}, {168, -1506}, {168, -956}}, color = {245, 121, 0}));
    connect(pipe.PiTemp_Out, nineRcoreNoDH.temp_In) annotation(
      Line(points = {{-234, -754}, {-1084, -754}, {-1084, -510}}));
    connect(nineRcoreNoDH.temp_Out, dhrs.DHRS_tempIN) annotation(
      Line(points = {{-1096, 76}, {-1072, 76}, {-1072, 714}, {-488, 714}}));
    connect(nineRcoreNoDH.decayHeat_Out, dhrs.DHRS_DecayHeat) annotation(
      Line(points = {{-1614, 76}, {-1600, 76}, {-1600, 520}, {-488, 520}}, color = {0, 225, 255}));
    connect(nineRcoreNoDH.decayHeat_Out, pipe.PiDecay_Heat) annotation(
      Line(points = {{-1614, 76}, {-1600, 76}, {-1600, 284}, {168, 284}, {168, -576}}, color = {0, 225, 255}));
    connect(nineRcoreNoDH.decayHeat_Out, PHX.P_decay) annotation(
      Line(points = {{-1614, 76}, {-1600, 76}, {-1600, 1106}, {1810, 1106}, {1810, 894}}, color = {0, 225, 255}));
    connect(primaryPump.primaryFlowFrac, nineRcoreNoDH.flowFraction_In) annotation(
      Line(points = {{-626, -1506}, {-1626, -1506}, {-1626, -518}}, color = {245, 121, 0}));
    annotation(
      Diagram(coordinateSystem(extent = {{-1480, 1160}, {3040, -1780}})));
  end nineE_MSREnoDH;
end Transient4;