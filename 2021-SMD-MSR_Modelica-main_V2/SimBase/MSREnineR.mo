package MSREnineR
  model nineRcore
    parameter SMD_MSR_Modelica.Units.ReactorPower P;
    parameter SMD_MSR_Modelica.Units.Mass TotalFuelMass;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_fuel;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_grap;
    parameter SMD_MSR_Modelica.Units.Convection hA[9];
    parameter SMD_MSR_Modelica.Units.VolumeImportance kFN1[9];
    parameter SMD_MSR_Modelica.Units.VolumeImportance kFN2[9];
    parameter SMD_MSR_Modelica.Units.HeatTransferFraction kHT1[9];
    parameter SMD_MSR_Modelica.Units.HeatTransferFraction kHT2[9];
    parameter SMD_MSR_Modelica.Units.Mass mMix;
    parameter SMD_MSR_Modelica.Units.Mass mFN1[9];
    parameter SMD_MSR_Modelica.Units.Mass mFN2[9];
    parameter SMD_MSR_Modelica.Units.Mass mGN[9];
    parameter SMD_MSR_Modelica.Units.MassFlowRate mDot;
    parameter SMD_MSR_Modelica.Units.Temperature TF1_0[9];
    parameter SMD_MSR_Modelica.Units.Temperature TF2_0[9];
    parameter SMD_MSR_Modelica.Units.Temperature TG_0[9];
    parameter SMD_MSR_Modelica.Units.Temperature Tmix_0;
    parameter Real IF1[9];
    parameter Real IF2[9];
    parameter Real IG[9];
    parameter SMD_MSR_Modelica.Units.TemperatureReactivityCoef aF;
    parameter SMD_MSR_Modelica.Units.TemperatureReactivityCoef aG;
    parameter SMD_MSR_Modelica.Units.Reactivity rho_0;
    parameter SMD_MSR_Modelica.Units.Reactivity stepReact;
    parameter SMD_MSR_Modelica.Units.InitiationTime stepTime;
    parameter SMD_MSR_Modelica.Units.Reactivity sinReact;
    parameter SMD_MSR_Modelica.Units.Frequency omega;
    parameter SMD_MSR_Modelica.Units.InitiationTime sinTime;
    parameter SMD_MSR_Modelica.Units.NeutronGenerationTime Lam;
    parameter SMD_MSR_Modelica.Units.PrecursorConc beta[6];
    parameter SMD_MSR_Modelica.Units.PrecursorDecayConstant lam[6];
    parameter SMD_MSR_Modelica.Units.ResidentTime tauCore;
    parameter SMD_MSR_Modelica.Units.ResidentTime tauLoop;
    parameter SMD_MSR_Modelica.Units.ResidentTime tauCoreDHRS;
    parameter SMD_MSR_Modelica.Units.FlowFraction flowFracRegions[4];
    parameter SMD_MSR_Modelica.Units.InitiationTime regionTripTime[4];
    parameter SMD_MSR_Modelica.Units.PumpConstant regionCoastDown;
    parameter SMD_MSR_Modelica.Units.FlowFraction regionFreeConvFF;
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = P, TotalFuelMass = TotalFuelMass) annotation(
      Placement(visible = true, transformation(origin = {-298, 228}, extent = {{-82, -82}, {82, 82}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region1(P = P, TF1_0 = TF1_0[1], TF2_0 = TF2_0[1], TG_0 = TG_0[1], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[1], kFN1 = kFN1[1], kFN2 = kFN2[1], kG = kHT1[1] + kHT2[1], kHT_FN1 = kHT1[1], kHT_FN2 = kHT2[1], m_FN1 = mFN1[1], m_FN2 = mFN2[1], m_GN = mGN[1], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[1], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {-288, 60}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region2(P = P, TF1_0 = TF1_0[2], TF2_0 = TF2_0[2], TG_0 = TG_0[2], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[2], kFN1 = kFN1[2], kFN2 = kFN2[2], kG = kHT1[2] + kHT2[2], kHT_FN1 = kHT1[2], kHT_FN2 = kHT2[2], m_FN1 = mFN1[2], m_FN2 = mFN2[2], m_GN = mGN[2], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[2], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {-26, -144}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region3(P = P, TF1_0 = TF1_0[3], TF2_0 = TF2_0[3], TG_0 = TG_0[3], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[3], kFN1 = kFN1[3], kFN2 = kFN2[3], kG = kHT1[3] + kHT2[3], kHT_FN1 = kHT1[3], kHT_FN2 = kHT2[3], m_FN1 = mFN1[3], m_FN2 = mFN2[3], m_GN = mGN[3], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[2], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {-28, 58}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region4(P = P, TF1_0 = TF1_0[4], TF2_0 = TF2_0[4], TG_0 = TG_0[4], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[4], kFN1 = kFN1[4], kFN2 = kFN2[4], kG = kHT1[4] + kHT2[4], kHT_FN1 = kHT1[4], kHT_FN2 = kHT2[4], m_FN1 = mFN1[4], m_FN2 = mFN2[4], m_GN = mGN[4], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[2], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {-27, 235}, extent = {{-57, -57}, {57, 57}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region5(P = P, TF1_0 = TF1_0[5], TF2_0 = TF2_0[5], TG_0 = TG_0[5], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[5], kFN1 = kFN1[5], kFN2 = kFN2[5], kG = kHT1[5] + kHT2[5], kHT_FN1 = kHT1[5], kHT_FN2 = kHT2[5], m_FN1 = mFN1[5], m_FN2 = mFN2[5], m_GN = mGN[5], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[3], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {208, -146}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region6(P = P, TF1_0 = TF1_0[6], TF2_0 = TF2_0[6], TG_0 = TG_0[6], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[6], kFN1 = kFN1[6], kFN2 = kFN2[6], kG = kHT1[6] + kHT2[6], kHT_FN1 = kHT1[6], kHT_FN2 = kHT2[6], m_FN1 = mFN1[6], m_FN2 = mFN2[6], m_GN = mGN[6], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[3], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {212, 54}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region7(P = P, TF1_0 = TF1_0[7], TF2_0 = TF2_0[7], TG_0 = TG_0[7], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[7], kFN1 = kFN1[7], kFN2 = kFN2[7], kG = kHT1[7] + kHT2[7], kHT_FN1 = kHT1[7], kHT_FN2 = kHT2[7], m_FN1 = mFN1[7], m_FN2 = mFN2[7], m_GN = mGN[7], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[3], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {212, 234}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region8(P = P, TF1_0 = TF1_0[8], TF2_0 = TF2_0[8], TG_0 = TG_0[8], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[8], kFN1 = kFN1[8], kFN2 = kFN2[8], kG = kHT1[8] + kHT2[8], kHT_FN1 = kHT1[8], kHT_FN2 = kHT2[8], m_FN1 = mFN1[8], m_FN2 = mFN2[8], m_GN = mGN[8], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[4], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {472, -18}, extent = {{-58, -58}, {58, 58}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region9(P = P, TF1_0 = TF1_0[9], TF2_0 = TF2_0[9], TG_0 = TG_0[9], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[9], kFN1 = kFN1[9], kFN2 = kFN2[9], kG = kHT1[9] + kHT2[9], kHT_FN1 = kHT1[9], kHT_FN2 = kHT2[9], m_FN1 = mFN1[9], m_FN2 = mFN2[9], m_GN = mGN[9], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[4], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {467, 159}, extent = {{-57, -57}, {57, 57}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region1FB(FuelTempSetPointNode1 = TF1_0[1], FuelTempSetPointNode2 = TF2_0[1], GrapTempSetPoint = TG_0[1], IF1 = IF1[1], IF2 = IF2[1], IG = IG[1], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {-188, 46}, extent = {{-56, 56}, {56, -56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region2FB(FuelTempSetPointNode1 = TF1_0[2], FuelTempSetPointNode2 = TF2_0[2], GrapTempSetPoint = TG_0[2], IF1 = IF1[2], IF2 = IF2[2], IG = IG[2], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {73, -155}, extent = {{-55, 55}, {55, -55}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region3FB(FuelTempSetPointNode1 = TF1_0[3], FuelTempSetPointNode2 = TF2_0[3], GrapTempSetPoint = TG_0[3], IF1 = IF1[3], IF2 = IF2[3], IG = IG[3], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {74, 40}, extent = {{-56, 56}, {56, -56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region4FB(FuelTempSetPointNode1 = TF1_0[4], FuelTempSetPointNode2 = TF2_0[4], GrapTempSetPoint = TG_0[4], IF1 = IF1[4], IF2 = IF2[4], IG = IG[4], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {76, 220}, extent = {{-54, 54}, {54, -54}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region5FB(FuelTempSetPointNode1 = TF1_0[5], FuelTempSetPointNode2 = TF2_0[5], GrapTempSetPoint = TG_0[5], IF1 = IF1[5], IF2 = IF2[5], IG = IG[5], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {302, -154}, extent = {{-56, 56}, {56, -56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region6FB(FuelTempSetPointNode1 = TF1_0[6], FuelTempSetPointNode2 = TF2_0[6], GrapTempSetPoint = TG_0[6], IF1 = IF1[6], IF2 = IF2[6], IG = IG[6], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {294, 34}, extent = {{-54, 54}, {54, -54}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region7FB(FuelTempSetPointNode1 = TF1_0[7], FuelTempSetPointNode2 = TF2_0[7], GrapTempSetPoint = TG_0[7], IF1 = IF1[7], IF2 = IF2[7], IG = IG[7], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {304, 222}, extent = {{-54, 54}, {54, -54}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region8FB(FuelTempSetPointNode1 = TF1_0[8], FuelTempSetPointNode2 = TF2_0[8], GrapTempSetPoint = TG_0[8], IF1 = IF1[8], IF2 = IF2[8], IG = IG[8], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {606, -34}, extent = {{-58, 58}, {58, -58}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region9FB(FuelTempSetPointNode1 = TF1_0[9], FuelTempSetPointNode2 = TF2_0[9], GrapTempSetPoint = TG_0[9], IF1 = IF1[9], IF2 = IF2[9], IG = IG[9], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {586, 146}, extent = {{-56, 56}, {56, -56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(P = P, DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5, TotalFuelMass = TotalFuelMass) annotation(
      Placement(visible = true, transformation(origin = {-431, -105}, extent = {{-79, -79}, {79, 79}}, rotation = 180)));
    SMD_MSR_Modelica.HeatTransport.FlowDistributor flowDistributor(coastDownK = 1, freeConvectionFF = regionFreeConvFF, numOutput = 4, regionTripTime = regionTripTime) annotation(
      Placement(visible = true, transformation(origin = {11, -455}, extent = {{-53, -53}, {53, 53}}, rotation = 180)));
    SMD_MSR_Modelica.HeatTransport.MixingPot mixingPot(Cp = cP_fuel, T_0 = Tmix_0, m = mMix, mDotNom = mDot, numInput = 4, regionFlowFrac = flowFracRegions, tauPotToDHRS = tauCoreDHRS) annotation(
      Placement(visible = true, transformation(origin = {80, 442}, extent = {{-110, -110}, {110, 110}}, rotation = 0)));
    input SMD_MSR_Modelica.PortsConnectors.Temp_In temp_In annotation(
      Placement(visible = true, transformation(origin = {180, -400}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {68, -74}, extent = {{-28, -28}, {28, 28}}, rotation = 0)));
    input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In flowFraction_In annotation(
      Placement(visible = true, transformation(origin = {80, -434}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-68, -76}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
    output SMD_MSR_Modelica.PortsConnectors.Temp_Out temp_Out annotation(
      Placement(visible = true, transformation(origin = {68, 510}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {65, 73}, extent = {{-25, -25}, {25, 25}}, rotation = 0)));
    output SMD_MSR_Modelica.PortsConnectors.DecayHeat_Out decayHeat_Out annotation(
      Placement(visible = true, transformation(origin = {10, 510}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-65, 73}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.SumReactivity sumReactivity(numInput = 9) annotation(
      Placement(visible = true, transformation(origin = {-259, -379}, extent = {{-105, -105}, {105, 105}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.mPKE mPKE(Lam = Lam, beta = beta, lam = lam, omega = omega, rho_0 = rho_0, sinInsertionTime = sinTime, sin_mag = sinReact, stepInsertionTime = stepTime, step_mag = stepReact, tauCoreNom = tauCore, tauLoopNom = tauLoop) annotation(
      Placement(visible = true, transformation(origin = {-432, -274}, extent = {{-116, -116}, {116, 116}}, rotation = 180)));
  initial equation
  //temp_In.T = 6.3222E+02;
  equation
    connect(Region2.temp_Out, Region3.temp_In) annotation(
      Line(points = {{2, -114}, {12, -114}, {12, -18}, {-73, -18}, {-73, 19}}));
    connect(Region3.temp_Out, Region4.temp_In) annotation(
      Line(points = {{0, 88}, {12, 88}, {12, 106}, {-73, 106}, {-73, 195}}));
    connect(Region5.temp_Out, Region6.temp_In) annotation(
      Line(points = {{236, -116}, {236, -26}, {167, -26}, {167, 15}}));
    connect(Region6.temp_Out, Region7.temp_In) annotation(
      Line(points = {{240, 84}, {240, 106}, {167, 106}, {167, 195}}));
    connect(Region8.temp_Out, Region9.temp_In) annotation(
      Line(points = {{501, 13}, {486, 13}, {486, 54}, {421, 54}, {421, 119}}));
    connect(Region1.fuelNode1, Region1FB.fuelNode1) annotation(
      Line(points = {{-260, 20}, {-210, 20}, {-210, 27}}));
    connect(Region1.fuelNode2, Region1FB.fuelNode2) annotation(
      Line(points = {{-260, 43}, {-234, 43}, {-234, 46}, {-210, 46}}));
    connect(Region1.grapNode, Region1FB.grapNode) annotation(
      Line(points = {{-260, 69}, {-236, 69}, {-236, 63}, {-210, 63}}));
    connect(Region2.fuelNode1, Region2FB.fuelNode1) annotation(
      Line(points = {{2, -184}, {30.5, -184}, {30.5, -174}, {51, -174}}));
    connect(Region2.fuelNode2, Region2FB.fuelNode2) annotation(
      Line(points = {{2, -161}, {29.5, -161}, {29.5, -155}, {51, -155}}));
    connect(Region2.grapNode, Region2FB.grapNode) annotation(
      Line(points = {{2, -135}, {20.5, -135}, {20.5, -138.5}, {51, -138.5}}));
    connect(Region3.fuelNode1, Region3FB.fuelNode1) annotation(
      Line(points = {{0, 18}, {52, 18}, {52, 21}}));
    connect(Region3.fuelNode2, Region3FB.fuelNode2) annotation(
      Line(points = {{0, 41}, {26, 41}, {26, 40}, {52, 40}}));
    connect(Region3.grapNode, Region3FB.grapNode) annotation(
      Line(points = {{0, 67}, {52, 67}, {52, 57}}));
    connect(Region4.fuelNode1, Region4FB.fuelNode1) annotation(
      Line(points = {{2, 194}, {29, 194}, {29, 202}, {54, 202}}));
    connect(Region4.fuelNode2, Region4FB.fuelNode2) annotation(
      Line(points = {{2, 218}, {33, 218}, {33, 220}, {54, 220}}));
    connect(Region4.grapNode, Region4FB.grapNode) annotation(
      Line(points = {{2, 244}, {28, 244}, {28, 236}, {54, 236}}));
    connect(Region5.fuelNode1, Region5FB.fuelNode1) annotation(
      Line(points = {{236, -186}, {280, -186}, {280, -173}}));
    connect(Region5.fuelNode2, Region5FB.fuelNode2) annotation(
      Line(points = {{236, -163}, {280, -163}, {280, -154}}));
    connect(Region5.grapNode, Region5FB.grapNode) annotation(
      Line(points = {{236, -137}, {280, -137}}));
    connect(Region6.fuelNode1, Region6FB.fuelNode1) annotation(
      Line(points = {{240, 14}, {272, 14}, {272, 16}}));
    connect(Region6.fuelNode2, Region6FB.fuelNode2) annotation(
      Line(points = {{240, 38}, {256, 38}, {256, 34}, {272, 34}}));
    connect(Region6.grapNode, Region6FB.grapNode) annotation(
      Line(points = {{240, 62}, {255, 62}, {255, 50}, {272, 50}}));
    connect(Region7.fuelNode1, Region7FB.fuelNode1) annotation(
      Line(points = {{240, 194}, {282, 194}, {282, 204}}));
    connect(Region7.fuelNode2, Region7FB.fuelNode2) annotation(
      Line(points = {{240, 217}, {282, 217}, {282, 222}}));
    connect(Region7.grapNode, Region7FB.grapNode) annotation(
      Line(points = {{240, 243}, {282, 243}, {282, 238}}));
    connect(Region8.fuelNode1, Region8FB.fuelNode1) annotation(
      Line(points = {{501, -60}, {542, -60}, {542, -54}, {583, -54}}));
    connect(Region8.fuelNode2, Region8FB.fuelNode2) annotation(
      Line(points = {{501, -35}, {508.5, -35}, {508.5, -34}, {583, -34}}));
    connect(Region8.grapNode, Region8FB.grapNode) annotation(
      Line(points = {{501, -9}, {539, -9}, {539, -17}, {583, -17}}));
    connect(Region9.fuelNode1, Region9FB.fuelNode1) annotation(
      Line(points = {{495.5, 118}, {507.75, 118}, {507.75, 127}, {564, 127}}));
    connect(Region9.fuelNode2, Region9FB.fuelNode2) annotation(
      Line(points = {{495.5, 142}, {530.75, 142}, {530.75, 146}, {564, 146}}));
    connect(Region9.grapNode, Region9FB.grapNode) annotation(
      Line(points = {{495.5, 168}, {515.75, 168}, {515.75, 163}, {564, 163}}));
    connect(temp_In, Region1.temp_In) annotation(
      Line(points = {{180, -400}, {-333, -400}, {-333, 21}}, color = {28, 113, 216}));
    connect(temp_In, Region2.temp_In) annotation(
      Line(points = {{180, -400}, {-71, -400}, {-71, -183}}, color = {28, 113, 216}));
    connect(temp_In, Region5.temp_In) annotation(
      Line(points = {{180, -400}, {180, -269}, {163, -269}, {163, -185}}, color = {28, 113, 216}));
    connect(temp_In, Region8.temp_In) annotation(
      Line(points = {{180, -400}, {426, -400}, {426, -59}}, color = {28, 113, 216}));
    connect(decayHeat.decayHeat_Out, Region1.P_decay) annotation(
      Line(points = {{-459, -81}, {-384, -81}, {-384, 57}, {-333, 57}}, color = {0, 225, 255}));
    connect(decayHeat.decayHeat_Out, Region2.P_decay) annotation(
      Line(points = {{-459, -81}, {-136, -81}, {-136, -147}, {-71, -147}}, color = {0, 225, 255}));
    connect(decayHeat.decayHeat_Out, Region3.P_decay) annotation(
      Line(points = {{-459, -81}, {-90, -81}, {-90, 55}, {-73, 55}}, color = {0, 225, 255}));
    connect(decayHeat.decayHeat_Out, Region4.P_decay) annotation(
      Line(points = {{-459, -81}, {-90, -81}, {-90, 232}, {-72, 232}}, color = {0, 225, 255}));
    connect(decayHeat.decayHeat_Out, Region5.P_decay) annotation(
      Line(points = {{-459, -81}, {152, -81}, {152, -149}, {163, -149}}, color = {0, 225, 255}));
    connect(decayHeat.decayHeat_Out, Region6.P_decay) annotation(
      Line(points = {{-459, -81}, {152, -81}, {152, 50}, {168, 50}}, color = {0, 225, 255}));
    connect(decayHeat.decayHeat_Out, Region7.P_decay) annotation(
      Line(points = {{-459, -81}, {152, -81}, {152, 230}, {168, 230}}, color = {0, 225, 255}));
    connect(decayHeat.decayHeat_Out, Region8.P_decay) annotation(
      Line(points = {{-459, -81}, {364, -81}, {364, -21}, {426, -21}}, color = {0, 225, 255}));
    connect(decayHeat.decayHeat_Out, Region9.P_decay) annotation(
      Line(points = {{-459, -81}, {364, -81}, {364, 156}, {421, 156}}, color = {0, 225, 255}));
    connect(decayHeat_Out, decayHeat.decayHeat_Out) annotation(
      Line(points = {{10, 510}, {-466, 510}, {-466, -81}, {-459, -81}}, color = {0, 225, 255}));
    connect(decayHeat.decayHeat_Out, powerBlock.P_decay) annotation(
      Line(points = {{-459, -81}, {-132, -81}, {-132, 196}, {-272, 196}}, color = {0, 225, 255}));
    connect(flowFraction_In, flowDistributor.flowFracIn) annotation(
      Line(points = {{80, -434}, {33, -434}}));
    connect(flowDistributor.flowFracOut[1], Region1.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-310, -434}, {-310, 15}}));
    connect(flowDistributor.flowFracOut[2], Region2.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -320.5}, {-48, -320.5}, {-48, -189}}));
    connect(flowDistributor.flowFracOut[2], Region3.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -220}, {-50, -220}, {-50, 13}}));
    connect(flowDistributor.flowFracOut[2], Region4.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -131}, {-50, -131}, {-50, 190}}));
    connect(flowDistributor.flowFracOut[3], Region5.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -264}, {186, -264}, {186, -191}}));
    connect(flowDistributor.flowFracOut[3], Region7.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -264}, {190, -264}, {190, 190}}));
    connect(flowDistributor.flowFracOut[3], Region6.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -266}, {190, -266}, {190, 10}}));
    connect(flowDistributor.flowFracOut[4], Region8.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -266}, {449, -266}, {449, -64}}));
    connect(flowDistributor.flowFracOut[4], Region9.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -268}, {444, -268}, {444, 113}}));
    connect(Region1.temp_Out, mixingPot.temp_In[1]) annotation(
      Line(points = {{-260, 90}, {-192, 90}, {-192, 440}, {14, 440}}));
    connect(Region4.temp_Out, mixingPot.temp_In[2]) annotation(
      Line(points = {{2, 266}, {2, 357}, {14, 357}, {14, 440}}));
    connect(Region7.temp_Out, mixingPot.temp_In[3]) annotation(
      Line(points = {{240, 264}, {286, 264}, {286, 314}, {14, 314}, {14, 440}}));
    connect(Region9.temp_Out, mixingPot.temp_In[4]) annotation(
      Line(points = {{495.5, 189}, {472, 189}, {472, 316}, {14, 316}, {14, 440}}));
    connect(mixingPot.temp_Out, temp_Out) annotation(
      Line(points = {{146, 442}, {142, 442}, {142, 510}, {68, 510}}));
    connect(decayHeat.decayHeat_Out, mixingPot.decayHeat) annotation(
      Line(points = {{-459, -81}, {-440, -81}, {-440, 375}, {145, 375}}, color = {0, 225, 255}));
    connect(flowDistributor.flowFracOut[1], mixingPot.flowFraction[1]) annotation(
      Line(points = {{-10, -434}, {-538, -434}, {-538, 376}, {16, 376}}, color = {245, 121, 0}, thickness = 0.5));
    connect(flowDistributor.flowFracOut[2], mixingPot.flowFraction[2]) annotation(
      Line(points = {{-10, -434}, {-564, -434}, {-564, 376}, {16, 376}}, color = {245, 121, 0}, thickness = 0.5));
    connect(flowDistributor.flowFracOut[3], mixingPot.flowFraction[3]) annotation(
      Line(points = {{-10, -434}, {-552, -434}, {-552, 376}, {16, 376}}, color = {245, 121, 0}, thickness = 0.5));
    connect(flowDistributor.flowFracOut[4], mixingPot.flowFraction[4]) annotation(
      Line(points = {{-10, -434}, {-528, -434}, {-528, 376}, {16, 376}}, color = {245, 121, 0}, thickness = 0.5));
    connect(Region1FB.feedback, sumReactivity.reactivityIn[1]) annotation(
      Line(points = {{-168, 36}, {-127, 36}, {-127, -322}, {-126, -322}, {-126, -358}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region2FB.feedback, sumReactivity.reactivityIn[2]) annotation(
      Line(points = {{92, -165}, {-238, -165}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region3FB.feedback, sumReactivity.reactivityIn[3]) annotation(
      Line(points = {{93, 30}, {39, 30}, {39, -358}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region4FB.feedback, sumReactivity.reactivityIn[4]) annotation(
      Line(points = {{94, 210}, {-238, 210}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region5FB.feedback, sumReactivity.reactivityIn[5]) annotation(
      Line(points = {{322, -164}, {316, -164}, {316, -358}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region6FB.feedback, sumReactivity.reactivityIn[6]) annotation(
      Line(points = {{312, 24}, {352, 24}, {352, -358}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region7FB.feedback, sumReactivity.reactivityIn[7]) annotation(
      Line(points = {{322, 212}, {-238, 212}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region8FB.feedback, sumReactivity.reactivityIn[8]) annotation(
      Line(points = {{626, -44}, {556, -44}, {556, -358}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region9FB.feedback, sumReactivity.reactivityIn[9]) annotation(
      Line(points = {{605, 136}, {664.5, 136}, {664.5, 44}, {663, 44}, {663, -358}, {-238, -358}}, color = {78, 154, 6}));
    connect(sumReactivity.reactivityOut, mPKE.feedback) annotation(
      Line(points = {{-280, -358}, {-432, -358}, {-432, -320}}, color = {78, 154, 6}));
    connect(flowFraction_In, mPKE.fuelFlowFrac) annotation(
      Line(points = {{80, -434}, {80, -480}, {-471, -480}, {-471, -320}}, color = {255, 120, 0}));
    connect(mPKE.n_population, decayHeat.nPop) annotation(
      Line(points = {{-432, -230}, {-406, -230}, {-406, -122}}, color = {20, 36, 248}));
    connect(mPKE.n_population, powerBlock.nPop) annotation(
      Line(points = {{-432, -230}, {-398, -230}, {-398, 196}, {-340, 196}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region1.nPop) annotation(
      Line(points = {{-432, -230}, {-368, -230}, {-368, 92}, {-333, 92}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region2.nPop) annotation(
      Line(points = {{-432, -230}, {-156, -230}, {-156, -112}, {-70, -112}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region3.nPop) annotation(
      Line(points = {{-432, -230}, {-116, -230}, {-116, 90}, {-73, 90}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region4.nPop) annotation(
      Line(points = {{-432, -230}, {-122, -230}, {-122, 268}, {-72, 268}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region5.nPop) annotation(
      Line(points = {{-432, -230}, {130, -230}, {130, -114}, {163, -114}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region6.nPop) annotation(
      Line(points = {{-432, -230}, {136, -230}, {136, 86}, {168, 86}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region7.nPop) annotation(
      Line(points = {{-432, -230}, {138, -230}, {138, 266}, {168, 266}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region8.nPop) annotation(
      Line(points = {{-432, -230}, {396, -230}, {396, 16}, {426, 16}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region9.nPop) annotation(
      Line(points = {{-432, -230}, {404, -230}, {404, 192}, {421, 192}}, color = {20, 36, 248}));
    annotation(
      Diagram(coordinateSystem(extent = {{-560, 580}, {600, -520}}), graphics = {Text(origin = {-360, 106}, extent = {{-26, 10}, {26, -10}}, textString = "Region1"), Text(origin = {-110, -96}, extent = {{-26, 10}, {26, -10}}, textString = "Region2"), Text(origin = {-112, 104}, extent = {{-26, 10}, {26, -10}}, textString = "Region3"), Text(origin = {-110, 286}, extent = {{-26, 10}, {26, -10}}, textString = "Region4"), Text(origin = {132, -96}, extent = {{-26, 10}, {26, -10}}, textString = "Region5"), Text(origin = {132, 104}, extent = {{-26, 10}, {26, -10}}, textString = "Region6"), Text(origin = {136, 286}, extent = {{-26, 10}, {26, -10}}, textString = "Region7"), Text(origin = {366, 32}, extent = {{-26, 10}, {26, -10}}, textString = "Region8"), Text(origin = {366, 196}, extent = {{-26, 10}, {26, -10}}, textString = "Region9")}),
      Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Text(origin = {-1, -2}, extent = {{-71, 28}, {71, -28}}, textString = "9Rcore")}));
  end nineRcore;

  model nineE_MSRE
    parameter Real CpCoefFuel = 1;
    parameter Real HTCcoefCore = 1;
    parameter Real HTCcoefPHX = 1;
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.05, primaryPumpK = 0.2, tripPrimaryPump = 2000000) annotation(
      Placement(visible = true, transformation(origin = {-626, -1274}, extent = {{-386, -386}, {386, 386}}, rotation = 0)));
    MSREnineR.nineRcore nineRcore(IF1 = {0.02168, 0.02197, 0.07897, 0.08249, 0.02254, 0.08255, 0.08623, 0.02745, 0.06936}, IF2 = {0.02678, 0.06519, 0.08438, 0.04124, 0.06801, 0.08823, 0.04290, 0.05529, 0.03473}, IG = {0.04443, 0.08835, 0.16671, 0.12077, 0.09181, 0.17429, 0.12612, 0.08408, 0.10343}, Lam = 2.400E-04, P = 8, TF1_0 = {641.0644, 636.1926, 662.2593, 692.4323, 634.6308, 650.7339, 669.3933, 633.4354, 638.9545}, TF2_0 = {651.9893, 647.9801, 677.5167, 699.8892, 641.9048, 660.1706, 673.9817, 635.8835, 640.4922}, TG_0 = {664.0462, 653.7083, 692.6790, 714.5485, 648.7869, 674.7077, 686.8285, 640.8657, 647.7290}, Tmix_0 = 663.5300, TotalFuelMass = 4.093499709644400e+03, aF = -8.71E-05, aG = -6.66E-05, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, cP_fuel = 0.0020*CpCoefFuel, cP_grap = 0.0018, flowFracRegions = {0.061409895150439, 0.138549617100260, 0.234231239311780, 0.565809248437521}, hA = {0.0007056, 0.0021672, 0.00162, 0.0021132, 0.0035586, 0.002745, 0.003573, 0.009801, 0.009648}, kFN1 = {0.01493, 0.02736, 0.04504, 0.05126, 0.03601, 0.06014, 0.06845, 0.06179, 0.09333}, kFN2 = {0.01721, 0.04550, 0.04656, 0.04261, 0.06069, 0.06218, 0.05664, 0.07707, 0.07311}, kHT1 = {0.000946, 0.001685, 0.003029, 0.003447, 0.002216, 0.004044, 0.004603, 0.003920, 0.006277}, kHT2 = {0.001081, 0.003060, 0.003131, 0.002395, 0.004081, 0.004182, 0.003184, 0.005183, 0.004305}, lam = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, mDot = 162.3880, mFN1 = {13.8215, 46.8650, 25.6293, 32.0366, 79.2677, 43.2952, 54.1876, 217.8490, 147.8261}, mFN2 = {14.4622, 31.9451, 25.6293, 62.4256, 54.1876, 43.2952, 105.4462, 126.6819, 248.0549}, mGN = {71.0660, 214.6193, 163.0457, 208.7310, 363.0457, 275.9391, 353.0964, 975.8376, 956.4467}, mMix = 324.7760, omega = 1, regionCoastDown = 0.02, regionFreeConvFF = 0.01, regionTripTime = {200000, 200000, 200000, 200000}, rho_0 = 0.002465140767843, sinReact = 0, sinTime = 20000000, stepReact = 100E-5, stepTime = 2000, tauCore = 8.46, tauCoreDHRS = 1.3850, tauLoop = 16.73) annotation(
      Placement(visible = true, transformation(origin = {-1349, -215}, extent = {{-399, -399}, {399, 399}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.05, secondaryPumpK = 0.2, tripSecondaryPump = 200000) annotation(
      Placement(visible = true, transformation(origin = {1204, -1400}, extent = {{-362, -362}, {362, 362}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.UHX uhx(SetDemand = 1, UHXP = 8, cP_pFluidUHX = 2.39E-3, m_PN1UHX = 3.924401289158446e+02, m_dot_pFluidUHXnom = 100.5793, tauUHXtoPHX = 8.24, tripUHX = 2000000) annotation(
      Placement(visible = true, transformation(origin = {2441, -231}, extent = {{-567, 567}, {567, -567}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.PrimaryHeatExchanger PHX(T_PN1_initial = 6.478750000000000e+02, cP_Tube = 5.778E-04, cP_pFluid = 0.0020*CpCoefFuel, cP_sFluid = 2.39E-3, hApnNom = 0.1620*HTCcoefPHX, hAsnNom = 0.0765, m_PN1 = 87.0959, m_PN2 = 87.0959, m_PN3 = 87.0959, m_PN4 = 87.0959, m_SN1 = 24.9200, m_SN2 = 24.9200, m_SN3 = 24.9200, m_SN4 = 24.9200, m_TN1 = 101.1662, m_TN2 = 101.1662, m_dot_pFluidNom = 162.3880, m_dot_sFluidNom = 100.5793, tauPHXtoPipe = 3.8350, tauPHXtoUHX = 4.71) annotation(
      Placement(visible = true, transformation(origin = {1691, 681}, extent = {{-427, -427}, {427, 427}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Cp_fluid = 0.0020*CpCoefFuel, DHRS_MaxPowerRm = 8, DHRS_PowerBleed = 0, DHRS_time = 20000, DHRS_timeConstant = 0.02, m_DHRS = 162.3880, m_Pi_C_PHX = 162.3880*3.7700, m_dotNom = 162.3880, tauDHRStoPHX = 1.3850)  annotation(
      Placement(visible = true, transformation(origin = {-183, 825}, extent = {{509, -509}, {-509, 509}}, rotation = 180)));
  SMD_MSR_Modelica.HeatTransport.Pipe pipe(Cp_fluidPi = 0.0020*CpCoefFuel, m_Pi_PHX_C = 162.3880*8.6700, m_dotPiNom = 162.3880, m_pi = 162.3880, tauPiToCore = 3.8350)  annotation(
      Placement(visible = true, transformation(origin = {-86, -892}, extent = {{528, -528}, {-528, 528}}, rotation = 0)));
  equation
    connect(primaryPump.primaryFlowFrac, nineRcore.flowFraction_In) annotation(
      Line(points = {{-626, -1506}, {-1620, -1506}, {-1620, -518}}, color = {245, 121, 0}));
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
      Line(points = {{-1608, 76}, {-1582, 76}, {-1582, 1104}, {1810, 1104}, {1810, 894}}, color = {0, 225, 255}));
  connect(dhrs.DHRS_TempOUT, PHX.T_in_pFluid) annotation(
      Line(points = {{40, 714}, {1256, 714}, {1256, 810}}));
  connect(nineRcore.decayHeat_Out, dhrs.DHRS_DecayHeat) annotation(
      Line(points = {{-1608, 76}, {-1594, 76}, {-1594, 520}, {-488, 520}}, color = {0, 225, 255}));
  connect(nineRcore.temp_Out, dhrs.DHRS_tempIN) annotation(
      Line(points = {{-1090, 76}, {-1026, 76}, {-1026, 714}, {-488, 714}}));
  connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
      Line(points = {{-626, -1506}, {-722, -1506}, {-722, 926}, {-488, 926}}, color = {245, 121, 0}));
  connect(nineRcore.decayHeat_Out, pipe.PiDecay_Heat) annotation(
      Line(points = {{-1608, 76}, {-1414, 76}, {-1414, 142}, {168, 142}, {168, -575}, {167, -575}}, color = {0, 225, 255}));
  connect(PHX.T_out_pFluid, pipe.PiTemp_IN) annotation(
      Line(points = {{2382, 810}, {3028, 810}, {3028, -820}, {167, -820}, {167, -755}}));
  connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
      Line(points = {{-626, -1506}, {168, -1506}, {168, -956}}, color = {245, 121, 0}));
  connect(pipe.PiTemp_Out, nineRcore.temp_In) annotation(
      Line(points = {{-234, -754}, {-1078, -754}, {-1078, -510}}));
    annotation(
      Diagram(coordinateSystem(extent = {{-1480, 1160}, {3040, -1780}})));
  end nineE_MSRE;

  model nineEcoreTest
    parameter Real CpCoefFuel = 1;
    parameter Real HTCcoefCore = 1;
    parameter Real HTCcoefPHX = 1;
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.05, primaryPumpK = 0.2, tripPrimaryPump = 200000) annotation(
      Placement(visible = true, transformation(origin = {-626, -1274}, extent = {{-386, -386}, {386, 386}}, rotation = 0)));
    MSREnineR.nineRcore nineRcore(IF1 = {0.02168, 0.02197, 0.07897, 0.08249, 0.02254, 0.08255, 0.08623, 0.02745, 0.06936}, IF2 = {0.02678, 0.06519, 0.08438, 0.04124, 0.06801, 0.08823, 0.04290, 0.05529, 0.03473}, IG = {0.04443, 0.08835, 0.16671, 0.12077, 0.09181, 0.17429, 0.12612, 0.08408, 0.10343}, Lam = 2.400E-04, P = 8, TF1_0 = {641.0644, 636.1926, 662.2593, 692.4323, 634.6308, 650.7339, 669.3933, 633.4354, 638.9545}, TF2_0 = {651.9893, 647.9801, 677.5167, 699.8892, 641.9048, 660.1706, 673.9817, 635.8835, 640.4922}, TG_0 = {664.0462, 653.7083, 692.6790, 714.5485, 648.7869, 674.7077, 686.8285, 640.8657, 647.7290}, Tmix_0 = 663.5300, TotalFuelMass = 4.093499709644400e+03, aF = -8.71E-05, aG = -6.66E-05, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, cP_fuel = 1.9665E-3*CpCoefFuel, cP_grap = 1.773E-3, flowFracRegions = {0.0614, 0.1385, 0.2342, 0.5658}, hA = {0.0007, 0.0022, 0.0016, 0.0021, 0.0036, 0.0027, 0.0036, 0.0098, 0.0096}, kFN1 = {0.01493, 0.02736, 0.04504, 0.05126, 0.03601, 0.06014, 0.06845, 0.06179, 0.09333}, kFN2 = {0.01721, 0.04550, 0.04656, 0.04261, 0.06069, 0.06218, 0.05664, 0.07707, 0.07311}, kHT1 = {0.000946, 0.001685, 0.003029, 0.003447, 0.002216, 0.004044, 0.004603, 0.003920, 0.006277}, kHT2 = {0.001081, 0.003060, 0.003131, 0.002395, 0.004081, 0.004182, 0.003184, 0.005183, 0.004305}, lam = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, mDot = 162.3880, mFN1 = {13.8215, 46.8650, 25.6293, 32.0366, 79.2677, 43.2952, 54.1876, 217.8490, 147.8261}, mFN2 = {14.4622, 31.9451, 25.6293, 62.4256, 54.1876, 43.2952, 105.4462, 126.6819, 248.0549}, mGN = {71.0660, 214.6193, 163.0457, 208.7310, 363.0457, 275.9391, 353.0964, 975.8376, 956.4467}, mMix = 324.7760, omega = 1, regionCoastDown = 0.02, regionFreeConvFF = 0.01, regionTripTime = {200000, 200000, 200000, 200000}, rho_0 = 0.002465140767843, sinReact = 0, sinTime = 20000000, stepReact = 100E-5, stepTime = 2000, tauCore = 8.460000000000001, tauCoreDHRS = 3.7700, tauLoop = 16.73) annotation(
      Placement(visible = true, transformation(origin = {-727, -27}, extent = {{-399, -399}, {399, 399}}, rotation = 0)));
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantHeatRemoval constantHeatRemoval(cP = 1.9665E-3, mDot = 1.623879934566580e+02, nodeMass = 800, qDotOut = 8, tempInSetPoint = 663.5300, tempoutSetPoint = 637) annotation(
      Placement(visible = true, transformation(origin = {541, 147}, extent = {{375, -375}, {-375, 375}}, rotation = 0)));
  equation
    connect(primaryPump.primaryFlowFrac, nineRcore.flowFraction_In) annotation(
      Line(points = {{-626, -1506}, {-1464, -1506}, {-1464, -330}, {-998, -330}}, color = {245, 121, 0}));
    connect(constantHeatRemoval.temp_Out, nineRcore.temp_In) annotation(
      Line(points = {{392, 148}, {0, 148}, {0, 90}, {-456, 90}, {-456, -322}}));
    connect(nineRcore.temp_Out, constantHeatRemoval.temp_In) annotation(
      Line(points = {{-468, 264}, {38, 264}, {38, 440}, {1030, 440}, {1030, 148}, {692, 148}}));
    annotation(
      Diagram(coordinateSystem(extent = {{-1720, 1020}, {2940, -1680}})));
  end nineEcoreTest;
  
  model nineRcoreNoDH
    parameter SMD_MSR_Modelica.Units.ReactorPower P;
    parameter SMD_MSR_Modelica.Units.Mass TotalFuelMass;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_fuel;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_grap;
    parameter SMD_MSR_Modelica.Units.Convection hA[9];
    parameter SMD_MSR_Modelica.Units.VolumeImportance kFN1[9];
    parameter SMD_MSR_Modelica.Units.VolumeImportance kFN2[9];
    parameter SMD_MSR_Modelica.Units.HeatTransferFraction kHT1[9];
    parameter SMD_MSR_Modelica.Units.HeatTransferFraction kHT2[9];
    parameter SMD_MSR_Modelica.Units.Mass mMix;
    parameter SMD_MSR_Modelica.Units.Mass mFN1[9];
    parameter SMD_MSR_Modelica.Units.Mass mFN2[9];
    parameter SMD_MSR_Modelica.Units.Mass mGN[9];
    parameter SMD_MSR_Modelica.Units.MassFlowRate mDot;
    parameter SMD_MSR_Modelica.Units.Temperature TF1_0[9];
    parameter SMD_MSR_Modelica.Units.Temperature TF2_0[9];
    parameter SMD_MSR_Modelica.Units.Temperature TG_0[9];
    parameter SMD_MSR_Modelica.Units.Temperature Tmix_0;
    parameter Real IF1[9];
    parameter Real IF2[9];
    parameter Real IG[9];
    parameter SMD_MSR_Modelica.Units.TemperatureReactivityCoef aF;
    parameter SMD_MSR_Modelica.Units.TemperatureReactivityCoef aG;
    parameter SMD_MSR_Modelica.Units.Reactivity rho_0;
    parameter SMD_MSR_Modelica.Units.Reactivity stepReact;
    parameter SMD_MSR_Modelica.Units.InitiationTime stepTime;
    parameter SMD_MSR_Modelica.Units.Reactivity sinReact;
    parameter SMD_MSR_Modelica.Units.Frequency omega;
    parameter SMD_MSR_Modelica.Units.InitiationTime sinTime;
    parameter SMD_MSR_Modelica.Units.NeutronGenerationTime Lam;
    parameter SMD_MSR_Modelica.Units.PrecursorConc beta[6];
    parameter SMD_MSR_Modelica.Units.PrecursorDecayConstant lam[6];
    parameter SMD_MSR_Modelica.Units.ResidentTime tauCore;
    parameter SMD_MSR_Modelica.Units.ResidentTime tauLoop;
    parameter SMD_MSR_Modelica.Units.ResidentTime tauCoreDHRS;
    parameter SMD_MSR_Modelica.Units.FlowFraction flowFracRegions[4];
    parameter SMD_MSR_Modelica.Units.InitiationTime regionTripTime[4];
    parameter SMD_MSR_Modelica.Units.PumpConstant regionCoastDown;
    parameter SMD_MSR_Modelica.Units.FlowFraction regionFreeConvFF;
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = P, TotalFuelMass = TotalFuelMass) annotation(
      Placement(visible = true, transformation(origin = {-298, 228}, extent = {{-82, -82}, {82, 82}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region1(P = P, TF1_0 = TF1_0[1], TF2_0 = TF2_0[1], TG_0 = TG_0[1], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[1], kFN1 = kFN1[1], kFN2 = kFN2[1], kG = kHT1[1] + kHT2[1], kHT_FN1 = kHT1[1], kHT_FN2 = kHT2[1], m_FN1 = mFN1[1], m_FN2 = mFN2[1], m_GN = mGN[1], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[1], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {-288, 60}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region2(P = P, TF1_0 = TF1_0[2], TF2_0 = TF2_0[2], TG_0 = TG_0[2], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[2], kFN1 = kFN1[2], kFN2 = kFN2[2], kG = kHT1[2] + kHT2[2], kHT_FN1 = kHT1[2], kHT_FN2 = kHT2[2], m_FN1 = mFN1[2], m_FN2 = mFN2[2], m_GN = mGN[2], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[2], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {-26, -144}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region3(P = P, TF1_0 = TF1_0[3], TF2_0 = TF2_0[3], TG_0 = TG_0[3], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[3], kFN1 = kFN1[3], kFN2 = kFN2[3], kG = kHT1[3] + kHT2[3], kHT_FN1 = kHT1[3], kHT_FN2 = kHT2[3], m_FN1 = mFN1[3], m_FN2 = mFN2[3], m_GN = mGN[3], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[2], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {-28, 58}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region4(P = P, TF1_0 = TF1_0[4], TF2_0 = TF2_0[4], TG_0 = TG_0[4], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[4], kFN1 = kFN1[4], kFN2 = kFN2[4], kG = kHT1[4] + kHT2[4], kHT_FN1 = kHT1[4], kHT_FN2 = kHT2[4], m_FN1 = mFN1[4], m_FN2 = mFN2[4], m_GN = mGN[4], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[2], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {-27, 235}, extent = {{-57, -57}, {57, 57}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region5(P = P, TF1_0 = TF1_0[5], TF2_0 = TF2_0[5], TG_0 = TG_0[5], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[5], kFN1 = kFN1[5], kFN2 = kFN2[5], kG = kHT1[5] + kHT2[5], kHT_FN1 = kHT1[5], kHT_FN2 = kHT2[5], m_FN1 = mFN1[5], m_FN2 = mFN2[5], m_GN = mGN[5], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[3], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {208, -146}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region6(P = P, TF1_0 = TF1_0[6], TF2_0 = TF2_0[6], TG_0 = TG_0[6], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[6], kFN1 = kFN1[6], kFN2 = kFN2[6], kG = kHT1[6] + kHT2[6], kHT_FN1 = kHT1[6], kHT_FN2 = kHT2[6], m_FN1 = mFN1[6], m_FN2 = mFN2[6], m_GN = mGN[6], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[3], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {212, 54}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region7(P = P, TF1_0 = TF1_0[7], TF2_0 = TF2_0[7], TG_0 = TG_0[7], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[7], kFN1 = kFN1[7], kFN2 = kFN2[7], kG = kHT1[7] + kHT2[7], kHT_FN1 = kHT1[7], kHT_FN2 = kHT2[7], m_FN1 = mFN1[7], m_FN2 = mFN2[7], m_GN = mGN[7], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[3], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {212, 234}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region8(P = P, TF1_0 = TF1_0[8], TF2_0 = TF2_0[8], TG_0 = TG_0[8], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[8], kFN1 = kFN1[8], kFN2 = kFN2[8], kG = kHT1[8] + kHT2[8], kHT_FN1 = kHT1[8], kHT_FN2 = kHT2[8], m_FN1 = mFN1[8], m_FN2 = mFN2[8], m_GN = mGN[8], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[4], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {472, -18}, extent = {{-58, -58}, {58, 58}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.FuelChannel Region9(P = P, TF1_0 = TF1_0[9], TF2_0 = TF2_0[9], TG_0 = TG_0[9], cP_fuel = cP_fuel, cP_graphite = cP_grap, hAnom = hA[9], kFN1 = kFN1[9], kFN2 = kFN2[9], kG = kHT1[9] + kHT2[9], kHT_FN1 = kHT1[9], kHT_FN2 = kHT2[9], m_FN1 = mFN1[9], m_FN2 = mFN2[9], m_GN = mGN[9], mdot_fuelNom = mDot, regionFlowFrac = flowFracRegions[4], tauCoreToDHRS = 0) annotation(
      Placement(visible = true, transformation(origin = {469, 159}, extent = {{-57, -57}, {57, 57}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region1FB(FuelTempSetPointNode1 = TF1_0[1], FuelTempSetPointNode2 = TF2_0[1], GrapTempSetPoint = TG_0[1], IF1 = IF1[1], IF2 = IF2[1], IG = IG[1], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {-188, 46}, extent = {{-56, 56}, {56, -56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region2FB(FuelTempSetPointNode1 = TF1_0[2], FuelTempSetPointNode2 = TF2_0[2], GrapTempSetPoint = TG_0[2], IF1 = IF1[2], IF2 = IF2[2], IG = IG[2], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {73, -155}, extent = {{-55, 55}, {55, -55}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region3FB(FuelTempSetPointNode1 = TF1_0[3], FuelTempSetPointNode2 = TF2_0[3], GrapTempSetPoint = TG_0[3], IF1 = IF1[3], IF2 = IF2[3], IG = IG[3], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {74, 40}, extent = {{-56, 56}, {56, -56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region4FB(FuelTempSetPointNode1 = TF1_0[4], FuelTempSetPointNode2 = TF2_0[4], GrapTempSetPoint = TG_0[4], IF1 = IF1[4], IF2 = IF2[4], IG = IG[4], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {76, 220}, extent = {{-54, 54}, {54, -54}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region5FB(FuelTempSetPointNode1 = TF1_0[5], FuelTempSetPointNode2 = TF2_0[5], GrapTempSetPoint = TG_0[5], IF1 = IF1[5], IF2 = IF2[5], IG = IG[5], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {302, -154}, extent = {{-56, 56}, {56, -56}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region6FB(FuelTempSetPointNode1 = TF1_0[6], FuelTempSetPointNode2 = TF2_0[6], GrapTempSetPoint = TG_0[6], IF1 = IF1[6], IF2 = IF2[6], IG = IG[6], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {294, 34}, extent = {{-54, 54}, {54, -54}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region7FB(FuelTempSetPointNode1 = TF1_0[7], FuelTempSetPointNode2 = TF2_0[7], GrapTempSetPoint = TG_0[7], IF1 = IF1[7], IF2 = IF2[7], IG = IG[7], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {304, 222}, extent = {{-54, 54}, {54, -54}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region8FB(FuelTempSetPointNode1 = TF1_0[8], FuelTempSetPointNode2 = TF2_0[8], GrapTempSetPoint = TG_0[8], IF1 = IF1[8], IF2 = IF2[8], IG = IG[8], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {606, -34}, extent = {{-58, 58}, {58, -58}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback Region9FB(FuelTempSetPointNode1 = TF1_0[9], FuelTempSetPointNode2 = TF2_0[9], GrapTempSetPoint = TG_0[9], IF1 = IF1[9], IF2 = IF2[9], IG = IG[9], a_F = aF, a_G = aG) annotation(
      Placement(visible = true, transformation(origin = {586, 146}, extent = {{-56, 56}, {56, -56}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.FlowDistributor flowDistributor(coastDownK = 1, freeConvectionFF = regionFreeConvFF, numOutput = 4, regionTripTime = regionTripTime) annotation(
      Placement(visible = true, transformation(origin = {11, -455}, extent = {{-53, -53}, {53, 53}}, rotation = 180)));
    SMD_MSR_Modelica.HeatTransport.MixingPot mixingPot(Cp = cP_fuel, T_0 = Tmix_0, m = mMix, mDotNom = mDot, numInput = 4, regionFlowFrac = flowFracRegions, tauPotToDHRS = tauCoreDHRS) annotation(
      Placement(visible = true, transformation(origin = {80, 442}, extent = {{-110, -110}, {110, 110}}, rotation = 0)));
    input SMD_MSR_Modelica.PortsConnectors.Temp_In temp_In annotation(
      Placement(visible = true, transformation(origin = {180, -400}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {68, -74}, extent = {{-28, -28}, {28, 28}}, rotation = 0)));
    input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In flowFraction_In annotation(
      Placement(visible = true, transformation(origin = {80, -434}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-68, -76}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
    output SMD_MSR_Modelica.PortsConnectors.Temp_Out temp_Out annotation(
      Placement(visible = true, transformation(origin = {68, 510}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {65, 73}, extent = {{-25, -25}, {25, 25}}, rotation = 0)));
    output SMD_MSR_Modelica.PortsConnectors.DecayHeat_Out decayHeat_Out annotation(
      Placement(visible = true, transformation(origin = {10, 510}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-65, 73}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.SumReactivity sumReactivity(numInput = 9) annotation(
      Placement(visible = true, transformation(origin = {-259, -379}, extent = {{-105, -105}, {105, 105}}, rotation = 0)));
    SMD_MSR_Modelica.Nuclear.mPKE mPKE(Lam = Lam, beta = beta, lam = lam, omega = omega, rho_0 = rho_0, sinInsertionTime = sinTime, sin_mag = sinReact, stepInsertionTime = stepTime, step_mag = stepReact, tauCoreNom = tauCore, tauLoopNom = tauLoop) annotation(
      Placement(visible = true, transformation(origin = {-432, -274}, extent = {{-116, -116}, {116, 116}}, rotation = 180)));
  SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
      Placement(visible = true, transformation(origin = {-449, -97}, extent = {{-87, -87}, {87, 87}}, rotation = 0)));
  initial equation
  //temp_In.T = 6.3222E+02;
  equation
    connect(Region2.temp_Out, Region3.temp_In) annotation(
      Line(points = {{2, -114}, {12, -114}, {12, -18}, {-73, -18}, {-73, 19}}));
    connect(Region3.temp_Out, Region4.temp_In) annotation(
      Line(points = {{0, 88}, {12, 88}, {12, 106}, {-73, 106}, {-73, 195}}));
    connect(Region5.temp_Out, Region6.temp_In) annotation(
      Line(points = {{236, -116}, {236, -26}, {167, -26}, {167, 15}}));
    connect(Region6.temp_Out, Region7.temp_In) annotation(
      Line(points = {{240, 84}, {240, 106}, {167, 106}, {167, 195}}));
    connect(Region8.temp_Out, Region9.temp_In) annotation(
      Line(points = {{501, 13}, {486, 13}, {486, 54}, {423, 54}, {423, 119}}));
    connect(Region1.fuelNode1, Region1FB.fuelNode1) annotation(
      Line(points = {{-260, 20}, {-210, 20}, {-210, 27}}));
    connect(Region1.fuelNode2, Region1FB.fuelNode2) annotation(
      Line(points = {{-260, 43}, {-234, 43}, {-234, 46}, {-210, 46}}));
    connect(Region1.grapNode, Region1FB.grapNode) annotation(
      Line(points = {{-260, 69}, {-236, 69}, {-236, 63}, {-210, 63}}));
    connect(Region2.fuelNode1, Region2FB.fuelNode1) annotation(
      Line(points = {{2, -184}, {30.5, -184}, {30.5, -174}, {51, -174}}));
    connect(Region2.fuelNode2, Region2FB.fuelNode2) annotation(
      Line(points = {{2, -161}, {29.5, -161}, {29.5, -155}, {51, -155}}));
    connect(Region2.grapNode, Region2FB.grapNode) annotation(
      Line(points = {{2, -135}, {20.5, -135}, {20.5, -138.5}, {51, -138.5}}));
    connect(Region3.fuelNode1, Region3FB.fuelNode1) annotation(
      Line(points = {{0, 18}, {52, 18}, {52, 21}}));
    connect(Region3.fuelNode2, Region3FB.fuelNode2) annotation(
      Line(points = {{0, 41}, {26, 41}, {26, 40}, {52, 40}}));
    connect(Region3.grapNode, Region3FB.grapNode) annotation(
      Line(points = {{0, 67}, {52, 67}, {52, 57}}));
    connect(Region4.fuelNode1, Region4FB.fuelNode1) annotation(
      Line(points = {{2, 194}, {29, 194}, {29, 202}, {54, 202}}));
    connect(Region4.fuelNode2, Region4FB.fuelNode2) annotation(
      Line(points = {{2, 218}, {33, 218}, {33, 220}, {54, 220}}));
    connect(Region4.grapNode, Region4FB.grapNode) annotation(
      Line(points = {{2, 244}, {28, 244}, {28, 236}, {54, 236}}));
    connect(Region5.fuelNode1, Region5FB.fuelNode1) annotation(
      Line(points = {{236, -186}, {280, -186}, {280, -173}}));
    connect(Region5.fuelNode2, Region5FB.fuelNode2) annotation(
      Line(points = {{236, -163}, {280, -163}, {280, -154}}));
    connect(Region5.grapNode, Region5FB.grapNode) annotation(
      Line(points = {{236, -137}, {280, -137}}));
    connect(Region6.fuelNode1, Region6FB.fuelNode1) annotation(
      Line(points = {{240, 14}, {272, 14}, {272, 16}}));
    connect(Region6.fuelNode2, Region6FB.fuelNode2) annotation(
      Line(points = {{240, 38}, {256, 38}, {256, 34}, {272, 34}}));
    connect(Region6.grapNode, Region6FB.grapNode) annotation(
      Line(points = {{240, 62}, {255, 62}, {255, 50}, {272, 50}}));
    connect(Region7.fuelNode1, Region7FB.fuelNode1) annotation(
      Line(points = {{240, 194}, {282, 194}, {282, 204}}));
    connect(Region7.fuelNode2, Region7FB.fuelNode2) annotation(
      Line(points = {{240, 217}, {282, 217}, {282, 222}}));
    connect(Region7.grapNode, Region7FB.grapNode) annotation(
      Line(points = {{240, 243}, {282, 243}, {282, 238}}));
    connect(Region8.fuelNode1, Region8FB.fuelNode1) annotation(
      Line(points = {{501, -60}, {542, -60}, {542, -54}, {583, -54}}));
    connect(Region8.fuelNode2, Region8FB.fuelNode2) annotation(
      Line(points = {{501, -35}, {508.5, -35}, {508.5, -34}, {583, -34}}));
    connect(Region8.grapNode, Region8FB.grapNode) annotation(
      Line(points = {{501, -9}, {539, -9}, {539, -17}, {583, -17}}));
    connect(Region9.fuelNode1, Region9FB.fuelNode1) annotation(
      Line(points = {{497.5, 118}, {507.75, 118}, {507.75, 127}, {564, 127}}));
    connect(Region9.fuelNode2, Region9FB.fuelNode2) annotation(
      Line(points = {{497.5, 142}, {531.75, 142}, {531.75, 146}, {564, 146}}));
    connect(Region9.grapNode, Region9FB.grapNode) annotation(
      Line(points = {{497.5, 168}, {515.75, 168}, {515.75, 163}, {564, 163}}));
    connect(temp_In, Region1.temp_In) annotation(
      Line(points = {{180, -400}, {-333, -400}, {-333, 21}}, color = {28, 113, 216}));
    connect(temp_In, Region2.temp_In) annotation(
      Line(points = {{180, -400}, {-71, -400}, {-71, -183}}, color = {28, 113, 216}));
    connect(temp_In, Region5.temp_In) annotation(
      Line(points = {{180, -400}, {180, -269}, {163, -269}, {163, -185}}, color = {28, 113, 216}));
    connect(temp_In, Region8.temp_In) annotation(
      Line(points = {{180, -400}, {426, -400}, {426, -59}}, color = {28, 113, 216}));
    connect(flowFraction_In, flowDistributor.flowFracIn) annotation(
      Line(points = {{80, -434}, {33, -434}}));
    connect(flowDistributor.flowFracOut[1], Region1.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-310, -434}, {-310, 15}}));
    connect(flowDistributor.flowFracOut[2], Region2.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -320.5}, {-48, -320.5}, {-48, -189}}));
    connect(flowDistributor.flowFracOut[2], Region3.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -220}, {-50, -220}, {-50, 13}}));
    connect(flowDistributor.flowFracOut[2], Region4.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -131}, {-50, -131}, {-50, 190}}));
    connect(flowDistributor.flowFracOut[3], Region5.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -264}, {186, -264}, {186, -191}}));
    connect(flowDistributor.flowFracOut[3], Region7.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -264}, {190, -264}, {190, 190}}));
    connect(flowDistributor.flowFracOut[3], Region6.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -266}, {190, -266}, {190, 10}}));
    connect(flowDistributor.flowFracOut[4], Region8.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -266}, {449, -266}, {449, -64}}));
    connect(flowDistributor.flowFracOut[4], Region9.fuelFlowFraction) annotation(
      Line(points = {{-10, -434}, {-10, -268}, {446, -268}, {446, 113}}));
    connect(Region1.temp_Out, mixingPot.temp_In[1]) annotation(
      Line(points = {{-260, 90}, {-192, 90}, {-192, 440}, {14, 440}}));
    connect(Region4.temp_Out, mixingPot.temp_In[2]) annotation(
      Line(points = {{2, 266}, {2, 357}, {14, 357}, {14, 440}}));
    connect(Region7.temp_Out, mixingPot.temp_In[3]) annotation(
      Line(points = {{240, 264}, {286, 264}, {286, 314}, {14, 314}, {14, 440}}));
    connect(Region9.temp_Out, mixingPot.temp_In[4]) annotation(
      Line(points = {{497.5, 189}, {472, 189}, {472, 316}, {14, 316}, {14, 440}}));
    connect(mixingPot.temp_Out, temp_Out) annotation(
      Line(points = {{146, 442}, {142, 442}, {142, 510}, {68, 510}}));
    connect(flowDistributor.flowFracOut[1], mixingPot.flowFraction[1]) annotation(
      Line(points = {{-10, -434}, {-538, -434}, {-538, 376}, {16, 376}}, color = {245, 121, 0}, thickness = 0.5));
    connect(flowDistributor.flowFracOut[2], mixingPot.flowFraction[2]) annotation(
      Line(points = {{-10, -434}, {-564, -434}, {-564, 376}, {16, 376}}, color = {245, 121, 0}, thickness = 0.5));
    connect(flowDistributor.flowFracOut[3], mixingPot.flowFraction[3]) annotation(
      Line(points = {{-10, -434}, {-552, -434}, {-552, 376}, {16, 376}}, color = {245, 121, 0}, thickness = 0.5));
    connect(flowDistributor.flowFracOut[4], mixingPot.flowFraction[4]) annotation(
      Line(points = {{-10, -434}, {-528, -434}, {-528, 376}, {16, 376}}, color = {245, 121, 0}, thickness = 0.5));
    connect(Region1FB.feedback, sumReactivity.reactivityIn[1]) annotation(
      Line(points = {{-168, 36}, {-127, 36}, {-127, -322}, {-126, -322}, {-126, -358}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region2FB.feedback, sumReactivity.reactivityIn[2]) annotation(
      Line(points = {{92, -165}, {-238, -165}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region3FB.feedback, sumReactivity.reactivityIn[3]) annotation(
      Line(points = {{93, 30}, {39, 30}, {39, -358}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region4FB.feedback, sumReactivity.reactivityIn[4]) annotation(
      Line(points = {{94, 210}, {-238, 210}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region5FB.feedback, sumReactivity.reactivityIn[5]) annotation(
      Line(points = {{322, -164}, {316, -164}, {316, -358}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region6FB.feedback, sumReactivity.reactivityIn[6]) annotation(
      Line(points = {{312, 24}, {352, 24}, {352, -358}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region7FB.feedback, sumReactivity.reactivityIn[7]) annotation(
      Line(points = {{322, 212}, {-238, 212}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region8FB.feedback, sumReactivity.reactivityIn[8]) annotation(
      Line(points = {{626, -44}, {556, -44}, {556, -358}, {-238, -358}}, color = {78, 154, 6}));
    connect(Region9FB.feedback, sumReactivity.reactivityIn[9]) annotation(
      Line(points = {{605, 136}, {664.5, 136}, {664.5, 44}, {663, 44}, {663, -358}, {-238, -358}}, color = {78, 154, 6}));
    connect(sumReactivity.reactivityOut, mPKE.feedback) annotation(
      Line(points = {{-280, -358}, {-432, -358}, {-432, -320}}, color = {78, 154, 6}));
    connect(flowFraction_In, mPKE.fuelFlowFrac) annotation(
      Line(points = {{80, -434}, {80, -480}, {-471, -480}, {-471, -320}}, color = {255, 120, 0}));
    connect(mPKE.n_population, powerBlock.nPop) annotation(
      Line(points = {{-432, -230}, {-398, -230}, {-398, 196}, {-340, 196}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region1.nPop) annotation(
      Line(points = {{-432, -230}, {-368, -230}, {-368, 92}, {-333, 92}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region2.nPop) annotation(
      Line(points = {{-432, -230}, {-156, -230}, {-156, -112}, {-70, -112}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region3.nPop) annotation(
      Line(points = {{-432, -230}, {-116, -230}, {-116, 90}, {-73, 90}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region4.nPop) annotation(
      Line(points = {{-432, -230}, {-122, -230}, {-122, 268}, {-72, 268}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region5.nPop) annotation(
      Line(points = {{-432, -230}, {130, -230}, {130, -114}, {163, -114}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region6.nPop) annotation(
      Line(points = {{-432, -230}, {136, -230}, {136, 86}, {168, 86}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region7.nPop) annotation(
      Line(points = {{-432, -230}, {138, -230}, {138, 266}, {168, 266}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region8.nPop) annotation(
      Line(points = {{-432, -230}, {396, -230}, {396, 16}, {426, 16}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Region9.nPop) annotation(
      Line(points = {{-432, -230}, {404, -230}, {404, 192}, {424, 192}}, color = {20, 36, 248}));
  connect(constantDecayPower.decayHeat_Out, powerBlock.P_decay) annotation(
      Line(points = {{-448, -62}, {-448, 67}, {-450, 67}, {-450, 214}, {-359, 214}, {-359, 196}, {-272, 196}}, color = {0, 225, 255}));
  connect(constantDecayPower.decayHeat_Out, mixingPot.decayHeat) annotation(
      Line(points = {{-448, -62}, {-446, -62}, {-446, 374}, {144, 374}}, color = {0, 225, 255}));
  connect(constantDecayPower.decayHeat_Out, Region1.P_decay) annotation(
      Line(points = {{-448, -62}, {-428, -62}, {-428, 56}, {-332, 56}}, color = {0, 225, 255}));
  connect(constantDecayPower.decayHeat_Out, Region2.P_decay) annotation(
      Line(points = {{-448, -62}, {-190, -62}, {-190, -148}, {-70, -148}}, color = {0, 225, 255}));
  connect(constantDecayPower.decayHeat_Out, Region3.P_decay) annotation(
      Line(points = {{-448, -62}, {-104, -62}, {-104, 54}, {-72, 54}}, color = {0, 225, 255}));
  connect(constantDecayPower.decayHeat_Out, Region4.P_decay) annotation(
      Line(points = {{-448, -62}, {-104, -62}, {-104, 232}, {-72, 232}}, color = {0, 225, 255}));
  connect(constantDecayPower.decayHeat_Out, Region5.P_decay) annotation(
      Line(points = {{-448, -62}, {116, -62}, {116, -150}, {164, -150}}, color = {0, 225, 255}));
  connect(constantDecayPower.decayHeat_Out, Region6.P_decay) annotation(
      Line(points = {{-448, -62}, {118, -62}, {118, 50}, {168, 50}}, color = {0, 225, 255}));
  connect(constantDecayPower.decayHeat_Out, Region7.P_decay) annotation(
      Line(points = {{-448, -62}, {126, -62}, {126, 230}, {168, 230}}, color = {0, 225, 255}));
  connect(constantDecayPower.decayHeat_Out, Region8.P_decay) annotation(
      Line(points = {{-448, -62}, {372, -62}, {372, -22}, {426, -22}}, color = {0, 225, 255}));
  connect(constantDecayPower.decayHeat_Out, Region9.P_decay) annotation(
      Line(points = {{-448, -62}, {378, -62}, {378, 156}, {424, 156}}, color = {0, 225, 255}));
  connect(constantDecayPower.decayHeat_Out, decayHeat_Out) annotation(
      Line(points = {{-448, -62}, {-444, -62}, {-444, 510}, {10, 510}}, color = {0, 225, 255}));
    annotation(
      Diagram(coordinateSystem(extent = {{-560, 580}, {600, -520}}), graphics = {Text(origin = {-360, 106}, extent = {{-26, 10}, {26, -10}}, textString = "Region1"), Text(origin = {-110, -96}, extent = {{-26, 10}, {26, -10}}, textString = "Region2"), Text(origin = {-112, 104}, extent = {{-26, 10}, {26, -10}}, textString = "Region3"), Text(origin = {-110, 286}, extent = {{-26, 10}, {26, -10}}, textString = "Region4"), Text(origin = {132, -96}, extent = {{-26, 10}, {26, -10}}, textString = "Region5"), Text(origin = {132, 104}, extent = {{-26, 10}, {26, -10}}, textString = "Region6"), Text(origin = {136, 286}, extent = {{-26, 10}, {26, -10}}, textString = "Region7"), Text(origin = {366, 32}, extent = {{-26, 10}, {26, -10}}, textString = "Region8"), Text(origin = {366, 196}, extent = {{-26, 10}, {26, -10}}, textString = "Region9")}),
      Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Text(origin = {-1, -2}, extent = {{-71, 28}, {71, -28}}, textString = "9Rcore")}));
  end nineRcoreNoDH;
  
  model nineE_MSREnoDH
    parameter Real CpCoefFuel = 1;
    parameter Real HTCcoefCore = 1;
    parameter Real HTCcoefPHX = 1;
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.05, primaryPumpK = 0.2, tripPrimaryPump = 2000000) annotation(
      Placement(visible = true, transformation(origin = {-626, -1274}, extent = {{-386, -386}, {386, 386}}, rotation = 0)));
    MSREnineR.nineRcoreNoDH nineRcoreNoDH(IF1 = {0.02168, 0.02197, 0.07897, 0.08249, 0.02254, 0.08255, 0.08623, 0.02745, 0.06936}, IF2 = {0.02678, 0.06519, 0.08438, 0.04124, 0.06801, 0.08823, 0.04290, 0.05529, 0.03473}, IG = {0.04443, 0.08835, 0.16671, 0.12077, 0.09181, 0.17429, 0.12612, 0.08408, 0.10343}, Lam = 2.400E-04, P = 8, TF1_0 = {641.0644, 636.1926, 662.2593, 692.4323, 634.6308, 650.7339, 669.3933, 633.4354, 638.9545}, TF2_0 = {651.9893, 647.9801, 677.5167, 699.8892, 641.9048, 660.1706, 673.9817, 635.8835, 640.4922}, TG_0 = {664.0462, 653.7083, 692.6790, 714.5485, 648.7869, 674.7077, 686.8285, 640.8657, 647.7290}, Tmix_0 = 663.5300, TotalFuelMass = 4.093499709644400e+03, aF = -8.71E-05, aG = -6.66E-05, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, cP_fuel = 0.0020*CpCoefFuel, cP_grap = 0.0018, flowFracRegions = {0.061409895150439, 0.138549617100260, 0.234231239311780, 0.565809248437521}, hA = {0.0007056, 0.0021672, 0.00162, 0.0021132, 0.0035586, 0.002745, 0.003573, 0.009801, 0.009648}, kFN1 = {0.01493, 0.02736, 0.04504, 0.05126, 0.03601, 0.06014, 0.06845, 0.06179, 0.09333}, kFN2 = {0.01721, 0.04550, 0.04656, 0.04261, 0.06069, 0.06218, 0.05664, 0.07707, 0.07311}, kHT1 = {0.000946, 0.001685, 0.003029, 0.003447, 0.002216, 0.004044, 0.004603, 0.003920, 0.006277}, kHT2 = {0.001081, 0.003060, 0.003131, 0.002395, 0.004081, 0.004182, 0.003184, 0.005183, 0.004305}, lam = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, mDot = 162.3880, mFN1 = {13.8215, 46.8650, 25.6293, 32.0366, 79.2677, 43.2952, 54.1876, 217.8490, 147.8261}, mFN2 = {14.4622, 31.9451, 25.6293, 62.4256, 54.1876, 43.2952, 105.4462, 126.6819, 248.0549}, mGN = {71.0660, 214.6193, 163.0457, 208.7310, 363.0457, 275.9391, 353.0964, 975.8376, 956.4467}, mMix = 324.7760, omega = 1, regionCoastDown = 0.02, regionFreeConvFF = 0.01, regionTripTime = {200000, 200000, 200000, 200000}, rho_0 = 0.002465140767843, sinReact = 0, sinTime = 20000000, stepReact = 100E-5, stepTime = 2000, tauCore = 8.46, tauCoreDHRS = 1.3850, tauLoop = 16.73) annotation(
      Placement(visible = true, transformation(origin = {-1349, -215}, extent = {{-399, -399}, {399, 399}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.05, secondaryPumpK = 0.2, tripSecondaryPump = 200000) annotation(
      Placement(visible = true, transformation(origin = {1204, -1400}, extent = {{-362, -362}, {362, 362}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.UHX uhx(SetDemand = 1, UHXP = 8, cP_pFluidUHX = 2.39E-3, m_PN1UHX = 3.924401289158446e+02, m_dot_pFluidUHXnom = 100.5793, tauUHXtoPHX = 8.24, tripUHX = 2000000) annotation(
      Placement(visible = true, transformation(origin = {2441, -231}, extent = {{-567, 567}, {567, -567}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.PrimaryHeatExchanger PHX(T_PN1_initial = 6.478750000000000e+02, cP_Tube = 5.778E-04, cP_pFluid = 0.0020*CpCoefFuel, cP_sFluid = 2.39E-3, hApnNom = 0.1620*HTCcoefPHX, hAsnNom = 0.0765, m_PN1 = 87.0959, m_PN2 = 87.0959, m_PN3 = 87.0959, m_PN4 = 87.0959, m_SN1 = 24.9200, m_SN2 = 24.9200, m_SN3 = 24.9200, m_SN4 = 24.9200, m_TN1 = 101.1662, m_TN2 = 101.1662, m_dot_pFluidNom = 162.3880, m_dot_sFluidNom = 100.5793, tauPHXtoPipe = 3.8350, tauPHXtoUHX = 4.71) annotation(
      Placement(visible = true, transformation(origin = {1691, 681}, extent = {{-427, -427}, {427, 427}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Cp_fluid = 0.0020*CpCoefFuel, DHRS_MaxPowerRm = 8, DHRS_PowerBleed = 0, DHRS_time = 20000, DHRS_timeConstant = 0.02, m_DHRS = 162.3880, m_Pi_C_PHX = 162.3880*3.7700, m_dotNom = 162.3880, tauDHRStoPHX = 1.3850)  annotation(
      Placement(visible = true, transformation(origin = {-183, 825}, extent = {{509, -509}, {-509, 509}}, rotation = 180)));
  SMD_MSR_Modelica.HeatTransport.Pipe pipe(Cp_fluidPi = 0.0020*CpCoefFuel, m_Pi_PHX_C = 162.3880*8.6700, m_dotPiNom = 162.3880, m_pi = 162.3880, tauPiToCore = 3.8350)  annotation(
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
      Line(points = {{-234, -754}, {-1078, -754}, {-1078, -510}}));
    connect(nineRcoreNoDH.temp_Out, dhrs.DHRS_tempIN) annotation(
      Line(points = {{-1090, 76}, {-1072, 76}, {-1072, 714}, {-488, 714}}));
    connect(nineRcoreNoDH.decayHeat_Out, dhrs.DHRS_DecayHeat) annotation(
      Line(points = {{-1608, 76}, {-1600, 76}, {-1600, 520}, {-488, 520}}, color = {0, 225, 255}));
    connect(nineRcoreNoDH.decayHeat_Out, pipe.PiDecay_Heat) annotation(
      Line(points = {{-1608, 76}, {-1600, 76}, {-1600, 284}, {168, 284}, {168, -576}}, color = {0, 225, 255}));
    connect(nineRcoreNoDH.decayHeat_Out, PHX.P_decay) annotation(
      Line(points = {{-1608, 76}, {-1600, 76}, {-1600, 1106}, {1810, 1106}, {1810, 894}}, color = {0, 225, 255}));
    connect(primaryPump.primaryFlowFrac, nineRcoreNoDH.flowFraction_In) annotation(
      Line(points = {{-626, -1506}, {-1620, -1506}, {-1620, -518}}, color = {245, 121, 0}));
    annotation(
      Diagram(coordinateSystem(extent = {{-1480, 1160}, {3040, -1780}})));
  end nineE_MSREnoDH;
  annotation(
    uses(Modelica(version = "4.0.0")));
end MSREnineR;