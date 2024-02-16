package MSRE
  package Components
    model MSRE1R
      parameter Integer numGroups;
      parameter SMD_MSR_Modelica.Units.DecayConstant lambda[numGroups];
      parameter SMD_MSR_Modelica.Units.DelayedNeutronFrac beta[numGroups];
      parameter SMD_MSR_Modelica.Units.NeutronGenerationTime LAMBDA;
      parameter SMD_MSR_Modelica.Units.TemperatureReactivityCoef a_F;
      parameter SMD_MSR_Modelica.Units.TemperatureReactivityCoef a_G;
      parameter SMD_MSR_Modelica.Units.NominalNeutronPopulation n_0;
      parameter Boolean addRho0;
      parameter SMD_MSR_Modelica.Units.Density rho_fuel;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_fuel;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel;
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel;
      parameter SMD_MSR_Modelica.Units.Density rho_grap;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_grap;
      parameter SMD_MSR_Modelica.Units.Temperature TF1_0;
      parameter SMD_MSR_Modelica.Units.Temperature TF2_0;
      parameter SMD_MSR_Modelica.Units.Temperature TG_0;
      parameter Boolean EnableRad;
      parameter SMD_MSR_Modelica.Units.Temperature T_inf;
      SMD_MSR_Modelica.Nuclear.mPKE mpke(numGroups = numGroups, lambda = lambda, beta = beta, LAMBDA = LAMBDA, n_0 = n_0, addRho0 = addRho0, nomTauLoop = 16.73, nomTauCore = 8.46) annotation(
        Placement(transformation(origin = {-205.6, 84.4}, extent = {{-20.4, -20.4}, {13.6, 13.6}})));
      SMD_MSR_Modelica.Nuclear.FuelChannel fuelchannel(rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, TF1_0 = TF1_0, TF2_0 = TF2_0, TG_0 = TG_0, regionFlowFrac = 1, KF = kFuel, Ac = 1.5882, LF1 = 0.7875, LF2 = 0.7875, OutterRegion = EnableRad, ArF1 = 3.5442, ArF2 = 3.5442, e = 0, Tinf = T_inf, vol_FN1 = 0.3202, vol_FN2 = 0.3202, vol_GN = 1.95386, hAnom = 18000) annotation(
        Placement(transformation(origin = {-90.3333, -104.278}, extent = {{-67.2222, -67.2222}, {40.3333, 53.7778}})));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback react(FuelTempSetPointNode1 = TF1_0, FuelTempSetPointNode2 = TF2_0, GrapTempSetPoint = TG_0, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = a_F, a_G = a_G) annotation(
        Placement(transformation(origin = {3.69334, 57.8548}, extent = {{-22.2896, 14.8598}, {14.8598, -22.2896}})));
      SMD_MSR_Modelica.Nuclear.PowerBlock powerblock(P(displayUnit = "MW") = 8000000, TotalFuelVol = 1.9071) annotation(
        Placement(transformation(origin = {-201.2, 5.6}, extent = {{-12.8, -25.6}, {19.2, 6.4}})));
      SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(numGroups = 3, DHYG = {2.3751e-03, 8.7763e-05, 1.9596e-06}, DHlamG = {0.09453, 0.004420, 8.6098E-5}) annotation(
        Placement(transformation(origin = {-140, 18}, extent = {{-16, -16}, {24, 24}})));
      SMD_MSR_Modelica.Signals.Constants.ConstantNeutronSource constantNeutronSource(S = 0) annotation(
        Placement(transformation(origin = {-293, 79}, extent = {{-27, -27}, {27, 27}})));
      input SMD_MSR_Modelica.PortsConnectors.FlowFractionIn flowFracIn annotation(
        Placement(transformation(origin = {-370, -124}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-200, -80}, extent = {{-18, -18}, {18, 18}})));
      input SMD_MSR_Modelica.PortsConnectors.TempIn tempIn annotation(
        Placement(transformation(origin = {-368, -156}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-200, -158}, extent = {{-18, -18}, {18, 18}})));
      output SMD_MSR_Modelica.PortsConnectors.TempOut tempOut annotation(
        Placement(transformation(origin = {158, -76}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-120, -80}, extent = {{-18, -18}, {18, 18}})));
      output SMD_MSR_Modelica.PortsConnectors.VolumetircPowerOut volumetircPowerOut annotation(
        Placement(transformation(origin = {160, -14}, extent = {{-6, -6}, {6, 6}}), iconTransformation(origin = {-119, -159}, extent = {{-17, -17}, {17, 17}})));
      input SMD_MSR_Modelica.PortsConnectors.RealIn realIn annotation(
        Placement(transformation(origin = {-372, 148}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-160, -80}, extent = {{-18, -18}, {18, 18}})));
    equation
      connect(mpke.n_population, powerblock.nPop) annotation(
        Line(points = {{-209, 71}, {-208, 71}, {-208, 6}}, color = {20, 36, 248}, thickness = 1));
      connect(powerblock.decayPowerM, fuelchannel.decayHeat) annotation(
        Line(points = {{-188, -14}, {-188, -96}, {-139, -96}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.fissionPower, fuelchannel.fissionPower) annotation(
        Line(points = {{-208, -14}, {-208, -67}, {-139, -67}}, color = {129, 61, 156}, thickness = 1));
      connect(mpke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-209, 71}, {-170.75, 71}, {-170.75, 34}, {-136, 34}}, color = {20, 36, 248}, thickness = 1));
      connect(decayHeat.decayHeat_Out, powerblock.decayNP) annotation(
        Line(points = {{-136, 10}, {-136, 6}, {-188, 6}}, color = {0, 225, 255}, thickness = 1));
      connect(react.feedback, mpke.feedback) annotation(
        Line(points = {{-4, 62}, {-4, 108}, {-206, 108}, {-206, 92}}, color = {78, 154, 6}, thickness = 1));
      connect(flowFracIn, fuelchannel.fuelFlowFraction) annotation(
        Line(points = {{-370, -124}, {-137, -124}}, color = {255, 120, 0}, thickness = 1));
      connect(fuelchannel.fuelNode2, react.fuelNode2) annotation(
        Line(points = {{-64, -76}, {0, -76}, {0, 39}, {7, 39}}, color = {204, 0, 0}, thickness = 1));
      connect(fuelchannel.grapNode, react.grapNode) annotation(
        Line(points = {{-64, -110}, {12, -110}, {12, 39}, {-4, 39}}, color = {204, 0, 0}, thickness = 1));
      connect(fuelchannel.fuelNode1, react.fuelNode1) annotation(
        Line(points = {{-64, -144}, {-12, -144}, {-12, 39}, {-15, 39}}, color = {204, 0, 0}, thickness = 1));
      connect(flowFracIn, mpke.fuelFlowFrac) annotation(
        Line(points = {{-370, -124}, {-366, -124}, {-366, 122}, {-198, 122}, {-198, 92}}, color = {255, 120, 0}, thickness = 1));
      connect(tempIn, fuelchannel.temp_In) annotation(
        Line(points = {{-368, -156}, {-138, -156}, {-138, -152}}, color = {204, 0, 0}, thickness = 1));
      connect(fuelchannel.fuelNode2, tempOut) annotation(
        Line(points = {{-64, -76}, {158, -76}}, color = {204, 0, 0}, thickness = 1));
      connect(constantNeutronSource.neutronEmsRateOut, mpke.S) annotation(
        Line(points = {{-294, 88}, {-294, 110}, {-212, 110}, {-212, 92}}, color = {191, 64, 188}, thickness = 1));
      connect(powerblock.decayPowerM, volumetircPowerOut) annotation(
        Line(points = {{-188, -14}, {160, -14}}, color = {220, 138, 221}, thickness = 1));
      connect(realIn, mpke.ReactivityIn) annotation(
        Line(points = {{-372, 148}, {-220, 148}, {-220, 92}}, thickness = 1));
      annotation(
        Diagram(coordinateSystem(extent = {{-400, 180}, {180, -180}})),
        Icon(coordinateSystem(extent = {{-220, -60}, {-100, -180}}), graphics = {Rectangle(origin = {-160, -120}, lineThickness = 1.5, extent = {{-60, 60}, {60, -60}}), Text(origin = {-157, -119}, extent = {{-53, 31}, {53, -31}}, textString = "Reactor")}));
    end MSRE1R;

    model MSRE9R
      parameter Integer numGroups = 6;
      parameter SMD_MSR_Modelica.Units.DecayConstant lambda[numGroups];
      parameter SMD_MSR_Modelica.Units.DelayedNeutronFrac beta[numGroups];
      parameter SMD_MSR_Modelica.Units.NeutronGenerationTime LAMBDA;
      parameter SMD_MSR_Modelica.Units.NominalNeutronPopulation n_0;
      parameter SMD_MSR_Modelica.Units.TemperatureReactivityCoef aF;
      parameter SMD_MSR_Modelica.Units.TemperatureReactivityCoef aG;
      parameter Boolean addRho0;
      parameter SMD_MSR_Modelica.Units.Density rho_fuel;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_fuel;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel;
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel;
      parameter SMD_MSR_Modelica.Units.Density rho_grap;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_grap;
      parameter SMD_MSR_Modelica.Units.Volume volF1[9];
      parameter SMD_MSR_Modelica.Units.Volume volF2[9];
      parameter SMD_MSR_Modelica.Units.Volume volG[9];
      parameter SMD_MSR_Modelica.Units.Volume volUP;
      parameter SMD_MSR_Modelica.Units.Convection hA[9];
      parameter SMD_MSR_Modelica.Units.VolumeImportance kFN1[9];
      parameter SMD_MSR_Modelica.Units.VolumeImportance kFN2[9];
      parameter SMD_MSR_Modelica.Units.HeatTransferFraction kHT1[9];
      parameter SMD_MSR_Modelica.Units.HeatTransferFraction kHT2[9];
      parameter SMD_MSR_Modelica.Units.Area Ac[4];
      parameter SMD_MSR_Modelica.Units.Length LF1[9];
      parameter SMD_MSR_Modelica.Units.Length LF2[9];
      parameter SMD_MSR_Modelica.Units.Area ArF1[9];
      parameter SMD_MSR_Modelica.Units.Area ArF2[9];
      parameter SMD_MSR_Modelica.Units.Temperature TF1_0[9];
      parameter SMD_MSR_Modelica.Units.Temperature TF2_0[9];
      parameter SMD_MSR_Modelica.Units.Temperature TG_0[9];
      parameter SMD_MSR_Modelica.Units.Temperature Tmix_0;
      parameter Real IF1[9];
      parameter Real IF2[9];
      parameter Real IG[9];
      parameter SMD_MSR_Modelica.Units.FlowFraction flowFracRegions[4];
      parameter SMD_MSR_Modelica.Units.InitiationTime regionTripTime[4];
      parameter SMD_MSR_Modelica.Units.PumpConstant regionCoastDownK;
      parameter SMD_MSR_Modelica.Units.FlowFraction freeConvFF;
      parameter Boolean EnableRad;
      parameter SMD_MSR_Modelica.Units.Emissivity e;
      parameter SMD_MSR_Modelica.Units.Temperature T_inf;
      SMD_MSR_Modelica.Nuclear.mPKE mpke(numGroups = numGroups, lambda = lambda, beta = beta, LAMBDA = LAMBDA, n_0 = n_0, addRho0 = addRho0, nomTauLoop = 16.73, nomTauCore = 8.46) annotation(
        Placement(transformation(origin = {624.6, 2100.2}, extent = {{-321.6, -321.6}, {214.4, 214.4}}, rotation = -90)));
      SMD_MSR_Modelica.Nuclear.PowerBlock powerblock(P(displayUnit = "W") = 8E6, TotalFuelVol = 1.894347837) annotation(
        Placement(transformation(origin = {-828.2, 1326.8}, extent = {{162.4, -324.8}, {-243.6, 81.2}}, rotation = -90)));
      SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(numGroups = 3, DHYG = {2.3751e-03, 8.7763e-05, 1.9596e-06}, DHlamG = {0.09453, 0.004420, 8.6098E-5}) annotation(
        Placement(transformation(origin = {-390.4, 2140.4}, extent = {{-141.6, 141.6}, {212.4, -212.4}}, rotation = 90)));
      SMD_MSR_Modelica.Signals.Constants.ConstantNeutronSource constantNeutronSource(S = 0) annotation(
        Placement(transformation(origin = {1416, 2506}, extent = {{-250, -250}, {250, 250}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.FuelChannel R1(vol_FN1 = volF1[1], vol_FN2 = volF2[1], vol_GN = volG[1], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[1], kFN2 = kFN2[1], kG = kHT1[1] + kHT2[1], hAnom = hA[1], kHT_FN1 = kHT1[1], kHT_FN2 = kHT2[1], TF1_0 = TF1_0[1], TF2_0 = TF2_0[1], TG_0 = TG_0[1], regionFlowFrac = flowFracRegions[1], KF = kFuel, Ac = Ac[1], LF1 = LF1[1], LF2 = LF2[1], OutterRegion = false, ArF1 = ArF1[1], ArF2 = ArF2[1], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {-2205.3, -773.305}, extent = {{-504.862, -504.862}, {302.916, 403.889}})));
      SMD_MSR_Modelica.Nuclear.FuelChannel R2(vol_FN1 = volF1[2], vol_FN2 = volF2[2], vol_GN = volG[2], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[2], kFN2 = kFN2[2], kG = kHT1[2] + kHT2[2], hAnom = hA[2], kHT_FN1 = kHT1[2], kHT_FN2 = kHT2[2], TF1_0 = TF1_0[2], TF2_0 = TF2_0[2], TG_0 = TG_0[2], regionFlowFrac = flowFracRegions[2], KF = kFuel, Ac = Ac[2], LF1 = LF1[2], LF2 = LF2[2], OutterRegion = false, ArF1 = ArF1[2], ArF2 = ArF2[2], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {-644.556, -1812.46}, extent = {{-495.416, -495.416}, {297.251, 396.334}})));
      SMD_MSR_Modelica.Nuclear.FuelChannel R3(vol_FN1 = volF1[3], vol_FN2 = volF2[3], vol_GN = volG[3], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[3], kFN2 = kFN2[3], kG = kHT1[3] + kHT2[3], hAnom = hA[3], kHT_FN1 = kHT1[3], kHT_FN2 = kHT2[3], TF1_0 = TF1_0[3], TF2_0 = TF2_0[3], TG_0 = TG_0[3], regionFlowFrac = flowFracRegions[2], KF = kFuel, Ac = Ac[2], LF1 = LF1[3], LF2 = LF2[3], OutterRegion = false, ArF1 = ArF1[3], ArF2 = ArF2[3], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {-684.666, -733.666}, extent = {{-494.446, -494.446}, {296.666, 395.555}})));
      SMD_MSR_Modelica.Nuclear.FuelChannel R4(vol_FN1 = volF1[4], vol_FN2 = volF2[4], vol_GN = volG[4], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[4], kFN2 = kFN2[4], kG = kHT1[4] + kHT2[4], hAnom = hA[4], kHT_FN1 = kHT1[4], kHT_FN2 = kHT2[4], TF1_0 = TF1_0[4], TF2_0 = TF2_0[4], TG_0 = TG_0[4], regionFlowFrac = flowFracRegions[2], KF = kFuel, Ac = Ac[2], LF1 = LF1[4], LF2 = LF2[4], OutterRegion = false, ArF1 = ArF1[4], ArF2 = ArF2[4], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {-677.445, 278.693}, extent = {{-495.832, -495.832}, {297.5, 396.668}})));
      SMD_MSR_Modelica.Nuclear.FuelChannel R5(vol_FN1 = volF1[5], vol_FN2 = volF2[5], vol_GN = volG[5], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[5], kFN2 = kFN2[5], kG = kHT1[5] + kHT2[5], hAnom = hA[5], kHT_FN1 = kHT1[5], kHT_FN2 = kHT2[5], TF1_0 = TF1_0[5], TF2_0 = TF2_0[5], TG_0 = TG_0[5], regionFlowFrac = flowFracRegions[3], KF = kFuel, Ac = Ac[3], LF1 = LF1[5], LF2 = LF2[5], OutterRegion = false, ArF1 = ArF1[5], ArF2 = ArF2[5], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {1000, -1754.96}, extent = {{-492.916, -492.916}, {295.75, 394.334}})));
      SMD_MSR_Modelica.Nuclear.FuelChannel R6(vol_FN1 = volF1[6], vol_FN2 = volF2[6], vol_GN = volG[6], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[6], kFN2 = kFN2[6], kG = kHT1[6] + kHT2[6], hAnom = hA[6], kHT_FN1 = kHT1[6], kHT_FN2 = kHT2[6], TF1_0 = TF1_0[6], TF2_0 = TF2_0[6], TG_0 = TG_0[6], regionFlowFrac = flowFracRegions[3], KF = kFuel, Ac = Ac[3], LF1 = LF1[6], LF2 = LF2[6], OutterRegion = false, ArF1 = ArF1[6], ArF2 = ArF2[6], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {1002.14, -714.194}, extent = {{-490.973, -490.973}, {294.583, 392.778}})));
      SMD_MSR_Modelica.Nuclear.FuelChannel R7(vol_FN1 = volF1[7], vol_FN2 = volF2[7], vol_GN = volG[7], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[7], kFN2 = kFN2[7], kG = kHT1[7] + kHT2[7], hAnom = hA[7], kHT_FN1 = kHT1[7], kHT_FN2 = kHT2[7], TF1_0 = TF1_0[7], TF2_0 = TF2_0[7], TG_0 = TG_0[7], regionFlowFrac = flowFracRegions[3], KF = kFuel, Ac = Ac[3], LF1 = LF1[7], LF2 = LF2[7], OutterRegion = false, ArF1 = ArF1[7], ArF2 = ArF2[7], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {1008.44, 264.791}, extent = {{-490.139, -490.139}, {294.084, 392.111}})));
      SMD_MSR_Modelica.Nuclear.FuelChannel R8(vol_FN1 = volF1[8], vol_FN2 = volF2[8], vol_GN = volG[8], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[8], kFN2 = kFN2[8], kG = kHT1[8] + kHT2[8], hAnom = hA[8], kHT_FN1 = kHT1[8], kHT_FN2 = kHT2[8], TF1_0 = TF1_0[8], TF2_0 = TF2_0[8], TG_0 = TG_0[8], regionFlowFrac = flowFracRegions[4], KF = kFuel, Ac = Ac[4], LF1 = LF1[8], LF2 = LF2[8], OutterRegion = EnableRad, ArF1 = ArF1[8], ArF2 = ArF2[8], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {2652.97, -1026.93}, extent = {{-489.166, -489.166}, {293.5, 391.333}})));
      SMD_MSR_Modelica.Nuclear.FuelChannel R9(vol_FN1 = volF1[9], vol_FN2 = volF2[9], vol_GN = volG[9], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[9], kFN2 = kFN2[9], kG = kHT1[9] + kHT2[9], hAnom = hA[9], kHT_FN1 = kHT1[9], kHT_FN2 = kHT2[9], TF1_0 = TF1_0[9], TF2_0 = TF2_0[9], TG_0 = TG_0[9], regionFlowFrac = flowFracRegions[4], KF = kFuel, Ac = Ac[4], LF1 = LF1[9], LF2 = LF2[9], OutterRegion = EnableRad, ArF1 = ArF1[9], ArF2 = ArF2[9], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {2672.17, 58.3894}, extent = {{-485.834, -485.834}, {291.5, 388.666}})));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF1(a_F = aF, a_G = aG, IF1 = IF1[1], IF2 = IF2[1], IG = IG[1], FuelTempSetPointNode1 = TF1_0[1], FuelTempSetPointNode2 = TF2_0[1], GrapTempSetPoint = TG_0[1]) annotation(
        Placement(transformation(origin = {-1680.4, -798.6}, extent = {{-231.6, -154.4}, {154.4, 231.6}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF2(a_F = aF, a_G = aG, IF1 = IF1[2], IF2 = IF2[2], IG = IG[2], FuelTempSetPointNode1 = TF1_0[2], FuelTempSetPointNode2 = TF2_0[2], GrapTempSetPoint = TG_0[2]) annotation(
        Placement(transformation(origin = {-38.8, -1910.8}, extent = {{-224.4, -149.6}, {149.6, 224.4}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF3(a_F = aF, a_G = aG, IF1 = IF1[3], IF2 = IF2[3], IG = IG[3], FuelTempSetPointNode1 = TF1_0[3], FuelTempSetPointNode2 = TF2_0[3], GrapTempSetPoint = TG_0[3]) annotation(
        Placement(transformation(origin = {-66.4, -802}, extent = {{-225.6, -150.4}, {150.4, 225.6}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF4(a_F = aF, a_G = aG, IF1 = IF1[4], IF2 = IF2[4], IG = IG[4], FuelTempSetPointNode1 = TF1_0[4], FuelTempSetPointNode2 = TF2_0[4], GrapTempSetPoint = TG_0[4]) annotation(
        Placement(transformation(origin = {-24.2, 200.2}, extent = {{-223.2, -148.8}, {148.8, 223.2}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF5(a_F = aF, a_G = aG, IF1 = IF1[5], IF2 = IF2[5], IG = IG[5], FuelTempSetPointNode1 = TF1_0[5], FuelTempSetPointNode2 = TF2_0[5], GrapTempSetPoint = TG_0[5]) annotation(
        Placement(transformation(origin = {1516.2, -1820.4}, extent = {{-218.4, -145.6}, {145.6, 218.4}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF6(a_F = aF, a_G = aG, IF1 = IF1[6], IF2 = IF2[6], IG = IG[6], FuelTempSetPointNode1 = TF1_0[6], FuelTempSetPointNode2 = TF2_0[6], GrapTempSetPoint = TG_0[6]) annotation(
        Placement(transformation(origin = {1518, -760.6}, extent = {{-214.8, -143.2}, {143.2, 214.8}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF7(a_F = aF, a_G = aG, IF1 = IF1[7], IF2 = IF2[7], IG = IG[7], FuelTempSetPointNode1 = TF1_0[7], FuelTempSetPointNode2 = TF2_0[7], GrapTempSetPoint = TG_0[7]) annotation(
        Placement(transformation(origin = {1552, 176.4}, extent = {{-214.8, -143.2}, {143.2, 214.8}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF8(a_F = aF, a_G = aG, IF1 = IF1[8], IF2 = IF2[8], IG = IG[8], FuelTempSetPointNode1 = TF1_0[8], FuelTempSetPointNode2 = TF2_0[8], GrapTempSetPoint = TG_0[8]) annotation(
        Placement(transformation(origin = {3178.2, -1056.4}, extent = {{-206.4, -137.6}, {137.6, 206.4}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF9(a_F = aF, a_G = aG, IF1 = IF1[9], IF2 = IF2[9], IG = IG[9], FuelTempSetPointNode1 = TF1_0[9], FuelTempSetPointNode2 = TF2_0[9], GrapTempSetPoint = TG_0[9]) annotation(
        Placement(transformation(origin = {3186.6, 12}, extent = {{-201.6, -134.4}, {134.4, 201.6}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.SumReactivity sumFB(numInput = 9) annotation(
        Placement(transformation(origin = {2024.34, 1758.74}, extent = {{-159.944, 239.916}, {239.916, -159.944}})));
      SMD_MSR_Modelica.HeatTransport.FlowDistributor flowDistributor(numOutput = 4, freeConvectionFF = freeConvFF, coastDownK = regionCoastDownK, regionTripTime = regionTripTime) annotation(
        Placement(transformation(origin = {1165.6, -3010.4}, extent = {{104.398, 313.195}, {-313.195, -104.398}})));
      SMD_MSR_Modelica.HeatTransport.MixingPot upperPlenum(numStreams = 4, vol = volDotFuel*2, VdotNom = volDotFuel, flowFractionsNom = flowFracRegions, rho = rho_fuel, Cp = cP_fuel, K = kFuel, Ac = 1.6118, L = 0.093874, Ar = 0.4225, e = e, T_0 = 663.53, Tinf = T_inf, EnableRad = EnableRad) annotation(
        Placement(transformation(origin = {1346.18, 1314.32}, extent = {{322.11, 241.582}, {-161.055, -241.582}}, rotation = -90)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFractionIn flowFracIn annotation(
        Placement(transformation(origin = {629, -3219}, extent = {{-103, -103}, {103, 103}}), iconTransformation(origin = {-40, 38}, extent = {{-18, -18}, {18, 18}})));
      input SMD_MSR_Modelica.PortsConnectors.TempIn tempIn annotation(
        Placement(transformation(origin = {-68, -3130}, extent = {{-104, -104}, {104, 104}}), iconTransformation(origin = {-40, -40}, extent = {{-18, -18}, {18, 18}})));
      output SMD_MSR_Modelica.PortsConnectors.TempOut tempOut annotation(
        Placement(transformation(origin = {1323, 1557}, extent = {{82, -82}, {-82, 82}}), iconTransformation(origin = {40, 38}, extent = {{-18, -18}, {18, 18}})));
      output SMD_MSR_Modelica.PortsConnectors.VolumetircPowerOut volumetircPowerOut annotation(
        Placement(transformation(origin = {-1187, 2067}, extent = {{68, 68}, {-68, -68}}), iconTransformation(origin = {41, -41}, extent = {{-17, -17}, {17, 17}})));
      input SMD_MSR_Modelica.PortsConnectors.RealIn realIn annotation(
        Placement(transformation(origin = {40, 2710.4}, extent = {{-72, 72}, {72, -72}}), iconTransformation(origin = {0, 38}, extent = {{-18, -18}, {18, 18}})));
    equation
      connect(realIn, mpke.ReactivityIn) annotation(
        Line(points = {{40, 2696}, {40, 2368}, {678, 2368}}, thickness = 1));
      connect(mpke.n_population, decayHeat.nPop) annotation(
        Line(points = {{357, 2207}, {72.6, 2207}, {72.6, 2208}, {-210, 2208}}, color = {20, 36, 248}, thickness = 1));
      connect(decayHeat.decayHeat_Out, powerblock.decayNP) annotation(
        Line(points = {{-426, 2211}, {-426, 1530}, {-950, 1530}}, color = {0, 225, 255}, thickness = 1));
      connect(mpke.n_population, powerblock.nPop) annotation(
        Line(points = {{357, 2207}, {72, 2207}, {72, 1286}, {-950, 1286}}, color = {20, 36, 248}, thickness = 1));
      connect(powerblock.decayPowerM, volumetircPowerOut) annotation(
        Line(points = {{-1194, 1530}, {-1194, 1837}, {-1187, 1837}, {-1187, 2067}}, color = {220, 138, 221}, thickness = 1));
      connect(constantNeutronSource.neutronEmsRateOut, mpke.S) annotation(
        Line(points = {{1331, 2501}, {995, 2501}, {995, 2261}, {678, 2261}}, color = {191, 64, 188}, thickness = 1));
      connect(upperPlenum.temp_Out, tempOut) annotation(
        Line(points = {{1346, 1306}, {1346, 1318}, {1323, 1318}, {1323, 1557}}, color = {204, 0, 0}, thickness = 1));
      connect(sumFB.reactivityOut, mpke.feedback) annotation(
        Line(points = {{2104, 1959}, {2104, 2154}, {678, 2154}}, color = {78, 154, 6}, thickness = 1));
      connect(R4.fuelNode2, RF4.fuelNode2) annotation(
        Line(points = {{-578, 437}, {-578, 237}, {-210, 237}}, color = {204, 0, 0}, thickness = 1));
      connect(R4.grapNode, RF4.grapNode) annotation(
        Line(points = {{-578, 180}, {-578, 126}, {-210, 126}}, color = {204, 0, 0}, thickness = 1));
      connect(R4.fuelNode1, RF4.fuelNode1) annotation(
        Line(points = {{-578, -68}, {-578, 14}, {-210, 14}}, color = {204, 0, 0}, thickness = 1));
      connect(R3.fuelNode2, RF3.fuelNode2) annotation(
        Line(points = {{-586, -575}, {-586, -764}, {-254, -764}}, color = {204, 0, 0}, thickness = 1));
      connect(R3.fuelNode1, RF3.fuelNode1) annotation(
        Line(points = {{-586, -1080}, {-586, -990}, {-254, -990}}, color = {204, 0, 0}, thickness = 1));
      connect(R3.grapNode, RF3.grapNode) annotation(
        Line(points = {{-586, -833}, {-586, -877}, {-254, -877}}, color = {204, 0, 0}, thickness = 1));
      connect(R2.grapNode, RF2.grapNode) annotation(
        Line(points = {{-545, -1912}, {-545, -1986}, {-226, -1986}}, color = {204, 0, 0}, thickness = 1));
      connect(R2.fuelNode1, RF2.fuelNode1) annotation(
        Line(points = {{-545, -2159}, {-545, -2098}, {-226, -2098}}, color = {204, 0, 0}, thickness = 1));
      connect(R2.fuelNode2, RF2.fuelNode2) annotation(
        Line(points = {{-545, -1654}, {-545, -1873}, {-226, -1873}}, color = {204, 0, 0}, thickness = 1));
      connect(R1.grapNode, RF1.grapNode) annotation(
        Line(points = {{-2104, -874}, {-1917, -874}, {-1917, -876}, {-1873, -876}}, color = {204, 0, 0}, thickness = 1));
      connect(R1.fuelNode1, RF1.fuelNode1) annotation(
        Line(points = {{-2104, -1127}, {-2104, -992}, {-1873, -992}}, color = {204, 0, 0}, thickness = 1));
      connect(R1.fuelNode2, RF1.fuelNode2) annotation(
        Line(points = {{-2104, -612}, {-2104, -760}, {-1873, -760}}, color = {204, 0, 0}, thickness = 1));
      connect(R5.grapNode, RF5.grapNode) annotation(
        Line(points = {{1099, -1854}, {1099, -1893}, {1334, -1893}}, color = {204, 0, 0}, thickness = 1));
      connect(R6.grapNode, RF6.grapNode) annotation(
        Line(points = {{1100, -812}, {1207.5, -812}, {1207.5, -832}, {1339, -832}}, color = {204, 0, 0}, thickness = 1));
      connect(R7.grapNode, RF7.grapNode) annotation(
        Line(points = {{1106, 167}, {1106, 105}, {1373, 105}}, color = {204, 0, 0}, thickness = 1));
      connect(R8.grapNode, RF8.grapNode) annotation(
        Line(points = {{2751, -1125}, {3006, -1125}}, color = {204, 0, 0}, thickness = 1));
      connect(R9.grapNode, RF9.grapNode) annotation(
        Line(points = {{2769, -39}, {2769, -55}, {3019, -55}}, color = {204, 0, 0}, thickness = 1));
      connect(R5.fuelNode1, RF5.fuelNode1) annotation(
        Line(points = {{1099, -2100}, {1099, -2002}, {1334, -2002}}, color = {204, 0, 0}, thickness = 1));
      connect(R5.fuelNode2, RF5.fuelNode2) annotation(
        Line(points = {{1099, -1597}, {1099, -1784}, {1334, -1784}}, color = {204, 0, 0}, thickness = 1));
      connect(R6.fuelNode1, RF6.fuelNode1) annotation(
        Line(points = {{1100, -1058}, {1100, -940}, {1339, -940}}, color = {204, 0, 0}, thickness = 1));
      connect(R6.fuelNode2, RF6.fuelNode2) annotation(
        Line(points = {{1100, -557}, {1100, -725}, {1339, -725}}, color = {204, 0, 0}, thickness = 1));
      connect(R7.fuelNode1, RF7.fuelNode1) annotation(
        Line(points = {{1106, -78}, {1106, -3}, {1373, -3}}, color = {204, 0, 0}, thickness = 1));
      connect(R7.fuelNode2, RF7.fuelNode2) annotation(
        Line(points = {{1106, 422}, {1106, 212}, {1373, 212}}, color = {204, 0, 0}, thickness = 1));
      connect(R9.fuelNode2, RF9.fuelNode2) annotation(
        Line(points = {{2769, 214}, {2769, 46}, {3019, 46}}, color = {204, 0, 0}, thickness = 1));
      connect(R9.fuelNode1, RF9.fuelNode1) annotation(
        Line(points = {{2769, -282}, {2769, -156}, {3019, -156}}, color = {204, 0, 0}, thickness = 1));
      connect(R8.fuelNode2, RF8.fuelNode2) annotation(
        Line(points = {{2751, -870}, {2751, -1022}, {3006, -1022}}, color = {204, 0, 0}, thickness = 1));
      connect(R8.fuelNode1, RF8.fuelNode1) annotation(
        Line(points = {{2751, -1369}, {2751, -1228}, {3006, -1228}}, color = {204, 0, 0}, thickness = 1));
      connect(R3.fuelNode2, R4.temp_In) annotation(
        Line(points = {{-586, -575}, {-492, -575}, {-492, -128}, {-1134, -128}}, color = {204, 0, 0}, thickness = 1));
      connect(R2.fuelNode2, R3.temp_In) annotation(
        Line(points = {{-545, -1654}, {-545, -1340}, {-1140, -1340}, {-1140, -1139}}, color = {204, 0, 0}, thickness = 1));
      connect(R5.fuelNode2, R6.temp_In) annotation(
        Line(points = {{1099, -1597}, {1099, -1333.23}, {550, -1333.23}, {550, -1117}}, color = {204, 0, 0}, thickness = 1));
      connect(R6.fuelNode2, R7.temp_In) annotation(
        Line(points = {{1100, -557}, {994.333, -557}, {994.333, -137}, {558, -137}}, color = {204, 0, 0}, thickness = 1));
      connect(R8.fuelNode2, R9.temp_In) annotation(
        Line(points = {{2751, -870}, {2656.8, -870}, {2656.8, -340}, {2225, -340}}, color = {204, 0, 0}, thickness = 1));
      connect(flowFracIn, flowDistributor.flowFracIn) annotation(
        Line(points = {{629, -3219}, {629, -2922}, {962, -2922}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[1], R1.fuelFlowFraction) annotation(
        Line(points = {{962, -2702}, {-2890, -2702}, {-2890, -975}, {-2660, -975}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[2], R2.fuelFlowFraction) annotation(
        Line(points = {{962, -2702}, {-1238, -2702}, {-1238, -2011}, {-1090, -2011}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[2], R3.fuelFlowFraction) annotation(
        Line(points = {{962, -2702}, {-1238, -2702}, {-1238, -931}, {-1130, -931}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[2], R4.fuelFlowFraction) annotation(
        Line(points = {{962, -2702}, {-1238, -2702}, {-1238, 80}, {-1124, 80}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[3], R5.fuelFlowFraction) annotation(
        Line(points = {{962, -2702}, {356, -2702}, {356, -1932}, {445, -1932}, {445, -1952}, {556, -1952}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[3], R6.fuelFlowFraction) annotation(
        Line(points = {{962, -2702}, {358, -2702}, {358, -911}, {560, -911}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[3], R7.fuelFlowFraction) annotation(
        Line(points = {{962, -2702}, {358, -2702}, {358, 69}, {567, 69}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[4], R8.fuelFlowFraction) annotation(
        Line(points = {{962, -2702}, {1992, -2702}, {1992, -1223}, {2213, -1223}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[4], R9.fuelFlowFraction) annotation(
        Line(points = {{962, -2702}, {1986, -2702}, {1986, -136}, {2235, -136}}, color = {255, 120, 0}, thickness = 1));
      connect(R1.fuelNode2, upperPlenum.temp_In[1]) annotation(
        Line(points = {{-2104, -612}, {-2104, 772}, {1354, 772}, {1354, 992}}, color = {204, 0, 0}, thickness = 1));
      connect(R4.fuelNode2, upperPlenum.temp_In[2]) annotation(
        Line(points = {{-578, 437}, {-578, 766}, {1354, 766}, {1354, 992}}, color = {204, 0, 0}, thickness = 1));
      connect(R7.fuelNode2, upperPlenum.temp_In[3]) annotation(
        Line(points = {{1106, 422}, {1082, 422}, {1082, 918}, {1354, 918}, {1354, 992}}, color = {204, 0, 0}, thickness = 1));
      connect(R9.fuelNode2, upperPlenum.temp_In[4]) annotation(
        Line(points = {{2769, 214}, {2702, 214}, {2702, 912}, {1354, 912}, {1354, 992}}, color = {204, 0, 0}, thickness = 1));
      connect(powerblock.fissionPower, R1.fissionPower) annotation(
        Line(points = {{-1194, 1286}, {-2788, 1286}, {-2788, -541}, {-2670, -541}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R4.fissionPower) annotation(
        Line(points = {{-1194, 1286}, {-1302, 1286}, {-1302, 507}, {-1134, 507}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R3.fissionPower) annotation(
        Line(points = {{-1194, 1286}, {-1302, 1286}, {-1302, -506}, {-1140, -506}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R2.fissionPower) annotation(
        Line(points = {{-1194, 1286}, {-1302, 1286}, {-1302, -1585}, {-1100, -1585}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R7.fissionPower) annotation(
        Line(points = {{-1194, 1286}, {-1302, 1286}, {-1302, 842}, {260, 842}, {260, 490}, {558, 490}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R6.fissionPower) annotation(
        Line(points = {{-1194, 1286}, {-1302, 1286}, {-1302, 838}, {260, 838}, {260, -488}, {550, -488}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R5.fissionPower) annotation(
        Line(points = {{-1194, 1286}, {-1298, 1286}, {-1298, 842}, {260, 842}, {260, -1528}, {547, -1528}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R9.fissionPower) annotation(
        Line(points = {{-1194, 1286}, {-1194, 842}, {1916, 842}, {1916, 282}, {2225, 282}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R8.fissionPower) annotation(
        Line(points = {{-1194, 1286}, {-1194, 838}, {1916, 838}, {1916, -802}, {2203, -802}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.decayPowerM, R1.decayHeat) annotation(
        Line(points = {{-1194, 1530}, {-2886, 1530}, {-2886, -763}, {-2670, -763}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R4.decayHeat) annotation(
        Line(points = {{-1194, 1530}, {-1390, 1530}, {-1390, 289}, {-1134, 289}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R3.decayHeat) annotation(
        Line(points = {{-1194, 1530}, {-1390, 1530}, {-1390, -724}, {-1140, -724}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R2.decayHeat) annotation(
        Line(points = {{-1194, 1530}, {-1394, 1530}, {-1394, -1803}, {-1100, -1803}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, upperPlenum.decayHeat) annotation(
        Line(points = {{-1194, 1530}, {-1194, 992}, {1185, 992}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R7.decayHeat) annotation(
        Line(points = {{-1194, 1530}, {-1394, 1530}, {-1394, 976}, {178, 976}, {178, 275}, {558, 275}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R6.decayHeat) annotation(
        Line(points = {{-1194, 1530}, {-1390, 1530}, {-1390, 970}, {184, 970}, {184, -704}, {550, -704}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R5.decayHeat) annotation(
        Line(points = {{-1194, 1530}, {-1394, 1530}, {-1394, 970}, {188, 970}, {188, -1745}, {547, -1745}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R9.decayHeat) annotation(
        Line(points = {{-1194, 1530}, {-1384, 1530}, {-1384, 690}, {1834, 690}, {1834, 68}, {2225, 68}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R8.decayHeat) annotation(
        Line(points = {{-1194, 1530}, {-1384, 1530}, {-1384, 684}, {1838, 684}, {1838, -1017}, {2203, -1017}}, color = {220, 138, 221}, thickness = 1));
      connect(tempIn, R1.temp_In) annotation(
        Line(points = {{-68, -3130}, {-2814, -3130}, {-2814, -1187}, {-2670, -1187}}, color = {204, 0, 0}, thickness = 1));
      connect(tempIn, R2.temp_In) annotation(
        Line(points = {{-68, -3130}, {-1368, -3130}, {-1368, -2219}, {-1100, -2219}}, color = {204, 0, 0}, thickness = 1));
      connect(tempIn, R5.temp_In) annotation(
        Line(points = {{-68, -3130}, {250, -3130}, {250, -2159}, {547, -2159}}, color = {204, 0, 0}, thickness = 1));
      connect(tempIn, R8.temp_In) annotation(
        Line(points = {{-68, -3130}, {490, -3130}, {490, -2410}, {1894, -2410}, {1894, -1428}, {2203, -1428}}, color = {204, 0, 0}, thickness = 1));
      connect(RF9.feedback, sumFB.reactivityIn[9]) annotation(
        Line(points = {{3220, -55}, {3340, -55}, {3340, 1719}, {2096, 1719}}, color = {78, 154, 6}, thickness = 1));
      connect(RF8.feedback, sumFB.reactivityIn[8]) annotation(
        Line(points = {{3213, -1125}, {3340, -1125}, {3340, 1719}, {2096, 1719}}, color = {78, 154, 6}, thickness = 1));
      connect(RF7.feedback, sumFB.reactivityIn[7]) annotation(
        Line(points = {{1588, 105}, {1752, 105}, {1752, 1719}, {2096, 1719}}, color = {78, 154, 6}, thickness = 1));
      connect(RF6.feedback, sumFB.reactivityIn[6]) annotation(
        Line(points = {{1554, -832}, {1752, -832}, {1752, 1719}, {2096, 1719}}, color = {78, 154, 6}, thickness = 1));
      connect(RF5.feedback, sumFB.reactivityIn[5]) annotation(
        Line(points = {{1553, -1893}, {1752, -1893}, {1752, 1719}, {2096, 1719}}, color = {78, 154, 6}, thickness = 1));
      connect(RF4.feedback, sumFB.reactivityIn[4]) annotation(
        Line(points = {{13, 126}, {132, 126}, {132, 1719}, {2096, 1719}}, color = {78, 154, 6}, thickness = 1));
      connect(RF3.feedback, sumFB.reactivityIn[3]) annotation(
        Line(points = {{-29, -877}, {132, -877}, {132, 1719}, {2096, 1719}}, color = {78, 154, 6}, thickness = 1));
      connect(RF2.feedback, sumFB.reactivityIn[2]) annotation(
        Line(points = {{-1, -1986}, {138, -1986}, {138, 1719}, {2096, 1719}}, color = {78, 154, 6}, thickness = 1));
      connect(RF1.feedback, sumFB.reactivityIn[1]) annotation(
        Line(points = {{-1642, -876}, {-1644, -876}, {-1644, 1719}, {2096, 1719}}, color = {78, 154, 6}, thickness = 1));
      connect(flowFracIn, mpke.fuelFlowFrac) annotation(
        Line(points = {{630, -3218}, {3484, -3218}, {3484, 2047}, {678, 2047}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut, upperPlenum.flowFraction) annotation(
        Line(points = {{962, -2702}, {3742, -2702}, {3742, 992}, {1515, 992}}, color = {245, 121, 0}, thickness = 1));
      annotation(
        Diagram(coordinateSystem(extent = {{-2900, 2880}, {3500, -3560}})),
        Icon(coordinateSystem(extent = {{-60, 60}, {60, -60}}), graphics = {Rectangle(lineThickness = 1.5, extent = {{-60, 60}, {60, -60}}), Text(origin = {3, -1}, extent = {{-53, 31}, {53, -31}}, textString = "Reactor")}));
    end MSRE9R;
  end Components;

  model R1MSREuhx
    parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
    parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 1966.5;
    parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 4.76;
    parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
    parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
    parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 4.76;
    parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
    parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
    SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 324000, hAsNom = 153000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0.01, Tinf = 500, TpIn_0 = 657.2700, TpOut_0 = 632, TsIn_0 = 546.4, TsOut_0 = 579.3) annotation(
      Placement(transformation(origin = {-36.2, 34.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
      Placement(transformation(origin = {115.462, -50.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
    SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0.08, Tinf = 644.4, T_0 = 579.3) annotation(
      Placement(transformation(origin = {-30.8, -57.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
    SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0.08, Tinf = 644.4, T_0 = 546.4) annotation(
      Placement(transformation(origin = {89.2, -0.2667}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = 0.075708, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0.08, Tinf = 644.5, T_0 = 657.2700) annotation(
      Placement(transformation(origin = {-261, 36.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
    SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.066953, volFracNode = 0.01, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
      Placement(transformation(origin = {-149.6, 13.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
    SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 0.6564, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0.08, Tinf = 644.4, T_0 = 632) annotation(
      Placement(transformation(origin = {-182.4, -135.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
    SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.066953, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
      Placement(transformation(origin = {-378.4, 5.2}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
    MSRE.Components.MSRE1R core1R(numGroups = 6, addRho0 = true, EnableRad = false, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, a_F = -8.71E-05, a_G = -6.66E-05, n_0 = 1, rho_fuel = rhoFuel, rho_grap = rhoGrap, cP_fuel = scpFuel, cP_grap = scpGrap, volDotFuel = volDotFuel, T_inf = 500, kFuel = kFuel, TF1_0 = 644.5, TF2_0 = 657.2700, TG_0 = 675.8561) annotation(
      Placement(transformation(origin = {-495.667, -31}, extent = {{-172.333, -141}, {-78.3333, -47}})));
    SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
      Placement(transformation(origin = {43, -145}, extent = {{-27, -27}, {27, 27}})));
    SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
      Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
    SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
      Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
    SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
      Placement(transformation(origin = {73, -225}, extent = {{-27, -27}, {27, 27}})));
    SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 0}) annotation(
      Placement(transformation(origin = {-662.198, 95.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
  equation
    connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
      Line(points = {{-77, 22}, {-94, 22}, {-94, -34}, {-53, -34}}, color = {204, 0, 0}, thickness = 1));
    connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
      Line(points = {{-17, -34}, {28.5, -34}, {28.5, -50}, {74, -50}}, color = {204, 0, 0}, thickness = 1));
    connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
      Line(points = {{129, -50}, {140, -50}, {140, 22}, {109, 22}}, color = {204, 0, 0}, thickness = 1));
    connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
      Line(points = {{79, 22}, {30, 22}}, color = {204, 0, 0}, thickness = 1));
    connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
      Line(points = {{-241, 49}, {-207, 49}, {-207, 45}, {-177, 45}}, color = {204, 0, 0}, thickness = 1));
    connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
      Line(points = {{-135, 45}, {-104, 45}, {-104, 47}, {-77, 47}}, color = {204, 0, 0}, thickness = 1));
    connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
      Line(points = {{30, 47}, {156, 47}, {156, -167}, {-151, -167}}, color = {204, 0, 0}, thickness = 1));
    connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
      Line(points = {{-359, 38}, {-325.5, 38}, {-325.5, 49}, {-294, 49}}, color = {204, 0, 0}, thickness = 1));
    connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
      Line(points = {{-715, -188}, {-715, 38}, {-411, 38}}, color = {204, 0, 0}, thickness = 1));
    connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
      Line(points = {{-201, -167}, {-201, -110.5}, {-235, -110.5}, {-235, -110}, {-543, -110}, {-543, -249}, {-778, -249}}, color = {204, 0, 0}, thickness = 1));
    connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
      Line(points = {{-714, -250}, {-448, -250}, {-448, 51}, {-411, 51}}, color = {220, 138, 221}, thickness = 1));
    connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
      Line(points = {{-714, -250}, {-448, -250}, {-448, 74}, {-294, 74}}, color = {220, 138, 221}, thickness = 1));
    connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
      Line(points = {{-714, -250}, {-448, -250}, {-448, 55}, {-177, 55}}, color = {220, 138, 221}, thickness = 1));
    connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
      Line(points = {{-714, -250}, {-714, 114}, {-23.5, 114}, {-23.5, 56}}, color = {220, 138, 221}, thickness = 1));
    connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
      Line(points = {{-714, -250}, {-714, -179}, {-151, -179}}, color = {220, 138, 221}, thickness = 1));
    connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
      Line(points = {{41, -129}, {-74, -129}, {-74, -25}, {-53, -25}}, color = {220, 138, 221}, thickness = 1));
    connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
      Line(points = {{41, -129}, {148, -129}, {148, 29}, {109, 29}}, color = {220, 138, 221}, thickness = 1));
    connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
      Line(points = {{-550, 92}, {-550, 35}, {-778, 35}, {-778, -188}}, color = {245, 121, 0}, thickness = 1));
    connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
      Line(points = {{-550, 92}, {-468, 92}, {-468, 25}, {-411, 25}}, color = {245, 121, 0}, thickness = 1));
    connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
      Line(points = {{-550, 92}, {-328, 92}, {-328, 24}, {-294, 24}}, color = {245, 121, 0}, thickness = 1));
    connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
      Line(points = {{-550, 92}, {-206, 92}, {-206, 36}, {-177, 36}}, color = {245, 121, 0}, thickness = 1));
    connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
      Line(points = {{-550, 92}, {-62, 92}, {-62, 56}, {-59, 56}}, color = {245, 121, 0}, thickness = 1));
    connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
      Line(points = {{-550, 92}, {-108, 92}, {-108, -154}, {-151, -154}}, color = {245, 121, 0}, thickness = 1));
    connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
      Line(points = {{175, 92}, {174, 92}, {174, 16}, {109, 16}}, color = {245, 121, 0}, thickness = 1));
    connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
      Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {12, -6}, {12, 13}}, color = {245, 121, 0}, thickness = 1));
    connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
      Line(points = {{175, 92}, {174, 92}, {174, -23}, {74, -23}}, color = {245, 121, 0}, thickness = 1));
    connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
      Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-53, -82}, {-53, -43}}, color = {245, 121, 0}, thickness = 1));
  connect(externalReact.step, core1R.realIn) annotation(
      Line(points = {{-660, 106}, {-660, -50}, {-746, -50}, {-746, -188}}, thickness = 1));
  connect(uhxDemand.realOut, uhx.powDemand) annotation(
      Line(points = {{74, -212}, {74, -78}}, thickness = 1));
    annotation(
      Diagram(coordinateSystem(extent = {{-800, 160}, {220, -300}})));
  end R1MSREuhx;

  model R9MSREuhx
    parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
    parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2000;
    parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
    parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
    parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
    parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
    parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
    parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
    SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 648000, hAsNom = 306000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0, Tinf = 663.53, TpIn_0 = 663.53, TpOut_0 = 632, TsIn_0 = 554.4300, TsOut_0 = 579.3) annotation(
      Placement(transformation(origin = {4.6, 24}, extent = {{-93.6, -46.8}, {140.4, 46.8}})));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
      Placement(transformation(origin = {80.7953, -102.346}, extent = {{-45.5285, -34.1463}, {22.7642, 34.1463}})));
    SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0, Tinf = 644.4, T_0 = 579.3) annotation(
      Placement(transformation(origin = {-32.1333, -61}, extent = {{-49.2, -16.4}, {32.7999, -65.6}})));
    SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0, Tinf = 644.4, T_0 = 546.4) annotation(
      Placement(transformation(origin = {164.4, -60.6667}, extent = {{-50, -16.6667}, {33.3333, -66.6667}})));
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = volDotFuel*1, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
      Placement(transformation(origin = {-259, 34.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
    SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
      Placement(transformation(origin = {-143.467, -2.53333}, extent = {{-59.2, 19.7333}, {39.4667, 78.9333}})));
    SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 8.67*volDotFuel, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0, Tinf = 632, T_0 = 632) annotation(
      Placement(transformation(origin = {-162, -104.8}, extent = {{-52.8, 17.6}, {35.2, 70.4}}, rotation = 180)));
    SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
      Placement(transformation(origin = {-374.4, 1.86667}, extent = {{-55.6, 18.5333}, {37.0667, 74.1333}})));
    SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
      Placement(transformation(origin = {-68, -180}, extent = {{-17, -17}, {17, 17}})));
    SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 1000000, freeConvFF = 0) annotation(
      Placement(transformation(origin = {-528, 122}, extent = {{-24, -32}, {24, 16}})));
    SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
      Placement(transformation(origin = {-33, -39.333}, extent = {{-17, -22.6667}, {17, 11.3333}})));
    SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
      Placement(transformation(origin = {36, -200}, extent = {{-17, -17}, {17, 17}})));
    SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 0}) annotation(
      Placement(transformation(origin = {-496.198, -16.1983}, extent = {{-17.7983, 17.7983}, {17.7983, -17.7983}}, rotation = -0)));
    MSRE.Components.MSRE9R msre9r(numGroups = 6, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, n_0 = 1, aF = -8.71E-05, aG = -6.66E-05, addRho0 = true, volF1 = {6.4392e-03, 2.1834e-02, 1.1940e-02, 1.4925e-02, 3.6929e-02, 2.0170e-02, 2.5245e-02, 1.0149e-01, 6.8869e-02}, volF2 = {6.7377e-03, 1.4883e-02, 1.1940e-02, 2.9083e-02, 2.5245e-02, 2.0170e-02, 4.9125e-02, 5.9019e-02, 1.1556e-01}, volG = {0.038208, 0.115387, 0.087659, 0.112221, 0.195186, 0.148354, 0.189837, 0.524644, 0.514219}, volUP = volDotFuel*2, hA = {705.60, 2167.20, 1620.00, 2113.20, 3558.60, 2745.00, 3573.00, 9801.00, 9648.00}, kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, TF1_0 = {641.0644, 636.1926, 662.2593, 692.4323, 634.6308, 650.7339, 669.3933, 633.4354, 638.9545}, TF2_0 = {651.9893, 647.9801, 677.5167, 699.8892, 641.9048, 660.1706, 673.9817, 635.8835, 640.4922}, TG_0 = {664.0462, 653.7083, 692.6790, 714.5485, 648.7869, 674.7077, 686.8285, 640.8657, 647.7290}, Tmix_0 = 663.53, IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, flowFracRegions = {0.061409895150439, 0.138549617100260, 0.234231239311780, 0.565809248437521}, rho_fuel = rhoFuel, cP_fuel = scpFuel, kFuel = kFuel, volDotFuel = volDotFuel, rho_grap = rhoGrap, cP_grap = scpGrap, regionTripTime = {200000, 200000, 200000, 200000}, regionCoastDownK = 0.02, freeConvFF = 0, EnableRad = false, T_inf = 663.53, LF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, LF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, Ac = {7.6800e-03, 4.4928e-02, 1.0368e-01, 1.7222e-01}, ArF1 = {0.5301, 0.5692, 0.3113, 0.3891, 0.9719, 0.5308, 0.6644, 2.0856, 1.4152}, ArF2 = {0.5546, 0.3880, 0.3113, 0.7582, 0.6644, 0.5308, 1.2929, 1.2128, 2.3747}, e = 0) annotation(
      Placement(transformation(origin = {-364.333, -15}, extent = {{-179.667, -147}, {-81.6667, -49}})));
  equation
    connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
      Line(points = {{-71, 1}, {-94, 1}, {-94, -102}, {-68, -102}}, color = {204, 0, 0}, thickness = 1));
    connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
      Line(points = {{-13, -102}, {35, -102}}, color = {204, 0, 0}, thickness = 1));
    connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
      Line(points = {{81, -102}, {128, -102}}, color = {204, 0, 0}, thickness = 1));
    connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
      Line(points = {{-239, 47}, {-186, 47}}, color = {204, 0, 0}, thickness = 1));
    connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
      Line(points = {{-120, 47}, {-71, 47}}, color = {204, 0, 0}, thickness = 1));
    connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
      Line(points = {{-353, 48}, {-328.5, 48}, {-328.5, 47}, {-292, 47}}, color = {204, 0, 0}, thickness = 1));
    connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
      Line(points = {{-528, 90}, {-468, 90}, {-468, 34}, {-415, 34}}, color = {245, 121, 0}, thickness = 1));
    connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
      Line(points = {{-528, 90}, {-328, 90}, {-328, 22}, {-292, 22}}, color = {245, 121, 0}, thickness = 1));
    connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
      Line(points = {{-528, 90}, {-206, 90}, {-206, 32}, {-186, 32}}, color = {245, 121, 0}, thickness = 1));
    connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
      Line(points = {{-528, 90}, {-299.5, 90}, {-299.5, 63}, {-38, 63}}, color = {245, 121, 0}, thickness = 1));
    connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
      Line(points = {{-528, 90}, {-108, 90}, {-108, -136}, {-124, -136}}, color = {245, 121, 0}, thickness = 1));
    connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
      Line(points = {{-462, -146}, {-293, -146}, {-293, -162}, {-124, -162}}, color = {220, 138, 221}, thickness = 1));
    connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
      Line(points = {{-462, -146}, {-434, -146}, {-434, 62}, {-415, 62}}, color = {220, 138, 221}, thickness = 1));
    connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
      Line(points = {{-462, -146}, {-318, -146}, {-318, 72}, {-292, 72}}, color = {220, 138, 221}, thickness = 1));
    connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
      Line(points = {{-462, -146}, {-200, -146}, {-200, 62}, {-186, 62}}, color = {220, 138, 221}, thickness = 1));
    connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
      Line(points = {{-462, -82}, {-460, -82}, {-460, 48}, {-415, 48}}, color = {204, 0, 0}, thickness = 1));
    connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
      Line(points = {{-528, 90}, {-528, -82}}, color = {245, 121, 0}, thickness = 1));
    connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
      Line(points = {{-462, -146}, {-434, -146}, {-434, 118}, {28, 118}, {28, 63}}, color = {220, 138, 221}, thickness = 1));
    connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
      Line(points = {{-183, -149}, {-183, -178}, {-528, -178}, {-528, -146}}, color = {204, 0, 0}, thickness = 1));
    connect(heatExchanger.T_in_sFluid, pipeUHXtoHX.PiTempOut) annotation(
      Line(points = {{128, 0}, {200, 0}, {200, -102}, {184, -102}}, color = {204, 0, 0}, thickness = 1));
    connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
      Line(points = {{-33, -62}, {-68, -62}, {-68, -90}}, color = {245, 121, 0}, thickness = 1));
    connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
      Line(points = {{-33, -62}, {-33, -80}, {35, -80}}, color = {245, 121, 0}, thickness = 1));
    connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
      Line(points = {{-33, -62}, {-33, -90}, {128, -90}}, color = {245, 121, 0}, thickness = 1));
    connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
      Line(points = {{-33, -62}, {-33, -16}, {94, -16}}, color = {245, 121, 0}, thickness = 1));
    connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
      Line(points = {{-70, -170}, {-70, -168}, {-68, -168}, {-68, -114}}, color = {220, 138, 221}, thickness = 1));
    connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
      Line(points = {{-70, -170}, {-70, -169}, {128, -169}, {128, -115}}, color = {220, 138, 221}, thickness = 1));
    connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
      Line(points = {{128, 48}, {224, 48}, {224, -149}, {-124, -149}}, color = {204, 0, 0}, thickness = 1));
  connect(externalReact.step, msre9r.realIn) annotation(
      Line(points = {{-496, -22}, {-494, -22}, {-494, -82}}, thickness = 1));
  connect(uhxDemand.realOut, uhx.powDemand) annotation(
      Line(points = {{36, -192}, {36, -126}}, thickness = 1));
    annotation(
      Diagram(coordinateSystem(extent = {{-580, 160}, {240, -240}})));
  end R9MSREuhx;

  package Transients
   
package R1fullSteps
      model R1MSRE2dol
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
        parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 1966.5;
        parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 4.76;
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
        parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
        parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 4.76;
        parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
        parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
        SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 324000, hAsNom = 153000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0.01, Tinf = 500, TpIn_0 = 657.2700, TpOut_0 = 632, TsIn_0 = 546.4, TsOut_0 = 579.3) annotation(
          Placement(transformation(origin = {-36.2, 34.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
        SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
          Placement(transformation(origin = {115.462, -50.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0.08, Tinf = 644.4, T_0 = 579.3) annotation(
          Placement(transformation(origin = {-30.8, -57.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0.08, Tinf = 644.4, T_0 = 546.4) annotation(
          Placement(transformation(origin = {89.2, 21.7333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = 0.075708, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0.08, Tinf = 644.5, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-261, 36.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.066953, volFracNode = 0.01, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-145.6, 5.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 0.6564, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0.08, Tinf = 644.4, T_0 = 632) annotation(
          Placement(transformation(origin = {-182.4, -135.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.066953, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-378.4, 5.2}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        MSRE.Components.MSRE1R core1R(numGroups = 6, addRho0 = true, EnableRad = false, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, a_F = -8.71E-05, a_G = -6.66E-05, n_0 = 1, rho_fuel = rhoFuel, rho_grap = rhoGrap, cP_fuel = scpFuel, cP_grap = scpGrap, volDotFuel = volDotFuel, T_inf = 500, kFuel = kFuel, TF1_0 = 644.5, TF2_0 = 657.2700, TG_0 = 675.8561) annotation(
          Placement(transformation(origin = {-495.667, -31}, extent = {{-172.333, -141}, {-78.3333, -47}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
          Placement(transformation(origin = {43, -145}, extent = {{-27, -27}, {27, 27}})));
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
          Placement(transformation(origin = {73, -225}, extent = {{-27, -27}, {27, 27}})));
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 829.172}) annotation(
          Placement(transformation(origin = {-662.198, 95.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-77, 22}, {-94, 22}, {-94, -34}, {-53, -34}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-17, -34}, {28.5, -34}, {28.5, -50}, {74, -50}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{129, -50}, {140, -50}, {140, 44}, {112, 44}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{76, 44}, {53.5, 44}, {53.5, 22}, {30, 22}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-241, 49}, {-207, 49}, {-207, 37}, {-178, 37}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-126, 37}, {-126, 47}, {-77, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{30, 47}, {156, 47}, {156, -167}, {-151, -167}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-359, 38}, {-325.5, 38}, {-325.5, 49}, {-294, 49}}, color = {204, 0, 0}, thickness = 1));
        connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-715, -188}, {-713, 38}, {-411, 38}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
          Line(points = {{-201, -167}, {-201, -110.5}, {-235, -110.5}, {-235, -110}, {-543, -110}, {-545, -249}, {-778, -249}}, color = {204, 0, 0}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 51}, {-411, 51}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 74}, {-294, 74}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 50}, {-178, 50}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-714, -250}, {-714, 114}, {-23.5, 114}, {-23.5, 56}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-712, -179}, {-151, -179}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{41, -129}, {-74, -129}, {-74, -25}, {-53, -25}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{41, -129}, {148, -129}, {148, 53}, {112, 53}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
          Line(points = {{-550, 92}, {-550, 35}, {-778, 35}, {-778, -188}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-550, 92}, {-468, 92}, {-468, 25}, {-411, 25}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-550, 92}, {-328, 92}, {-328, 24}, {-294, 24}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-550, 92}, {-206, 92}, {-206, 24}, {-178, 24}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-550, 92}, {-62, 92}, {-62, 56}, {-59, 56}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-550, 92}, {-108, 92}, {-108, -154}, {-151, -154}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, 35}, {112, 35}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {12, -6}, {12, 13}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, -23}, {74, -23}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-53, -82}, {-53, -43}}, color = {245, 121, 0}, thickness = 1));
      connect(externalReact.step, core1R.realIn) annotation(
          Line(points = {{-660, 106}, {-660, -50}, {-746, -50}, {-746, -188}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{74, -212}, {74, -78}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-800, 160}, {220, -300}})));
      end R1MSRE2dol;
      
      model R1MSRE1dol
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
        parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 1966.5;
        parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 4.76;
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
        parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
        parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 4.76;
        parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
        parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
        SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 324000, hAsNom = 153000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0.01, Tinf = 500, TpIn_0 = 657.2700, TpOut_0 = 632, TsIn_0 = 546.4, TsOut_0 = 579.3) annotation(
          Placement(transformation(origin = {-36.2, 34.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
        SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
          Placement(transformation(origin = {115.462, -50.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0.08, Tinf = 644.4, T_0 = 579.3) annotation(
          Placement(transformation(origin = {-30.8, -57.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0.08, Tinf = 644.4, T_0 = 546.4) annotation(
          Placement(transformation(origin = {89.2, 21.7333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = 0.075708, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0.08, Tinf = 644.5, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-261, 36.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.066953, volFracNode = 0.01, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-145.6, 5.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 0.6564, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0.08, Tinf = 644.4, T_0 = 632) annotation(
          Placement(transformation(origin = {-182.4, -135.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.066953, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-378.4, 5.2}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        MSRE.Components.MSRE1R core1R(numGroups = 6, addRho0 = true, EnableRad = false, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, a_F = -8.71E-05, a_G = -6.66E-05, n_0 = 1, rho_fuel = rhoFuel, rho_grap = rhoGrap, cP_fuel = scpFuel, cP_grap = scpGrap, volDotFuel = volDotFuel, T_inf = 500, kFuel = kFuel, TF1_0 = 644.5, TF2_0 = 657.2700, TG_0 = 675.8561) annotation(
          Placement(transformation(origin = {-495.667, -31}, extent = {{-172.333, -141}, {-78.3333, -47}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
          Placement(transformation(origin = {43, -145}, extent = {{-27, -27}, {27, 27}})));
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
          Placement(transformation(origin = {73, -225}, extent = {{-27, -27}, {27, 27}})));
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 414.586}) annotation(
          Placement(transformation(origin = {-662.198, 95.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-77, 22}, {-94, 22}, {-94, -34}, {-53, -34}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-17, -34}, {28.5, -34}, {28.5, -50}, {74, -50}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{129, -50}, {140, -50}, {140, 44}, {112, 44}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{76, 44}, {53.5, 44}, {53.5, 22}, {30, 22}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-241, 49}, {-207, 49}, {-207, 37}, {-178, 37}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-126, 37}, {-126, 47}, {-77, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{30, 47}, {156, 47}, {156, -167}, {-151, -167}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-359, 38}, {-325.5, 38}, {-325.5, 49}, {-294, 49}}, color = {204, 0, 0}, thickness = 1));
        connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-715, -188}, {-713, 38}, {-411, 38}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
          Line(points = {{-201, -167}, {-201, -110.5}, {-235, -110.5}, {-235, -110}, {-543, -110}, {-545, -249}, {-778, -249}}, color = {204, 0, 0}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 51}, {-411, 51}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 74}, {-294, 74}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 50}, {-178, 50}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-714, -250}, {-714, 114}, {-23.5, 114}, {-23.5, 56}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-712, -179}, {-151, -179}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{41, -129}, {-74, -129}, {-74, -25}, {-53, -25}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{41, -129}, {148, -129}, {148, 53}, {112, 53}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
          Line(points = {{-550, 92}, {-550, 35}, {-778, 35}, {-778, -188}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-550, 92}, {-468, 92}, {-468, 25}, {-411, 25}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-550, 92}, {-328, 92}, {-328, 24}, {-294, 24}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-550, 92}, {-206, 92}, {-206, 24}, {-178, 24}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-550, 92}, {-62, 92}, {-62, 56}, {-59, 56}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-550, 92}, {-108, 92}, {-108, -154}, {-151, -154}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, 35}, {112, 35}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {12, -6}, {12, 13}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, -23}, {74, -23}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-53, -82}, {-53, -43}}, color = {245, 121, 0}, thickness = 1));
      connect(externalReact.step, core1R.realIn) annotation(
          Line(points = {{-660, 106}, {-660, -50}, {-746, -50}, {-746, -188}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{74, -212}, {74, -78}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-800, 160}, {220, -300}})));
      end R1MSRE1dol;
      
      model R1MSREhalfDol
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
        parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 1966.5;
        parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 4.76;
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
        parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
        parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 4.76;
        parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
        parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
        SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 324000, hAsNom = 153000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0.01, Tinf = 500, TpIn_0 = 657.2700, TpOut_0 = 632, TsIn_0 = 546.4, TsOut_0 = 579.3) annotation(
          Placement(transformation(origin = {-36.2, 34.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
        SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
          Placement(transformation(origin = {115.462, -50.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0.08, Tinf = 644.4, T_0 = 579.3) annotation(
          Placement(transformation(origin = {-30.8, -57.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0.08, Tinf = 644.4, T_0 = 546.4) annotation(
          Placement(transformation(origin = {89.2, 21.7333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = 0.075708, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0.08, Tinf = 644.5, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-261, 36.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.066953, volFracNode = 0.01, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-145.6, 5.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 0.6564, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0.08, Tinf = 644.4, T_0 = 632) annotation(
          Placement(transformation(origin = {-182.4, -135.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.066953, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-378.4, 5.2}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        MSRE.Components.MSRE1R core1R(numGroups = 6, addRho0 = true, EnableRad = false, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, a_F = -8.71E-05, a_G = -6.66E-05, n_0 = 1, rho_fuel = rhoFuel, rho_grap = rhoGrap, cP_fuel = scpFuel, cP_grap = scpGrap, volDotFuel = volDotFuel, T_inf = 500, kFuel = kFuel, TF1_0 = 644.5, TF2_0 = 657.2700, TG_0 = 675.8561) annotation(
          Placement(transformation(origin = {-495.667, -31}, extent = {{-172.333, -141}, {-78.3333, -47}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
          Placement(transformation(origin = {43, -145}, extent = {{-27, -27}, {27, 27}})));
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
          Placement(transformation(origin = {73, -225}, extent = {{-27, -27}, {27, 27}})));
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 207.293}) annotation(
          Placement(transformation(origin = {-662.198, 95.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-77, 22}, {-94, 22}, {-94, -34}, {-53, -34}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-17, -34}, {28.5, -34}, {28.5, -50}, {74, -50}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{129, -50}, {140, -50}, {140, 44}, {112, 44}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{76, 44}, {53.5, 44}, {53.5, 22}, {30, 22}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-241, 49}, {-207, 49}, {-207, 37}, {-178, 37}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-126, 37}, {-126, 47}, {-77, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{30, 47}, {156, 47}, {156, -167}, {-151, -167}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-359, 38}, {-325.5, 38}, {-325.5, 49}, {-294, 49}}, color = {204, 0, 0}, thickness = 1));
        connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-715, -188}, {-713, 38}, {-411, 38}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
          Line(points = {{-201, -167}, {-201, -110.5}, {-235, -110.5}, {-235, -110}, {-543, -110}, {-545, -249}, {-778, -249}}, color = {204, 0, 0}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 51}, {-411, 51}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 74}, {-294, 74}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 50}, {-178, 50}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-714, -250}, {-714, 114}, {-23.5, 114}, {-23.5, 56}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-712, -179}, {-151, -179}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{41, -129}, {-74, -129}, {-74, -25}, {-53, -25}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{41, -129}, {148, -129}, {148, 53}, {112, 53}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
          Line(points = {{-550, 92}, {-550, 35}, {-778, 35}, {-778, -188}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-550, 92}, {-468, 92}, {-468, 25}, {-411, 25}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-550, 92}, {-328, 92}, {-328, 24}, {-294, 24}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-550, 92}, {-206, 92}, {-206, 24}, {-178, 24}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-550, 92}, {-62, 92}, {-62, 56}, {-59, 56}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-550, 92}, {-108, 92}, {-108, -154}, {-151, -154}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, 35}, {112, 35}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {12, -6}, {12, 13}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, -23}, {74, -23}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-53, -82}, {-53, -43}}, color = {245, 121, 0}, thickness = 1));
      connect(externalReact.step, core1R.realIn) annotation(
          Line(points = {{-660, 106}, {-660, -50}, {-746, -50}, {-746, -188}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{74, -212}, {74, -78}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-800, 160}, {220, -300}})));
      end R1MSREhalfDol;
      
      model R1MSREpOneDol
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
        parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 1966.5;
        parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 4.76;
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
        parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
        parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 4.76;
        parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
        parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
        SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 324000, hAsNom = 153000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0.01, Tinf = 500, TpIn_0 = 657.2700, TpOut_0 = 632, TsIn_0 = 546.4, TsOut_0 = 579.3) annotation(
          Placement(transformation(origin = {-36.2, 34.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
        SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
          Placement(transformation(origin = {115.462, -50.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0.08, Tinf = 644.4, T_0 = 579.3) annotation(
          Placement(transformation(origin = {-30.8, -57.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0.08, Tinf = 644.4, T_0 = 546.4) annotation(
          Placement(transformation(origin = {89.2, 21.7333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = 0.075708, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0.08, Tinf = 644.5, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-261, 36.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.066953, volFracNode = 0.01, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-145.6, 5.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 0.6564, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0.08, Tinf = 644.4, T_0 = 632) annotation(
          Placement(transformation(origin = {-182.4, -135.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.066953, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-378.4, 5.2}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        MSRE.Components.MSRE1R core1R(numGroups = 6, addRho0 = true, EnableRad = false, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, a_F = -8.71E-05, a_G = -6.66E-05, n_0 = 1, rho_fuel = rhoFuel, rho_grap = rhoGrap, cP_fuel = scpFuel, cP_grap = scpGrap, volDotFuel = volDotFuel, T_inf = 500, kFuel = kFuel, TF1_0 = 644.5, TF2_0 = 657.2700, TG_0 = 675.8561) annotation(
          Placement(transformation(origin = {-495.667, -31}, extent = {{-172.333, -141}, {-78.3333, -47}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
          Placement(transformation(origin = {43, -145}, extent = {{-27, -27}, {27, 27}})));
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
          Placement(transformation(origin = {73, -225}, extent = {{-27, -27}, {27, 27}})));
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 41.459}) annotation(
          Placement(transformation(origin = {-662.198, 95.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-77, 22}, {-94, 22}, {-94, -34}, {-53, -34}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-17, -34}, {28.5, -34}, {28.5, -50}, {74, -50}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{129, -50}, {140, -50}, {140, 44}, {112, 44}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{76, 44}, {53.5, 44}, {53.5, 22}, {30, 22}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-241, 49}, {-207, 49}, {-207, 37}, {-178, 37}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-126, 37}, {-126, 47}, {-77, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{30, 47}, {156, 47}, {156, -167}, {-151, -167}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-359, 38}, {-325.5, 38}, {-325.5, 49}, {-294, 49}}, color = {204, 0, 0}, thickness = 1));
        connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-715, -188}, {-713, 38}, {-411, 38}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
          Line(points = {{-201, -167}, {-201, -110.5}, {-235, -110.5}, {-235, -110}, {-543, -110}, {-545, -249}, {-778, -249}}, color = {204, 0, 0}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 51}, {-411, 51}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 74}, {-294, 74}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 50}, {-178, 50}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-714, -250}, {-714, 114}, {-23.5, 114}, {-23.5, 56}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-712, -179}, {-151, -179}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{41, -129}, {-74, -129}, {-74, -25}, {-53, -25}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{41, -129}, {148, -129}, {148, 53}, {112, 53}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
          Line(points = {{-550, 92}, {-550, 35}, {-778, 35}, {-778, -188}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-550, 92}, {-468, 92}, {-468, 25}, {-411, 25}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-550, 92}, {-328, 92}, {-328, 24}, {-294, 24}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-550, 92}, {-206, 92}, {-206, 24}, {-178, 24}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-550, 92}, {-62, 92}, {-62, 56}, {-59, 56}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-550, 92}, {-108, 92}, {-108, -154}, {-151, -154}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, 35}, {112, 35}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {12, -6}, {12, 13}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, -23}, {74, -23}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-53, -82}, {-53, -43}}, color = {245, 121, 0}, thickness = 1));
      connect(externalReact.step, core1R.realIn) annotation(
          Line(points = {{-660, 106}, {-660, -50}, {-746, -50}, {-746, -188}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{74, -212}, {74, -78}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-800, 160}, {220, -300}})));
      end R1MSREpOneDol;
      
      model R1MSRE100pcm
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
        parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 1966.5;
        parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 4.76;
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
        parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
        parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 4.76;
        parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
        parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
        SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 324000, hAsNom = 153000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0.01, Tinf = 500, TpIn_0 = 657.2700, TpOut_0 = 632, TsIn_0 = 546.4, TsOut_0 = 579.3) annotation(
          Placement(transformation(origin = {-36.2, 34.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
        SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
          Placement(transformation(origin = {115.462, -50.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0.08, Tinf = 644.4, T_0 = 579.3) annotation(
          Placement(transformation(origin = {-30.8, -57.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0.08, Tinf = 644.4, T_0 = 546.4) annotation(
          Placement(transformation(origin = {89.2, 21.7333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = 0.075708, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0.08, Tinf = 644.5, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-261, 36.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.066953, volFracNode = 0.01, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-145.6, 5.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 0.6564, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0.08, Tinf = 644.4, T_0 = 632) annotation(
          Placement(transformation(origin = {-182.4, -135.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.066953, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-378.4, 5.2}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        MSRE.Components.MSRE1R core1R(numGroups = 6, addRho0 = true, EnableRad = false, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, a_F = -8.71E-05, a_G = -6.66E-05, n_0 = 1, rho_fuel = rhoFuel, rho_grap = rhoGrap, cP_fuel = scpFuel, cP_grap = scpGrap, volDotFuel = volDotFuel, T_inf = 500, kFuel = kFuel, TF1_0 = 644.5, TF2_0 = 657.2700, TG_0 = 675.8561) annotation(
          Placement(transformation(origin = {-495.667, -31}, extent = {{-172.333, -141}, {-78.3333, -47}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
          Placement(transformation(origin = {43, -145}, extent = {{-27, -27}, {27, 27}})));
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
          Placement(transformation(origin = {73, -225}, extent = {{-27, -27}, {27, 27}})));
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 100}) annotation(
          Placement(transformation(origin = {-662.198, 95.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-77, 22}, {-94, 22}, {-94, -34}, {-53, -34}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-17, -34}, {28.5, -34}, {28.5, -50}, {74, -50}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{129, -50}, {140, -50}, {140, 44}, {112, 44}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{76, 44}, {53.5, 44}, {53.5, 22}, {30, 22}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-241, 49}, {-207, 49}, {-207, 37}, {-178, 37}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-126, 37}, {-126, 47}, {-77, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{30, 47}, {156, 47}, {156, -167}, {-151, -167}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-359, 38}, {-325.5, 38}, {-325.5, 49}, {-294, 49}}, color = {204, 0, 0}, thickness = 1));
        connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-715, -188}, {-713, 38}, {-411, 38}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
          Line(points = {{-201, -167}, {-201, -110.5}, {-235, -110.5}, {-235, -110}, {-543, -110}, {-545, -249}, {-778, -249}}, color = {204, 0, 0}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 51}, {-411, 51}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 74}, {-294, 74}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 50}, {-178, 50}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-714, -250}, {-714, 114}, {-23.5, 114}, {-23.5, 56}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-712, -179}, {-151, -179}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{41, -129}, {-74, -129}, {-74, -25}, {-53, -25}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{41, -129}, {148, -129}, {148, 53}, {112, 53}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
          Line(points = {{-550, 92}, {-550, 35}, {-778, 35}, {-778, -188}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-550, 92}, {-468, 92}, {-468, 25}, {-411, 25}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-550, 92}, {-328, 92}, {-328, 24}, {-294, 24}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-550, 92}, {-206, 92}, {-206, 24}, {-178, 24}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-550, 92}, {-62, 92}, {-62, 56}, {-59, 56}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-550, 92}, {-108, 92}, {-108, -154}, {-151, -154}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, 35}, {112, 35}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {12, -6}, {12, 13}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, -23}, {74, -23}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-53, -82}, {-53, -43}}, color = {245, 121, 0}, thickness = 1));
      connect(externalReact.step, core1R.realIn) annotation(
          Line(points = {{-660, 106}, {-660, -50}, {-746, -50}, {-746, -188}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{74, -212}, {74, -78}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-800, 160}, {220, -300}})));
      end R1MSRE100pcm;
      
      model R1MSRE10pcm
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
        parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 1966.5;
        parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 4.76;
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
        parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
        parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 4.76;
        parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
        parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
        SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 324000, hAsNom = 153000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0.01, Tinf = 500, TpIn_0 = 657.2700, TpOut_0 = 632, TsIn_0 = 546.4, TsOut_0 = 579.3) annotation(
          Placement(transformation(origin = {-36.2, 34.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
        SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
          Placement(transformation(origin = {115.462, -50.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0.08, Tinf = 644.4, T_0 = 579.3) annotation(
          Placement(transformation(origin = {-30.8, -57.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0.08, Tinf = 644.4, T_0 = 546.4) annotation(
          Placement(transformation(origin = {89.2, 21.7333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = 0.075708, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0.08, Tinf = 644.5, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-261, 36.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.066953, volFracNode = 0.01, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-145.6, 5.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 0.6564, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0.08, Tinf = 644.4, T_0 = 632) annotation(
          Placement(transformation(origin = {-182.4, -135.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.066953, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.2776, EnableRad = false, Ar = 3.3039, e = 0.08, Tinf = 644.4, T_0 = 657.2700) annotation(
          Placement(transformation(origin = {-378.4, 5.2}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        MSRE.Components.MSRE1R core1R(numGroups = 6, addRho0 = true, EnableRad = false, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, a_F = -8.71E-05, a_G = -6.66E-05, n_0 = 1, rho_fuel = rhoFuel, rho_grap = rhoGrap, cP_fuel = scpFuel, cP_grap = scpGrap, volDotFuel = volDotFuel, T_inf = 500, kFuel = kFuel, TF1_0 = 644.5, TF2_0 = 657.2700, TG_0 = 675.8561) annotation(
          Placement(transformation(origin = {-495.667, -31}, extent = {{-172.333, -141}, {-78.3333, -47}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
          Placement(transformation(origin = {43, -145}, extent = {{-27, -27}, {27, 27}})));
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
          Placement(transformation(origin = {73, -225}, extent = {{-27, -27}, {27, 27}})));
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 10}) annotation(
          Placement(transformation(origin = {-662.198, 95.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-77, 22}, {-94, 22}, {-94, -34}, {-53, -34}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-17, -34}, {28.5, -34}, {28.5, -50}, {74, -50}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{129, -50}, {140, -50}, {140, 44}, {112, 44}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{76, 44}, {53.5, 44}, {53.5, 22}, {30, 22}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-241, 49}, {-207, 49}, {-207, 37}, {-178, 37}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-126, 37}, {-126, 47}, {-77, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{30, 47}, {156, 47}, {156, -167}, {-151, -167}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-359, 38}, {-325.5, 38}, {-325.5, 49}, {-294, 49}}, color = {204, 0, 0}, thickness = 1));
        connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-715, -188}, {-713, 38}, {-411, 38}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
          Line(points = {{-201, -167}, {-201, -110.5}, {-235, -110.5}, {-235, -110}, {-543, -110}, {-545, -249}, {-778, -249}}, color = {204, 0, 0}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 51}, {-411, 51}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 74}, {-294, 74}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-448, -250}, {-448, 50}, {-178, 50}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-714, -250}, {-714, 114}, {-23.5, 114}, {-23.5, 56}}, color = {220, 138, 221}, thickness = 1));
        connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-714, -250}, {-712, -179}, {-151, -179}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{41, -129}, {-74, -129}, {-74, -25}, {-53, -25}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{41, -129}, {148, -129}, {148, 53}, {112, 53}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
          Line(points = {{-550, 92}, {-550, 35}, {-778, 35}, {-778, -188}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-550, 92}, {-468, 92}, {-468, 25}, {-411, 25}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-550, 92}, {-328, 92}, {-328, 24}, {-294, 24}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-550, 92}, {-206, 92}, {-206, 24}, {-178, 24}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-550, 92}, {-62, 92}, {-62, 56}, {-59, 56}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-550, 92}, {-108, 92}, {-108, -154}, {-151, -154}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, 35}, {112, 35}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {12, -6}, {12, 13}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, -23}, {74, -23}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-53, -82}, {-53, -43}}, color = {245, 121, 0}, thickness = 1));
      connect(externalReact.step, core1R.realIn) annotation(
          Line(points = {{-660, 106}, {-660, -50}, {-746, -50}, {-746, -188}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{74, -212}, {74, -78}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-800, 160}, {220, -300}})));
      end R1MSRE10pcm;
    end R1fullSteps;

    package R9fullSteps
      model R9MSRE2dol
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
        parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2000;
        parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
        parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
        parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
        parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
        parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
        SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 648000, hAsNom = 306000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0, Tinf = 663.53, TpIn_0 = 663.53, TpOut_0 = 632, TsIn_0 = 554.4300, TsOut_0 = 579.3) annotation(
          Placement(transformation(origin = {4.6, 24}, extent = {{-93.6, -46.8}, {140.4, 46.8}})));
        SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
          Placement(transformation(origin = {80.7953, -102.346}, extent = {{-45.5285, -34.1463}, {22.7642, 34.1463}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0, Tinf = 644.4, T_0 = 579.3) annotation(
          Placement(transformation(origin = {-32.1333, -61}, extent = {{-49.2, -16.4}, {32.7999, -65.6}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0, Tinf = 644.4, T_0 = 546.4) annotation(
          Placement(transformation(origin = {164.4, -60.6667}, extent = {{-50, -16.6667}, {33.3333, -66.6667}})));
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = volDotFuel*1, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-259, 34.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-143.467, -2.53333}, extent = {{-59.2, 19.7333}, {39.4667, 78.9333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 8.67*volDotFuel, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0, Tinf = 632, T_0 = 632) annotation(
          Placement(transformation(origin = {-162, -104.8}, extent = {{-52.8, 17.6}, {35.2, 70.4}}, rotation = 180)));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-374.4, 1.86667}, extent = {{-55.6, 18.5333}, {37.0667, 74.1333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
          Placement(transformation(origin = {-68, -180}, extent = {{-17, -17}, {17, 17}})));
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 1000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-528, 122}, extent = {{-24, -32}, {24, 16}})));
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-33, -39.333}, extent = {{-17, -22.6667}, {17, 11.3333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
          Placement(transformation(origin = {36, -200}, extent = {{-17, -17}, {17, 17}})));
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 829.172}) annotation(
          Placement(transformation(origin = {-496.198, -16.1983}, extent = {{-17.7983, 17.7983}, {17.7983, -17.7983}}, rotation = -0)));
        MSRE.Components.MSRE9R msre9r(numGroups = 6, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, n_0 = 1, aF = -8.71E-05, aG = -6.66E-05, addRho0 = true, volF1 = {6.4392e-03, 2.1834e-02, 1.1940e-02, 1.4925e-02, 3.6929e-02, 2.0170e-02, 2.5245e-02, 1.0149e-01, 6.8869e-02}, volF2 = {6.7377e-03, 1.4883e-02, 1.1940e-02, 2.9083e-02, 2.5245e-02, 2.0170e-02, 4.9125e-02, 5.9019e-02, 1.1556e-01}, volG = {0.038208, 0.115387, 0.087659, 0.112221, 0.195186, 0.148354, 0.189837, 0.524644, 0.514219}, volUP = volDotFuel*2, hA = {705.60, 2167.20, 1620.00, 2113.20, 3558.60, 2745.00, 3573.00, 9801.00, 9648.00}, kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, TF1_0 = {641.0644, 636.1926, 662.2593, 692.4323, 634.6308, 650.7339, 669.3933, 633.4354, 638.9545}, TF2_0 = {651.9893, 647.9801, 677.5167, 699.8892, 641.9048, 660.1706, 673.9817, 635.8835, 640.4922}, TG_0 = {664.0462, 653.7083, 692.6790, 714.5485, 648.7869, 674.7077, 686.8285, 640.8657, 647.7290}, Tmix_0 = 663.53, IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, flowFracRegions = {0.061409895150439, 0.138549617100260, 0.234231239311780, 0.565809248437521}, rho_fuel = rhoFuel, cP_fuel = scpFuel, kFuel = kFuel, volDotFuel = volDotFuel, rho_grap = rhoGrap, cP_grap = scpGrap, regionTripTime = {200000, 200000, 200000, 200000}, regionCoastDownK = 0.02, freeConvFF = 0, EnableRad = false, T_inf = 663.53, LF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, LF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, Ac = {7.6800e-03, 4.4928e-02, 1.0368e-01, 1.7222e-01}, ArF1 = {0.5301, 0.5692, 0.3113, 0.3891, 0.9719, 0.5308, 0.6644, 2.0856, 1.4152}, ArF2 = {0.5546, 0.3880, 0.3113, 0.7582, 0.6644, 0.5308, 1.2929, 1.2128, 2.3747}, e = 0) annotation(
          Placement(transformation(origin = {-364.333, -15}, extent = {{-179.667, -147}, {-81.6667, -49}})));
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-71, 1}, {-94, 1}, {-94, -102}, {-68, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-13, -102}, {35, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{81, -102}, {128, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-239, 47}, {-186, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-120, 47}, {-71, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-353, 48}, {-328.5, 48}, {-328.5, 47}, {-292, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-528, 90}, {-468, 90}, {-468, 34}, {-415, 34}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-528, 90}, {-328, 90}, {-328, 22}, {-292, 22}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-528, 90}, {-206, 90}, {-206, 32}, {-186, 32}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-528, 90}, {-299.5, 90}, {-299.5, 63}, {-38, 63}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-528, 90}, {-108, 90}, {-108, -136}, {-124, -136}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-293, -146}, {-293, -162}, {-124, -162}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-434, -146}, {-434, 62}, {-415, 62}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-462, -146}, {-318, -146}, {-318, 72}, {-292, 72}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-200, -146}, {-200, 62}, {-186, 62}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-462, -82}, {-460, -82}, {-460, 48}, {-415, 48}}, color = {204, 0, 0}, thickness = 1));
        connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
          Line(points = {{-528, 90}, {-528, -82}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-462, -146}, {-434, -146}, {-434, 118}, {28, 118}, {28, 63}}, color = {220, 138, 221}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
          Line(points = {{-183, -149}, {-183, -178}, {-528, -178}, {-528, -146}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_in_sFluid, pipeUHXtoHX.PiTempOut) annotation(
          Line(points = {{128, 0}, {200, 0}, {200, -102}, {184, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{-33, -62}, {-68, -62}, {-68, -90}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{-33, -62}, {-33, -80}, {35, -80}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{-33, -62}, {-33, -90}, {128, -90}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{-33, -62}, {-33, -16}, {94, -16}}, color = {245, 121, 0}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-70, -170}, {-70, -168}, {-68, -168}, {-68, -114}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-70, -170}, {-70, -169}, {128, -169}, {128, -115}}, color = {220, 138, 221}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{128, 48}, {224, 48}, {224, -149}, {-124, -149}}, color = {204, 0, 0}, thickness = 1));
      connect(externalReact.step, msre9r.realIn) annotation(
          Line(points = {{-496, -22}, {-494, -22}, {-494, -82}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{36, -192}, {36, -126}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-580, 160}, {240, -240}})));
      end R9MSRE2dol;
      
      model R9MSRE1dol
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
        parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2000;
        parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
        parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
        parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
        parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
        parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
        SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 648000, hAsNom = 306000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0, Tinf = 663.53, TpIn_0 = 663.53, TpOut_0 = 632, TsIn_0 = 554.4300, TsOut_0 = 579.3) annotation(
          Placement(transformation(origin = {4.6, 24}, extent = {{-93.6, -46.8}, {140.4, 46.8}})));
        SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
          Placement(transformation(origin = {80.7953, -102.346}, extent = {{-45.5285, -34.1463}, {22.7642, 34.1463}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0, Tinf = 644.4, T_0 = 579.3) annotation(
          Placement(transformation(origin = {-32.1333, -61}, extent = {{-49.2, -16.4}, {32.7999, -65.6}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0, Tinf = 644.4, T_0 = 546.4) annotation(
          Placement(transformation(origin = {164.4, -60.6667}, extent = {{-50, -16.6667}, {33.3333, -66.6667}})));
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = volDotFuel*1, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-259, 34.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-143.467, -2.53333}, extent = {{-59.2, 19.7333}, {39.4667, 78.9333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 8.67*volDotFuel, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0, Tinf = 632, T_0 = 632) annotation(
          Placement(transformation(origin = {-162, -104.8}, extent = {{-52.8, 17.6}, {35.2, 70.4}}, rotation = 180)));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-374.4, 1.86667}, extent = {{-55.6, 18.5333}, {37.0667, 74.1333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
          Placement(transformation(origin = {-68, -180}, extent = {{-17, -17}, {17, 17}})));
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 1000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-528, 122}, extent = {{-24, -32}, {24, 16}})));
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-33, -39.333}, extent = {{-17, -22.6667}, {17, 11.3333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
          Placement(transformation(origin = {36, -200}, extent = {{-17, -17}, {17, 17}})));
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 414.586}) annotation(
          Placement(transformation(origin = {-496.198, -16.1983}, extent = {{-17.7983, 17.7983}, {17.7983, -17.7983}}, rotation = -0)));
        MSRE.Components.MSRE9R msre9r(numGroups = 6, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, n_0 = 1, aF = -8.71E-05, aG = -6.66E-05, addRho0 = true, volF1 = {6.4392e-03, 2.1834e-02, 1.1940e-02, 1.4925e-02, 3.6929e-02, 2.0170e-02, 2.5245e-02, 1.0149e-01, 6.8869e-02}, volF2 = {6.7377e-03, 1.4883e-02, 1.1940e-02, 2.9083e-02, 2.5245e-02, 2.0170e-02, 4.9125e-02, 5.9019e-02, 1.1556e-01}, volG = {0.038208, 0.115387, 0.087659, 0.112221, 0.195186, 0.148354, 0.189837, 0.524644, 0.514219}, volUP = volDotFuel*2, hA = {705.60, 2167.20, 1620.00, 2113.20, 3558.60, 2745.00, 3573.00, 9801.00, 9648.00}, kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, TF1_0 = {641.0644, 636.1926, 662.2593, 692.4323, 634.6308, 650.7339, 669.3933, 633.4354, 638.9545}, TF2_0 = {651.9893, 647.9801, 677.5167, 699.8892, 641.9048, 660.1706, 673.9817, 635.8835, 640.4922}, TG_0 = {664.0462, 653.7083, 692.6790, 714.5485, 648.7869, 674.7077, 686.8285, 640.8657, 647.7290}, Tmix_0 = 663.53, IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, flowFracRegions = {0.061409895150439, 0.138549617100260, 0.234231239311780, 0.565809248437521}, rho_fuel = rhoFuel, cP_fuel = scpFuel, kFuel = kFuel, volDotFuel = volDotFuel, rho_grap = rhoGrap, cP_grap = scpGrap, regionTripTime = {200000, 200000, 200000, 200000}, regionCoastDownK = 0.02, freeConvFF = 0, EnableRad = false, T_inf = 663.53, LF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, LF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, Ac = {7.6800e-03, 4.4928e-02, 1.0368e-01, 1.7222e-01}, ArF1 = {0.5301, 0.5692, 0.3113, 0.3891, 0.9719, 0.5308, 0.6644, 2.0856, 1.4152}, ArF2 = {0.5546, 0.3880, 0.3113, 0.7582, 0.6644, 0.5308, 1.2929, 1.2128, 2.3747}, e = 0) annotation(
          Placement(transformation(origin = {-364.333, -15}, extent = {{-179.667, -147}, {-81.6667, -49}})));
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-71, 1}, {-94, 1}, {-94, -102}, {-68, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-13, -102}, {35, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{81, -102}, {128, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-239, 47}, {-186, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-120, 47}, {-71, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-353, 48}, {-328.5, 48}, {-328.5, 47}, {-292, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-528, 90}, {-468, 90}, {-468, 34}, {-415, 34}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-528, 90}, {-328, 90}, {-328, 22}, {-292, 22}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-528, 90}, {-206, 90}, {-206, 32}, {-186, 32}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-528, 90}, {-299.5, 90}, {-299.5, 63}, {-38, 63}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-528, 90}, {-108, 90}, {-108, -136}, {-124, -136}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-293, -146}, {-293, -162}, {-124, -162}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-434, -146}, {-434, 62}, {-415, 62}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-462, -146}, {-318, -146}, {-318, 72}, {-292, 72}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-200, -146}, {-200, 62}, {-186, 62}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-462, -82}, {-460, -82}, {-460, 48}, {-415, 48}}, color = {204, 0, 0}, thickness = 1));
        connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
          Line(points = {{-528, 90}, {-528, -82}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-462, -146}, {-434, -146}, {-434, 118}, {28, 118}, {28, 63}}, color = {220, 138, 221}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
          Line(points = {{-183, -149}, {-183, -178}, {-528, -178}, {-528, -146}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_in_sFluid, pipeUHXtoHX.PiTempOut) annotation(
          Line(points = {{128, 0}, {200, 0}, {200, -102}, {184, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{-33, -62}, {-68, -62}, {-68, -90}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{-33, -62}, {-33, -80}, {35, -80}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{-33, -62}, {-33, -90}, {128, -90}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{-33, -62}, {-33, -16}, {94, -16}}, color = {245, 121, 0}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-70, -170}, {-70, -168}, {-68, -168}, {-68, -114}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-70, -170}, {-70, -169}, {128, -169}, {128, -115}}, color = {220, 138, 221}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{128, 48}, {224, 48}, {224, -149}, {-124, -149}}, color = {204, 0, 0}, thickness = 1));
      connect(externalReact.step, msre9r.realIn) annotation(
          Line(points = {{-496, -22}, {-494, -22}, {-494, -82}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{36, -192}, {36, -126}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-580, 160}, {240, -240}})));
      end R9MSRE1dol;
      
      model R9MSREhalfDol
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
        parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2000;
        parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
        parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
        parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
        parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
        parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
        SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 648000, hAsNom = 306000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0, Tinf = 663.53, TpIn_0 = 663.53, TpOut_0 = 632, TsIn_0 = 554.4300, TsOut_0 = 579.3) annotation(
          Placement(transformation(origin = {4.6, 24}, extent = {{-93.6, -46.8}, {140.4, 46.8}})));
        SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
          Placement(transformation(origin = {80.7953, -102.346}, extent = {{-45.5285, -34.1463}, {22.7642, 34.1463}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0, Tinf = 644.4, T_0 = 579.3) annotation(
          Placement(transformation(origin = {-32.1333, -61}, extent = {{-49.2, -16.4}, {32.7999, -65.6}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0, Tinf = 644.4, T_0 = 546.4) annotation(
          Placement(transformation(origin = {164.4, -60.6667}, extent = {{-50, -16.6667}, {33.3333, -66.6667}})));
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = volDotFuel*1, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-259, 34.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-143.467, -2.53333}, extent = {{-59.2, 19.7333}, {39.4667, 78.9333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 8.67*volDotFuel, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0, Tinf = 632, T_0 = 632) annotation(
          Placement(transformation(origin = {-162, -104.8}, extent = {{-52.8, 17.6}, {35.2, 70.4}}, rotation = 180)));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-374.4, 1.86667}, extent = {{-55.6, 18.5333}, {37.0667, 74.1333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
          Placement(transformation(origin = {-68, -180}, extent = {{-17, -17}, {17, 17}})));
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 1000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-528, 122}, extent = {{-24, -32}, {24, 16}})));
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-33, -39.333}, extent = {{-17, -22.6667}, {17, 11.3333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
          Placement(transformation(origin = {36, -200}, extent = {{-17, -17}, {17, 17}})));
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 207.293}) annotation(
          Placement(transformation(origin = {-496.198, -16.1983}, extent = {{-17.7983, 17.7983}, {17.7983, -17.7983}}, rotation = -0)));
        MSRE.Components.MSRE9R msre9r(numGroups = 6, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, n_0 = 1, aF = -8.71E-05, aG = -6.66E-05, addRho0 = true, volF1 = {6.4392e-03, 2.1834e-02, 1.1940e-02, 1.4925e-02, 3.6929e-02, 2.0170e-02, 2.5245e-02, 1.0149e-01, 6.8869e-02}, volF2 = {6.7377e-03, 1.4883e-02, 1.1940e-02, 2.9083e-02, 2.5245e-02, 2.0170e-02, 4.9125e-02, 5.9019e-02, 1.1556e-01}, volG = {0.038208, 0.115387, 0.087659, 0.112221, 0.195186, 0.148354, 0.189837, 0.524644, 0.514219}, volUP = volDotFuel*2, hA = {705.60, 2167.20, 1620.00, 2113.20, 3558.60, 2745.00, 3573.00, 9801.00, 9648.00}, kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, TF1_0 = {641.0644, 636.1926, 662.2593, 692.4323, 634.6308, 650.7339, 669.3933, 633.4354, 638.9545}, TF2_0 = {651.9893, 647.9801, 677.5167, 699.8892, 641.9048, 660.1706, 673.9817, 635.8835, 640.4922}, TG_0 = {664.0462, 653.7083, 692.6790, 714.5485, 648.7869, 674.7077, 686.8285, 640.8657, 647.7290}, Tmix_0 = 663.53, IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, flowFracRegions = {0.061409895150439, 0.138549617100260, 0.234231239311780, 0.565809248437521}, rho_fuel = rhoFuel, cP_fuel = scpFuel, kFuel = kFuel, volDotFuel = volDotFuel, rho_grap = rhoGrap, cP_grap = scpGrap, regionTripTime = {200000, 200000, 200000, 200000}, regionCoastDownK = 0.02, freeConvFF = 0, EnableRad = false, T_inf = 663.53, LF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, LF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, Ac = {7.6800e-03, 4.4928e-02, 1.0368e-01, 1.7222e-01}, ArF1 = {0.5301, 0.5692, 0.3113, 0.3891, 0.9719, 0.5308, 0.6644, 2.0856, 1.4152}, ArF2 = {0.5546, 0.3880, 0.3113, 0.7582, 0.6644, 0.5308, 1.2929, 1.2128, 2.3747}, e = 0) annotation(
          Placement(transformation(origin = {-364.333, -15}, extent = {{-179.667, -147}, {-81.6667, -49}})));
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-71, 1}, {-94, 1}, {-94, -102}, {-68, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-13, -102}, {35, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{81, -102}, {128, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-239, 47}, {-186, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-120, 47}, {-71, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-353, 48}, {-328.5, 48}, {-328.5, 47}, {-292, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-528, 90}, {-468, 90}, {-468, 34}, {-415, 34}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-528, 90}, {-328, 90}, {-328, 22}, {-292, 22}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-528, 90}, {-206, 90}, {-206, 32}, {-186, 32}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-528, 90}, {-299.5, 90}, {-299.5, 63}, {-38, 63}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-528, 90}, {-108, 90}, {-108, -136}, {-124, -136}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-293, -146}, {-293, -162}, {-124, -162}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-434, -146}, {-434, 62}, {-415, 62}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-462, -146}, {-318, -146}, {-318, 72}, {-292, 72}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-200, -146}, {-200, 62}, {-186, 62}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-462, -82}, {-460, -82}, {-460, 48}, {-415, 48}}, color = {204, 0, 0}, thickness = 1));
        connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
          Line(points = {{-528, 90}, {-528, -82}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-462, -146}, {-434, -146}, {-434, 118}, {28, 118}, {28, 63}}, color = {220, 138, 221}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
          Line(points = {{-183, -149}, {-183, -178}, {-528, -178}, {-528, -146}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_in_sFluid, pipeUHXtoHX.PiTempOut) annotation(
          Line(points = {{128, 0}, {200, 0}, {200, -102}, {184, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{-33, -62}, {-68, -62}, {-68, -90}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{-33, -62}, {-33, -80}, {35, -80}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{-33, -62}, {-33, -90}, {128, -90}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{-33, -62}, {-33, -16}, {94, -16}}, color = {245, 121, 0}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-70, -170}, {-70, -168}, {-68, -168}, {-68, -114}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-70, -170}, {-70, -169}, {128, -169}, {128, -115}}, color = {220, 138, 221}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{128, 48}, {224, 48}, {224, -149}, {-124, -149}}, color = {204, 0, 0}, thickness = 1));
      connect(externalReact.step, msre9r.realIn) annotation(
          Line(points = {{-496, -22}, {-494, -22}, {-494, -82}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{36, -192}, {36, -126}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-580, 160}, {240, -240}})));
      end R9MSREhalfDol;
      
      model R9MSREpOneDol
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
        parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2000;
        parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
        parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
        parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
        parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
        parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
        SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 648000, hAsNom = 306000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0, Tinf = 663.53, TpIn_0 = 663.53, TpOut_0 = 632, TsIn_0 = 554.4300, TsOut_0 = 579.3) annotation(
          Placement(transformation(origin = {4.6, 24}, extent = {{-93.6, -46.8}, {140.4, 46.8}})));
        SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
          Placement(transformation(origin = {80.7953, -102.346}, extent = {{-45.5285, -34.1463}, {22.7642, 34.1463}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0, Tinf = 644.4, T_0 = 579.3) annotation(
          Placement(transformation(origin = {-32.1333, -61}, extent = {{-49.2, -16.4}, {32.7999, -65.6}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0, Tinf = 644.4, T_0 = 546.4) annotation(
          Placement(transformation(origin = {164.4, -60.6667}, extent = {{-50, -16.6667}, {33.3333, -66.6667}})));
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = volDotFuel*1, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-259, 34.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-143.467, -2.53333}, extent = {{-59.2, 19.7333}, {39.4667, 78.9333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 8.67*volDotFuel, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0, Tinf = 632, T_0 = 632) annotation(
          Placement(transformation(origin = {-162, -104.8}, extent = {{-52.8, 17.6}, {35.2, 70.4}}, rotation = 180)));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-374.4, 1.86667}, extent = {{-55.6, 18.5333}, {37.0667, 74.1333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
          Placement(transformation(origin = {-68, -180}, extent = {{-17, -17}, {17, 17}})));
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 1000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-528, 122}, extent = {{-24, -32}, {24, 16}})));
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-33, -39.333}, extent = {{-17, -22.6667}, {17, 11.3333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
          Placement(transformation(origin = {36, -200}, extent = {{-17, -17}, {17, 17}})));
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 41.459}) annotation(
          Placement(transformation(origin = {-496.198, -16.1983}, extent = {{-17.7983, 17.7983}, {17.7983, -17.7983}}, rotation = -0)));
        MSRE.Components.MSRE9R msre9r(numGroups = 6, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, n_0 = 1, aF = -8.71E-05, aG = -6.66E-05, addRho0 = true, volF1 = {6.4392e-03, 2.1834e-02, 1.1940e-02, 1.4925e-02, 3.6929e-02, 2.0170e-02, 2.5245e-02, 1.0149e-01, 6.8869e-02}, volF2 = {6.7377e-03, 1.4883e-02, 1.1940e-02, 2.9083e-02, 2.5245e-02, 2.0170e-02, 4.9125e-02, 5.9019e-02, 1.1556e-01}, volG = {0.038208, 0.115387, 0.087659, 0.112221, 0.195186, 0.148354, 0.189837, 0.524644, 0.514219}, volUP = volDotFuel*2, hA = {705.60, 2167.20, 1620.00, 2113.20, 3558.60, 2745.00, 3573.00, 9801.00, 9648.00}, kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, TF1_0 = {641.0644, 636.1926, 662.2593, 692.4323, 634.6308, 650.7339, 669.3933, 633.4354, 638.9545}, TF2_0 = {651.9893, 647.9801, 677.5167, 699.8892, 641.9048, 660.1706, 673.9817, 635.8835, 640.4922}, TG_0 = {664.0462, 653.7083, 692.6790, 714.5485, 648.7869, 674.7077, 686.8285, 640.8657, 647.7290}, Tmix_0 = 663.53, IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, flowFracRegions = {0.061409895150439, 0.138549617100260, 0.234231239311780, 0.565809248437521}, rho_fuel = rhoFuel, cP_fuel = scpFuel, kFuel = kFuel, volDotFuel = volDotFuel, rho_grap = rhoGrap, cP_grap = scpGrap, regionTripTime = {200000, 200000, 200000, 200000}, regionCoastDownK = 0.02, freeConvFF = 0, EnableRad = false, T_inf = 663.53, LF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, LF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, Ac = {7.6800e-03, 4.4928e-02, 1.0368e-01, 1.7222e-01}, ArF1 = {0.5301, 0.5692, 0.3113, 0.3891, 0.9719, 0.5308, 0.6644, 2.0856, 1.4152}, ArF2 = {0.5546, 0.3880, 0.3113, 0.7582, 0.6644, 0.5308, 1.2929, 1.2128, 2.3747}, e = 0) annotation(
          Placement(transformation(origin = {-364.333, -15}, extent = {{-179.667, -147}, {-81.6667, -49}})));
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-71, 1}, {-94, 1}, {-94, -102}, {-68, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-13, -102}, {35, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{81, -102}, {128, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-239, 47}, {-186, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-120, 47}, {-71, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-353, 48}, {-328.5, 48}, {-328.5, 47}, {-292, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-528, 90}, {-468, 90}, {-468, 34}, {-415, 34}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-528, 90}, {-328, 90}, {-328, 22}, {-292, 22}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-528, 90}, {-206, 90}, {-206, 32}, {-186, 32}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-528, 90}, {-299.5, 90}, {-299.5, 63}, {-38, 63}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-528, 90}, {-108, 90}, {-108, -136}, {-124, -136}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-293, -146}, {-293, -162}, {-124, -162}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-434, -146}, {-434, 62}, {-415, 62}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-462, -146}, {-318, -146}, {-318, 72}, {-292, 72}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-200, -146}, {-200, 62}, {-186, 62}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-462, -82}, {-460, -82}, {-460, 48}, {-415, 48}}, color = {204, 0, 0}, thickness = 1));
        connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
          Line(points = {{-528, 90}, {-528, -82}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-462, -146}, {-434, -146}, {-434, 118}, {28, 118}, {28, 63}}, color = {220, 138, 221}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
          Line(points = {{-183, -149}, {-183, -178}, {-528, -178}, {-528, -146}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_in_sFluid, pipeUHXtoHX.PiTempOut) annotation(
          Line(points = {{128, 0}, {200, 0}, {200, -102}, {184, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{-33, -62}, {-68, -62}, {-68, -90}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{-33, -62}, {-33, -80}, {35, -80}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{-33, -62}, {-33, -90}, {128, -90}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{-33, -62}, {-33, -16}, {94, -16}}, color = {245, 121, 0}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-70, -170}, {-70, -168}, {-68, -168}, {-68, -114}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-70, -170}, {-70, -169}, {128, -169}, {128, -115}}, color = {220, 138, 221}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{128, 48}, {224, 48}, {224, -149}, {-124, -149}}, color = {204, 0, 0}, thickness = 1));
      connect(externalReact.step, msre9r.realIn) annotation(
          Line(points = {{-496, -22}, {-494, -22}, {-494, -82}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{36, -192}, {36, -126}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-580, 160}, {240, -240}})));
      end R9MSREpOneDol;
      
      model R9MSRE100pcm
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
        parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2000;
        parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
        parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
        parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
        parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
        parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
        SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 648000, hAsNom = 306000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0, Tinf = 663.53, TpIn_0 = 663.53, TpOut_0 = 632, TsIn_0 = 554.4300, TsOut_0 = 579.3) annotation(
          Placement(transformation(origin = {4.6, 24}, extent = {{-93.6, -46.8}, {140.4, 46.8}})));
        SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
          Placement(transformation(origin = {80.7953, -102.346}, extent = {{-45.5285, -34.1463}, {22.7642, 34.1463}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0, Tinf = 644.4, T_0 = 579.3) annotation(
          Placement(transformation(origin = {-32.1333, -61}, extent = {{-49.2, -16.4}, {32.7999, -65.6}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0, Tinf = 644.4, T_0 = 546.4) annotation(
          Placement(transformation(origin = {164.4, -60.6667}, extent = {{-50, -16.6667}, {33.3333, -66.6667}})));
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = volDotFuel*1, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-259, 34.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-143.467, -2.53333}, extent = {{-59.2, 19.7333}, {39.4667, 78.9333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 8.67*volDotFuel, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0, Tinf = 632, T_0 = 632) annotation(
          Placement(transformation(origin = {-162, -104.8}, extent = {{-52.8, 17.6}, {35.2, 70.4}}, rotation = 180)));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-374.4, 1.86667}, extent = {{-55.6, 18.5333}, {37.0667, 74.1333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
          Placement(transformation(origin = {-68, -180}, extent = {{-17, -17}, {17, 17}})));
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 1000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-528, 122}, extent = {{-24, -32}, {24, 16}})));
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-33, -39.333}, extent = {{-17, -22.6667}, {17, 11.3333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
          Placement(transformation(origin = {36, -200}, extent = {{-17, -17}, {17, 17}})));
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 100}) annotation(
          Placement(transformation(origin = {-496.198, -16.1983}, extent = {{-17.7983, 17.7983}, {17.7983, -17.7983}}, rotation = -0)));
        MSRE.Components.MSRE9R msre9r(numGroups = 6, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, n_0 = 1, aF = -8.71E-05, aG = -6.66E-05, addRho0 = true, volF1 = {6.4392e-03, 2.1834e-02, 1.1940e-02, 1.4925e-02, 3.6929e-02, 2.0170e-02, 2.5245e-02, 1.0149e-01, 6.8869e-02}, volF2 = {6.7377e-03, 1.4883e-02, 1.1940e-02, 2.9083e-02, 2.5245e-02, 2.0170e-02, 4.9125e-02, 5.9019e-02, 1.1556e-01}, volG = {0.038208, 0.115387, 0.087659, 0.112221, 0.195186, 0.148354, 0.189837, 0.524644, 0.514219}, volUP = volDotFuel*2, hA = {705.60, 2167.20, 1620.00, 2113.20, 3558.60, 2745.00, 3573.00, 9801.00, 9648.00}, kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, TF1_0 = {641.0644, 636.1926, 662.2593, 692.4323, 634.6308, 650.7339, 669.3933, 633.4354, 638.9545}, TF2_0 = {651.9893, 647.9801, 677.5167, 699.8892, 641.9048, 660.1706, 673.9817, 635.8835, 640.4922}, TG_0 = {664.0462, 653.7083, 692.6790, 714.5485, 648.7869, 674.7077, 686.8285, 640.8657, 647.7290}, Tmix_0 = 663.53, IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, flowFracRegions = {0.061409895150439, 0.138549617100260, 0.234231239311780, 0.565809248437521}, rho_fuel = rhoFuel, cP_fuel = scpFuel, kFuel = kFuel, volDotFuel = volDotFuel, rho_grap = rhoGrap, cP_grap = scpGrap, regionTripTime = {200000, 200000, 200000, 200000}, regionCoastDownK = 0.02, freeConvFF = 0, EnableRad = false, T_inf = 663.53, LF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, LF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, Ac = {7.6800e-03, 4.4928e-02, 1.0368e-01, 1.7222e-01}, ArF1 = {0.5301, 0.5692, 0.3113, 0.3891, 0.9719, 0.5308, 0.6644, 2.0856, 1.4152}, ArF2 = {0.5546, 0.3880, 0.3113, 0.7582, 0.6644, 0.5308, 1.2929, 1.2128, 2.3747}, e = 0) annotation(
          Placement(transformation(origin = {-364.333, -15}, extent = {{-179.667, -147}, {-81.6667, -49}})));
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-71, 1}, {-94, 1}, {-94, -102}, {-68, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-13, -102}, {35, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{81, -102}, {128, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-239, 47}, {-186, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-120, 47}, {-71, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-353, 48}, {-328.5, 48}, {-328.5, 47}, {-292, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-528, 90}, {-468, 90}, {-468, 34}, {-415, 34}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-528, 90}, {-328, 90}, {-328, 22}, {-292, 22}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-528, 90}, {-206, 90}, {-206, 32}, {-186, 32}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-528, 90}, {-299.5, 90}, {-299.5, 63}, {-38, 63}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-528, 90}, {-108, 90}, {-108, -136}, {-124, -136}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-293, -146}, {-293, -162}, {-124, -162}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-434, -146}, {-434, 62}, {-415, 62}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-462, -146}, {-318, -146}, {-318, 72}, {-292, 72}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-200, -146}, {-200, 62}, {-186, 62}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-462, -82}, {-460, -82}, {-460, 48}, {-415, 48}}, color = {204, 0, 0}, thickness = 1));
        connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
          Line(points = {{-528, 90}, {-528, -82}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-462, -146}, {-434, -146}, {-434, 118}, {28, 118}, {28, 63}}, color = {220, 138, 221}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
          Line(points = {{-183, -149}, {-183, -178}, {-528, -178}, {-528, -146}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_in_sFluid, pipeUHXtoHX.PiTempOut) annotation(
          Line(points = {{128, 0}, {200, 0}, {200, -102}, {184, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{-33, -62}, {-68, -62}, {-68, -90}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{-33, -62}, {-33, -80}, {35, -80}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{-33, -62}, {-33, -90}, {128, -90}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{-33, -62}, {-33, -16}, {94, -16}}, color = {245, 121, 0}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-70, -170}, {-70, -168}, {-68, -168}, {-68, -114}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-70, -170}, {-70, -169}, {128, -169}, {128, -115}}, color = {220, 138, 221}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{128, 48}, {224, 48}, {224, -149}, {-124, -149}}, color = {204, 0, 0}, thickness = 1));
      connect(externalReact.step, msre9r.realIn) annotation(
          Line(points = {{-496, -22}, {-494, -22}, {-494, -82}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{36, -192}, {36, -126}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-580, 160}, {240, -240}})));
      end R9MSRE100pcm;
      
      model R9MSRE10pcm
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.075708;
        parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2146.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2000;
        parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
        parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.053627;
        parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1922;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
        parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
        parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1860;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
        parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 577.80;
        SMD_MSR_Modelica.HeatTransport.HeatExchanger heatExchanger(vol_P = 0.1623, vol_T = 0.023059, vol_S = 0.051863, rhoP = rhoFuel, rhoT = rhoHXtube, rhoS = rhoCoolant, cP_P = scpFuel, cP_T = scpHXtube, cP_S = scpCoolant, VdotPnom = volDotFuel, VdotSnom = volDotCoolant, hApNom = 648000, hAsNom = 306000, EnableRad = false, AcShell = 0.1298, AcTube = 0.020276, ArShell = 3.1165, Kp = kFuel, Ks = kCoolant, L_shell = 2.4400, L_tube = 3.7677, e = 0, Tinf = 663.53, TpIn_0 = 663.53, TpOut_0 = 632, TsIn_0 = 554.4300, TsOut_0 = 579.3) annotation(
          Placement(transformation(origin = {4.6, 24}, extent = {{-93.6, -46.8}, {140.4, 46.8}})));
        SMD_MSR_Modelica.HeatTransport.UHX uhx(Tp_0 = 546.4, vDot = volDotCoolant, cP = scpCoolant, rho = rhoCoolant, vol = 0.204183209633634) annotation(
          Placement(transformation(origin = {80.7953, -102.346}, extent = {{-45.5285, -34.1463}, {22.7642, 34.1463}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(vol = 0.2526, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 19.932, EnableRad = false, Ar = 7.9559, e = 0, Tinf = 644.4, T_0 = 579.3) annotation(
          Placement(transformation(origin = {-32.1333, -61}, extent = {{-49.2, -16.4}, {32.7999, -65.6}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(vol = 0.4419, volFracNode = 0.1, vDotNom = volDotCoolant, rho = rhoCoolant, cP = scpCoolant, K = kCoolant, Ac = 0.012673, L = 34.870, EnableRad = false, Ar = 13.918, e = 0, Tinf = 644.4, T_0 = 546.4) annotation(
          Placement(transformation(origin = {164.4, -60.6667}, extent = {{-50, -16.6667}, {33.3333, -66.6667}})));
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(vol = volDotFuel*1, rho = rhoFuel, cP = scpFuel, vDotNom = volDotFuel, DHRS_tK = 10, DHRS_MaxP_Rm(displayUnit = "MW") = 320000, DHRS_P_Bleed = 0, DHRS_time = 1000000, K = kFuel, Ac = 0.012673, L = 5.9741, Ar = 0.075708, EnableRad = false, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-259, 34.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-143.467, -2.53333}, extent = {{-59.2, 19.7333}, {39.4667, 78.9333}})));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(vol = 8.67*volDotFuel, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 51.796, EnableRad = false, Ar = 20.674, e = 0, Tinf = 632, T_0 = 632) annotation(
          Placement(transformation(origin = {-162, -104.8}, extent = {{-52.8, 17.6}, {35.2, 70.4}}, rotation = 180)));
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(vol = 0.104780118, volFracNode = 0.1, vDotNom = volDotFuel, rho = rhoFuel, cP = scpFuel, K = kFuel, Ac = 0.012673, L = 8.268120397, EnableRad = false, Ar = 3.300161199, e = 0, Tinf = 663.53, T_0 = 663.53) annotation(
          Placement(transformation(origin = {-374.4, 1.86667}, extent = {{-55.6, 18.5333}, {37.0667, 74.1333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0) annotation(
          Placement(transformation(origin = {-68, -180}, extent = {{-17, -17}, {17, 17}})));
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(numRampUp = 1, rampUpK = {1}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 1000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-528, 122}, extent = {{-24, -32}, {24, 16}})));
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(numRampUp = 1, rampUpK = {100}, rampUpTo = {1}, rampUpTime = {0}, tripK = 50, tripTime = 10000000, freeConvFF = 0) annotation(
          Placement(transformation(origin = {-33, -39.333}, extent = {{-17, -22.6667}, {17, 11.3333}})));
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 8E6) annotation(
          Placement(transformation(origin = {36, -200}, extent = {{-17, -17}, {17, 17}})));
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(numSteps = 2, stepTime = {0, 4000}, amplitude = {0, 10}) annotation(
          Placement(transformation(origin = {-496.198, -16.1983}, extent = {{-17.7983, 17.7983}, {17.7983, -17.7983}}, rotation = -0)));
        MSRE.Components.MSRE9R msre9r(numGroups = 6, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, LAMBDA = 2.400E-04, n_0 = 1, aF = -8.71E-05, aG = -6.66E-05, addRho0 = true, volF1 = {6.4392e-03, 2.1834e-02, 1.1940e-02, 1.4925e-02, 3.6929e-02, 2.0170e-02, 2.5245e-02, 1.0149e-01, 6.8869e-02}, volF2 = {6.7377e-03, 1.4883e-02, 1.1940e-02, 2.9083e-02, 2.5245e-02, 2.0170e-02, 4.9125e-02, 5.9019e-02, 1.1556e-01}, volG = {0.038208, 0.115387, 0.087659, 0.112221, 0.195186, 0.148354, 0.189837, 0.524644, 0.514219}, volUP = volDotFuel*2, hA = {705.60, 2167.20, 1620.00, 2113.20, 3558.60, 2745.00, 3573.00, 9801.00, 9648.00}, kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, TF1_0 = {641.0644, 636.1926, 662.2593, 692.4323, 634.6308, 650.7339, 669.3933, 633.4354, 638.9545}, TF2_0 = {651.9893, 647.9801, 677.5167, 699.8892, 641.9048, 660.1706, 673.9817, 635.8835, 640.4922}, TG_0 = {664.0462, 653.7083, 692.6790, 714.5485, 648.7869, 674.7077, 686.8285, 640.8657, 647.7290}, Tmix_0 = 663.53, IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, flowFracRegions = {0.061409895150439, 0.138549617100260, 0.234231239311780, 0.565809248437521}, rho_fuel = rhoFuel, cP_fuel = scpFuel, kFuel = kFuel, volDotFuel = volDotFuel, rho_grap = rhoGrap, cP_grap = scpGrap, regionTripTime = {200000, 200000, 200000, 200000}, regionCoastDownK = 0.02, freeConvFF = 0, EnableRad = false, T_inf = 663.53, LF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, LF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, Ac = {7.6800e-03, 4.4928e-02, 1.0368e-01, 1.7222e-01}, ArF1 = {0.5301, 0.5692, 0.3113, 0.3891, 0.9719, 0.5308, 0.6644, 2.0856, 1.4152}, ArF2 = {0.5546, 0.3880, 0.3113, 0.7582, 0.6644, 0.5308, 1.2929, 1.2128, 2.3747}, e = 0) annotation(
          Placement(transformation(origin = {-364.333, -15}, extent = {{-179.667, -147}, {-81.6667, -49}})));
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-71, 1}, {-94, 1}, {-94, -102}, {-68, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-13, -102}, {35, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{81, -102}, {128, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-239, 47}, {-186, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-120, 47}, {-71, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-353, 48}, {-328.5, 48}, {-328.5, 47}, {-292, 47}}, color = {204, 0, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-528, 90}, {-468, 90}, {-468, 34}, {-415, 34}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-528, 90}, {-328, 90}, {-328, 22}, {-292, 22}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-528, 90}, {-206, 90}, {-206, 32}, {-186, 32}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-528, 90}, {-299.5, 90}, {-299.5, 63}, {-38, 63}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-528, 90}, {-108, 90}, {-108, -136}, {-124, -136}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-293, -146}, {-293, -162}, {-124, -162}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-434, -146}, {-434, 62}, {-415, 62}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-462, -146}, {-318, -146}, {-318, 72}, {-292, 72}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-462, -146}, {-200, -146}, {-200, 62}, {-186, 62}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-462, -82}, {-460, -82}, {-460, 48}, {-415, 48}}, color = {204, 0, 0}, thickness = 1));
        connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
          Line(points = {{-528, 90}, {-528, -82}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-462, -146}, {-434, -146}, {-434, 118}, {28, 118}, {28, 63}}, color = {220, 138, 221}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
          Line(points = {{-183, -149}, {-183, -178}, {-528, -178}, {-528, -146}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_in_sFluid, pipeUHXtoHX.PiTempOut) annotation(
          Line(points = {{128, 0}, {200, 0}, {200, -102}, {184, -102}}, color = {204, 0, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{-33, -62}, {-68, -62}, {-68, -90}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{-33, -62}, {-33, -80}, {35, -80}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{-33, -62}, {-33, -90}, {128, -90}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{-33, -62}, {-33, -16}, {94, -16}}, color = {245, 121, 0}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-70, -170}, {-70, -168}, {-68, -168}, {-68, -114}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-70, -170}, {-70, -169}, {128, -169}, {128, -115}}, color = {220, 138, 221}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{128, 48}, {224, 48}, {224, -149}, {-124, -149}}, color = {204, 0, 0}, thickness = 1));
      connect(externalReact.step, msre9r.realIn) annotation(
          Line(points = {{-496, -22}, {-494, -22}, {-494, -82}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{36, -192}, {36, -126}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-580, 160}, {240, -240}})));
      end R9MSRE10pcm;
    end R9fullSteps;
  end Transients;
end MSRE;