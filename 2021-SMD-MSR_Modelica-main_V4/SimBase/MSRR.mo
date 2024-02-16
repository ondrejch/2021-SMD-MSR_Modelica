package MSRR
   
package Components
    model FuelChannel
      parameter SMD_MSR_Modelica.Units.Volume vol_FN1;
      parameter SMD_MSR_Modelica.Units.Volume vol_FN2;
      parameter SMD_MSR_Modelica.Units.Volume vol_GN;
      parameter SMD_MSR_Modelica.Units.Density rho_fuel;
      parameter SMD_MSR_Modelica.Units.Density rho_grap;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_fuel;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_grap;
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate Vdot_fuelNom;
      parameter SMD_MSR_Modelica.Units.VolumeImportance kFN1;
      parameter SMD_MSR_Modelica.Units.VolumeImportance kFN2;
      parameter SMD_MSR_Modelica.Units.VolumeImportance kG;
      parameter SMD_MSR_Modelica.Units.Convection hAnom;
      parameter SMD_MSR_Modelica.Units.HeatTransferFraction kHT_FN1;
      parameter SMD_MSR_Modelica.Units.HeatTransferFraction kHT_FN2;
      parameter SMD_MSR_Modelica.Units.Temperature TF1_0;
      parameter SMD_MSR_Modelica.Units.Temperature TF2_0;
      parameter SMD_MSR_Modelica.Units.Temperature TG_0;
      parameter SMD_MSR_Modelica.Units.FlowFraction regionFlowFrac;
      parameter SMD_MSR_Modelica.Units.Conductivity KF;
      parameter SMD_MSR_Modelica.Units.Area Ac;
      parameter SMD_MSR_Modelica.Units.Length LF1;
      parameter SMD_MSR_Modelica.Units.Length LF2;
      parameter Boolean OutterRegion;
      parameter SMD_MSR_Modelica.Units.Area ArF1;
      parameter SMD_MSR_Modelica.Units.Area ArF2;
      parameter SMD_MSR_Modelica.Units.Emissivity e;
      parameter SMD_MSR_Modelica.Units.Temperature Tinf;
      SMD_MSR_Modelica.Units.MassFlowRate mdot_fuel;
      SMD_MSR_Modelica.Units.Convection hA;
      SMD_MSR_Modelica.Units.Mass m_FN1;
      SMD_MSR_Modelica.Units.Mass m_FN2;
      SMD_MSR_Modelica.Units.Mass m_GN;
      SMD_MSR_Modelica.Units.Power powFN1;
      SMD_MSR_Modelica.Units.Power powFN2;
      SMD_MSR_Modelica.Units.Power powGN;
      SMD_MSR_Modelica.Units.Power flowPowFN1;
      SMD_MSR_Modelica.Units.Power flowPowFN2;
      SMD_MSR_Modelica.Units.Power condPowFN1;
      SMD_MSR_Modelica.Units.Power condPowFN2;
      SMD_MSR_Modelica.Units.Power convPowFN1;
      SMD_MSR_Modelica.Units.Power convPowFN2;
      SMD_MSR_Modelica.Units.Power radPowFN1;
      SMD_MSR_Modelica.Units.Power radPowFN2;
      SMD_MSR_Modelica.Units.Power fissPowFN1;
      SMD_MSR_Modelica.Units.Power fissPowFN2;
      SMD_MSR_Modelica.Units.Power fissPowGN;
      SMD_MSR_Modelica.Units.Power decayPowFN1;
      SMD_MSR_Modelica.Units.Power decayPowFN2;
      input SMD_MSR_Modelica.PortsConnectors.VolumetircPowerIn decayHeat annotation(
        Placement(visible = true, transformation(origin = {-72, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-72, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.TempOut grapNode annotation(
        Placement(visible = true, transformation(origin = {40, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {40, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.TempOut fuelNode1 annotation(
        Placement(visible = true, transformation(origin = {40, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {40, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.TempIn temp_In annotation(
        Placement(visible = true, transformation(origin = {-72, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-72, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.PowerIn fissionPower annotation(
        Placement(visible = true, transformation(origin = {-72, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-72, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.TempOut fuelNode2 annotation(
        Placement(visible = true, transformation(origin = {40, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {40, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFractionIn fuelFlowFraction annotation(
        Placement(visible = true, transformation(origin = {-70, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      fuelNode1.T = TF1_0;
      fuelNode2.T = TF2_0;
      grapNode.T = TG_0;
    equation
      hA = hAnom*fuelFlowFraction.FF^(0.33);
      mdot_fuel = Vdot_fuelNom*rho_fuel*fuelFlowFraction.FF*regionFlowFrac;
      m_FN1 = vol_FN1*rho_fuel;
      m_FN2 = vol_FN2*rho_fuel;
      m_GN = vol_GN*rho_grap;
      powFN1 = m_FN1*cP_fuel*der(fuelNode1.T);
      powFN2 = m_FN2*cP_fuel*der(fuelNode2.T);
      powGN = m_GN*cP_grap*der(grapNode.T);
      flowPowFN1 = mdot_fuel*cP_fuel*(temp_In.T - fuelNode1.T);
      flowPowFN2 = mdot_fuel*cP_fuel*(fuelNode1.T - fuelNode2.T);
      condPowFN1 = ((KF*Ac)/LF1)*(temp_In.T - fuelNode1.T);
      condPowFN2 = ((KF*Ac)/LF2)*(fuelNode1.T - fuelNode2.T);
      convPowFN1 = hA*(kHT_FN1/(kHT_FN1 + kHT_FN2))*(fuelNode1.T - grapNode.T);
      convPowFN2 = hA*(kHT_FN2/(kHT_FN1 + kHT_FN2))*(fuelNode1.T - grapNode.T);
      if OutterRegion == true then
        radPowFN1 = e*SMD_MSR_Modelica.Constants.SigSBK*ArF1*((Tinf + 273.15)^4 - (fuelNode1.T + 273.15)^4);
        radPowFN2 = e*SMD_MSR_Modelica.Constants.SigSBK*ArF2*((Tinf + 273.15)^4 - (fuelNode1.T + 273.15)^4);
      else
        radPowFN1 = 0;
        radPowFN2 = 0;
      end if;
      fissPowFN1 = kFN1*fissionPower.P;
      fissPowFN2 = kFN2*fissionPower.P;
      fissPowGN = kG*fissionPower.P;
      decayPowFN1 = decayHeat.Q*vol_FN1;
      decayPowFN2 = decayHeat.Q*vol_FN2;
      powFN1 = flowPowFN1 + condPowFN1 - convPowFN1 + radPowFN1 + fissPowFN1 + decayPowFN1;
      powFN2 = flowPowFN2 + condPowFN2 - convPowFN2 + radPowFN2 + fissPowFN2 + decayPowFN2;
      powGN = convPowFN1 + convPowFN2 + fissPowGN;
      annotation(
        Diagram(coordinateSystem(extent = {{-100, 80}, {60, -100}}), graphics = {Text(origin = {9, 68}, extent = {{-47, 12}, {47, -12}}, textString = "CoreRegion"), Rectangle(origin = {-15.43, -10.01}, lineThickness = 2, extent = {{-75.81, 90.1}, {75.81, -90.1}}), Text(origin = {44, -76}, extent = {{-13, 13}, {13, -13}}, textString = "T_FN1"), Rectangle(origin = {-34, 20}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}}), Text(origin = {44, 24}, extent = {{-13, 13}, {13, -13}}, textString = "T_FN2"), Rectangle(origin = {-15.43, -10.01}, lineThickness = 2, extent = {{-75.81, 90.1}, {75.81, -90.1}}), Rectangle(origin = {-15.43, -10.01}, lineThickness = 2, extent = {{-75.81, 90.1}, {75.81, -90.1}}), Text(origin = {-71, 36}, extent = {{-6, 8}, {6, -8}}, textString = "FH"), Text(origin = {-71, -45}, extent = {{-6, 8}, {6, -8}}, textString = "FF"), Rectangle(origin = {11, -10}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-15, 50}, {15, -50}}), Rectangle(origin = {-34, -40}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}}), Text(origin = {-71, -5}, extent = {{-6, 8}, {6, -8}}, textString = "DH"), Text(origin = {-70, -86}, extent = {{-9, 10}, {9, -10}}, textString = "T_In"), Text(origin = {42, -25}, extent = {{-11, 12}, {11, -12}}, textString = "T_G1")}),
        Icon(coordinateSystem(extent = {{-100, 80}, {60, -100}}), graphics = {Text(origin = {-73, 36}, extent = {{-6, 8}, {6, -8}}, textString = "FH"), Text(origin = {-70, -86}, extent = {{-9, 10}, {9, -10}}, textString = "T_In"), Rectangle(origin = {-34, 20}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}}), Text(origin = {42, -25}, extent = {{-11, 12}, {11, -12}}, textString = "T_G1"), Text(origin = {44, -76}, extent = {{-13, 13}, {13, -13}}, textString = "T_FN1"), Text(origin = {44, 24}, extent = {{-13, 13}, {13, -13}}, textString = "T_FN2"), Rectangle(origin = {11, -10}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-15, 50}, {15, -50}}), Text(origin = {9, 68}, extent = {{-47, 12}, {47, -12}}, textString = "CoreRegion"), Rectangle(origin = {-15.43, -10.01}, lineThickness = 2, extent = {{-75.81, 90.1}, {75.81, -90.1}}), Text(origin = {-71, -45}, extent = {{-6, 8}, {6, -8}}, textString = "FF"), Text(origin = {-71, -5}, extent = {{-6, 8}, {6, -8}}, textString = "DH"), Rectangle(origin = {-15.43, -10.01}, lineThickness = 2, extent = {{-75.81, 90.1}, {75.81, -90.1}}), Rectangle(origin = {-15.43, -10.01}, lineThickness = 2, extent = {{-75.81, 90.1}, {75.81, -90.1}}), Rectangle(origin = {-34, -40}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}})}));
    end FuelChannel;

    model HeatExchanger
      parameter SMD_MSR_Modelica.Units.Volume vol_P;
      parameter SMD_MSR_Modelica.Units.Volume vol_T;
      parameter SMD_MSR_Modelica.Units.Volume vol_S;
      parameter SMD_MSR_Modelica.Units.Density rhoP;
      parameter SMD_MSR_Modelica.Units.Density rhoT;
      parameter SMD_MSR_Modelica.Units.Density rhoS;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_P;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_T;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_S;
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate VdotPnom;
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate VdotSnom;
      parameter SMD_MSR_Modelica.Units.Convection hApNom;
      parameter SMD_MSR_Modelica.Units.Convection hAsNom;
      parameter SMD_MSR_Modelica.Units.Conductivity Kp;
      parameter SMD_MSR_Modelica.Units.Conductivity Ks;
      parameter SMD_MSR_Modelica.Units.Area AcShell;
      parameter SMD_MSR_Modelica.Units.Area AcTube;
      parameter SMD_MSR_Modelica.Units.Length L_shell;
      parameter SMD_MSR_Modelica.Units.Length L_tube;
      parameter Boolean EnableRad;
      parameter SMD_MSR_Modelica.Units.Area ArShell;
      parameter SMD_MSR_Modelica.Units.Emissivity e;
      parameter SMD_MSR_Modelica.Units.Temperature Tinf;
      parameter SMD_MSR_Modelica.Units.Temperature TpIn_0;
      parameter SMD_MSR_Modelica.Units.Temperature TpOut_0;
      parameter SMD_MSR_Modelica.Units.Temperature TsIn_0;
      parameter SMD_MSR_Modelica.Units.Temperature TsOut_0;
      SMD_MSR_Modelica.Units.MassFlowRate mDotP;
      SMD_MSR_Modelica.Units.MassFlowRate mDotS;
      SMD_MSR_Modelica.Units.Mass m_PN;
      SMD_MSR_Modelica.Units.Mass m_TN;
      SMD_MSR_Modelica.Units.Mass m_SN;
      SMD_MSR_Modelica.Units.Volume volPN;
      SMD_MSR_Modelica.Units.Volume volTN;
      SMD_MSR_Modelica.Units.Volume volSN;
      SMD_MSR_Modelica.Units.Area Ar_PN;
      SMD_MSR_Modelica.Units.Length LpN;
      SMD_MSR_Modelica.Units.Length LsN;
      SMD_MSR_Modelica.Units.Convection hApn;
      SMD_MSR_Modelica.Units.Convection hAsn;
      SMD_MSR_Modelica.Units.Temperature T_PN1;
      SMD_MSR_Modelica.Units.Temperature T_PN2;
      SMD_MSR_Modelica.Units.Temperature T_PN3;
      SMD_MSR_Modelica.Units.Temperature T_TN1;
      SMD_MSR_Modelica.Units.Temperature T_TN2;
      SMD_MSR_Modelica.Units.Temperature T_SN1;
      SMD_MSR_Modelica.Units.Temperature T_SN2;
      SMD_MSR_Modelica.Units.Temperature T_SN3;
      SMD_MSR_Modelica.Units.Power powPN1;
      SMD_MSR_Modelica.Units.Power powPN2;
      SMD_MSR_Modelica.Units.Power powPN3;
      SMD_MSR_Modelica.Units.Power powPN4;
      SMD_MSR_Modelica.Units.Power powTN1;
      SMD_MSR_Modelica.Units.Power powTN2;
      SMD_MSR_Modelica.Units.Power powSN1;
      SMD_MSR_Modelica.Units.Power powSN2;
      SMD_MSR_Modelica.Units.Power powSN3;
      SMD_MSR_Modelica.Units.Power powSN4;
      SMD_MSR_Modelica.Units.Power flowPowPN1;
      SMD_MSR_Modelica.Units.Power flowPowPN2;
      SMD_MSR_Modelica.Units.Power flowPowPN3;
      SMD_MSR_Modelica.Units.Power flowPowPN4;
      SMD_MSR_Modelica.Units.Power flowPowSN1;
      SMD_MSR_Modelica.Units.Power flowPowSN2;
      SMD_MSR_Modelica.Units.Power flowPowSN3;
      SMD_MSR_Modelica.Units.Power flowPowSN4;
      SMD_MSR_Modelica.Units.Power condPowPN1;
      SMD_MSR_Modelica.Units.Power condPowPN2;
      SMD_MSR_Modelica.Units.Power condPowPN3;
      SMD_MSR_Modelica.Units.Power condPowPN4;
      SMD_MSR_Modelica.Units.Power condPowSN1;
      SMD_MSR_Modelica.Units.Power condPowSN2;
      SMD_MSR_Modelica.Units.Power condPowSN3;
      SMD_MSR_Modelica.Units.Power condPowSN4;
      SMD_MSR_Modelica.Units.Power convPowPN1;
      SMD_MSR_Modelica.Units.Power convPowPN2;
      SMD_MSR_Modelica.Units.Power convPowPN3;
      SMD_MSR_Modelica.Units.Power convPowPN4;
      SMD_MSR_Modelica.Units.Power convPowSN1;
      SMD_MSR_Modelica.Units.Power convPowSN2;
      SMD_MSR_Modelica.Units.Power convPowSN3;
      SMD_MSR_Modelica.Units.Power convPowSN4;
      SMD_MSR_Modelica.Units.Power radPowPN1;
      SMD_MSR_Modelica.Units.Power radPowPN2;
      SMD_MSR_Modelica.Units.Power radPowPN3;
      SMD_MSR_Modelica.Units.Power radPowPN4;
      SMD_MSR_Modelica.Units.Power decayPowPN1;
      SMD_MSR_Modelica.Units.Power decayPowPN2;
      SMD_MSR_Modelica.Units.Power decayPowPN3;
      SMD_MSR_Modelica.Units.Power decayPowPN4;
      input SMD_MSR_Modelica.PortsConnectors.TempIn T_in_pFluid annotation(
        Placement(visible = true, transformation(origin = {-100, 32}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.TempIn T_in_sFluid annotation(
        Placement(visible = true, transformation(origin = {163, -31}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {160, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.TempOut T_out_sFluid annotation(
        Placement(visible = true, transformation(origin = {-100, -38}, extent = {{-14, -14}, {14, 14}}, rotation = 0), iconTransformation(origin = {-100, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.TempOut T_out_pFluid annotation(
        Placement(visible = true, transformation(origin = {158, 28}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {160, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFractionIn primaryFF annotation(
        Placement(visible = true, transformation(origin = {-69, 49}, extent = {{-9, -9}, {9, 9}}, rotation = 0), iconTransformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFractionIn secondaryFF annotation(
        Placement(visible = true, transformation(origin = {123, -47}, extent = {{-11, -11}, {11, 11}}, rotation = 0), iconTransformation(origin = {130, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.VolumetircPowerIn P_decay annotation(
        Placement(visible = true, transformation(origin = {30, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {30, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      T_PN1 = TpIn_0 - (TpIn_0 - TpOut_0)/4;
      T_PN2 = TpIn_0 - 2*(TpIn_0 - TpOut_0)/4;
      T_PN3 = TpIn_0 - 3*(TpIn_0 - TpOut_0)/4;
      T_out_pFluid.T = TpOut_0;
      T_TN1 = (T_PN1*hApn + T_SN3*hAsn)/(hApn + hAsn);
      T_TN2 = (T_PN3*hApn + T_SN1*hAsn)/(hApn + hAsn);
      T_SN1 = TpIn_0 + (TsOut_0 - TpIn_0);
      T_SN2 = TpIn_0 + 2*(TsOut_0 - TpIn_0);
      T_SN3 = TpIn_0 + 3*(TsOut_0 - TpIn_0);
      T_out_sFluid.T = TsOut_0;
    equation
      mDotP = VdotPnom*rhoP*primaryFF.FF;
      mDotS = VdotSnom*rhoS*secondaryFF.FF;
      hApn = 0.07797E6*primaryFF.FF^(0.33);
      hAsn = (1 - 0.01)*1E6*hAsNom*secondaryFF.FF + 0.01E6*hAsNom;
      volPN = vol_P/4;
      volTN = vol_T/2;
      volSN = vol_S/4;
      Ar_PN = ArShell/4;
      LpN = L_shell/4;
      LsN = L_tube/4;
      m_PN = volPN*rhoP;
      m_TN = volTN*rhoT;
      m_SN = volSN*rhoS;
      powPN1 = m_PN*cP_P*der(T_PN1);
      powPN2 = m_PN*cP_P*der(T_PN2);
      powPN3 = m_PN*cP_P*der(T_PN3);
      powPN4 = m_PN*cP_P*der(T_out_pFluid.T);
      powTN1 = m_TN*cP_T*der(T_TN1);
      powTN2 = m_TN*cP_T*der(T_TN2);
      powSN1 = m_SN*cP_S*der(T_SN1);
      powSN2 = m_SN*cP_S*der(T_SN2);
      powSN3 = m_SN*cP_S*der(T_SN3);
      powSN4 = m_SN*cP_S*der(T_out_sFluid.T);
      flowPowPN1 = mDotP*cP_P*(T_in_pFluid.T - T_PN1);
      flowPowPN2 = mDotP*cP_P*(T_PN1 - T_PN2);
      flowPowPN3 = mDotP*cP_P*(T_PN2 - T_PN3);
      flowPowPN4 = mDotP*cP_P*(T_PN3 - T_out_pFluid.T);
      flowPowSN1 = mDotS*cP_S*(T_in_sFluid.T - T_SN1);
      flowPowSN2 = mDotS*cP_S*(T_SN1 - T_SN2);
      flowPowSN3 = mDotS*cP_S*(T_SN2 - T_SN3);
      flowPowSN4 = mDotS*cP_S*(T_SN3 - T_out_sFluid.T);
      condPowPN1 = ((Kp*AcShell)/LpN)*(T_in_pFluid.T - T_PN1);
      condPowPN2 = ((Kp*AcShell)/LpN)*(T_PN1 - T_PN2);
      condPowPN3 = ((Kp*AcShell)/LpN)*(T_PN2 - T_PN3);
      condPowPN4 = ((Kp*AcShell)/LpN)*(T_PN3 - T_out_pFluid.T);
      condPowSN1 = ((Ks*AcTube)/LsN)*(T_in_sFluid.T - T_SN1);
      condPowSN2 = ((Ks*AcTube)/LsN)*(T_SN1 - T_SN2);
      condPowSN3 = ((Ks*AcTube)/LsN)*(T_SN2 - T_SN3);
      condPowSN4 = ((Ks*AcTube)/LsN)*(T_SN3 - T_out_sFluid.T);
      convPowPN1 = hApn*(T_PN1 - T_TN1);
      convPowPN2 = hApn*(T_PN1 - T_TN1);
      convPowPN3 = hApn*(T_PN3 - T_TN2);
      convPowPN4 = hApn*(T_PN3 - T_TN2);
      convPowSN1 = hAsn*(T_TN2 - T_SN1);
      convPowSN2 = hAsn*(T_TN2 - T_SN1);
      convPowSN3 = hAsn*(T_TN1 - T_SN3);
      convPowSN4 = hAsn*(T_TN1 - T_SN3);
      if EnableRad == true then
        radPowPN1 = e*SMD_MSR_Modelica.Constants.SigSBK*Ar_PN*((Tinf + 273.15)^4 - (T_PN1 + 273.15)^4);
        radPowPN2 = e*SMD_MSR_Modelica.Constants.SigSBK*Ar_PN*((Tinf + 273.15)^4 - (T_PN2 + 273.15)^4);
        radPowPN3 = e*SMD_MSR_Modelica.Constants.SigSBK*Ar_PN*((Tinf + 273.15)^4 - (T_PN3 + 273.15)^4);
        radPowPN4 = e*SMD_MSR_Modelica.Constants.SigSBK*Ar_PN*((Tinf + 273.15)^4 - (T_out_pFluid.T + 273.15)^4);
      else
        radPowPN1 = 0;
        radPowPN2 = 0;
        radPowPN3 = 0;
        radPowPN4 = 0;
      end if;
      decayPowPN1 = P_decay.Q*volPN;
      decayPowPN2 = P_decay.Q*volPN;
      decayPowPN3 = P_decay.Q*volPN;
      decayPowPN4 = P_decay.Q*volPN;
      powPN1 = flowPowPN1 + condPowPN1 - convPowPN1 + radPowPN1 + decayPowPN1;
      powPN2 = flowPowPN2 + condPowPN2 - convPowPN2 + radPowPN2 + decayPowPN2;
      powPN3 = flowPowPN3 + condPowPN3 - convPowPN3 + radPowPN3 + decayPowPN3;
      powPN4 = flowPowPN4 + condPowPN4 - convPowPN4 + radPowPN4 + decayPowPN4;
      powTN1 = convPowPN1 + convPowPN2 - convPowSN3 - convPowSN4;
      powTN2 = convPowPN3 + convPowPN4 - convPowSN1 - convPowSN2;
      powSN1 = flowPowSN1 + convPowSN1;
      powSN2 = flowPowSN2 + convPowSN2;
      powSN3 = flowPowSN3 + convPowSN3;
      powSN4 = flowPowSN4 + convPowSN4;
      annotation(
        Diagram(graphics = {Rectangle(origin = {30.2864, 0.0340212}, extent = {{-149.242, 59.993}, {149.242, -59.993}}), Rectangle(origin = {-60, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {120, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {0, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-59.95, -29.61}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.63, -9.61}, {19.63, 9.61}}), Rectangle(origin = {-0.12, -29.82}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.88, 9.82}, {19.88, -9.82}}), Rectangle(origin = {59.8, -29.85}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.8, 9.85}, {19.8, -9.85}}), Rectangle(origin = {119.85, -29.89}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.85, 9.89}, {19.85, -9.89}}), Rectangle(origin = {89.92, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {60, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-30.08, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Text(origin = {150, 49}, extent = {{-28, 9}, {28, -9}}, textString = "HX")}, coordinateSystem(extent = {{-120, 60}, {180, -60}})),
        Icon(graphics = {Rectangle(origin = {30.29, 0.03}, lineThickness = 2, extent = {{-149.24, 59.99}, {149.24, -59.99}}), Rectangle(origin = {-30.08, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {60, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {89.92, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {-0.12, -29.82}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.88, 9.82}, {19.88, -9.82}}), Rectangle(origin = {0, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-60, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {120, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {59.8, -29.85}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.8, 9.85}, {19.8, -9.85}}), Rectangle(origin = {119.85, -29.89}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.85, 9.89}, {19.85, -9.89}}), Rectangle(origin = {-59.95, -29.61}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.63, -9.61}, {19.63, 9.61}}), Text(origin = {150, 49}, extent = {{-28, 9}, {28, -9}}, textString = "HX")}, coordinateSystem(extent = {{-120, 60}, {180, -60}})));
    end HeatExchanger;

    model MSRR9R
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
      SMD_MSR_Modelica.Nuclear.mPKE mpke(numGroups = numGroups, lambda = lambda, beta = beta, LAMBDA = LAMBDA, n_0 = n_0, addRho0 = addRho0, nomTauLoop = 8.7006, nomTauCore = 34.8025) annotation(
        Placement(transformation(origin = {174.6, 322.2}, extent = {{-63.6, -63.6}, {42.4, 42.4}}, rotation = -90)));
      SMD_MSR_Modelica.Nuclear.PowerBlock powerblock(P(displayUnit = "W") = 1E6, TotalFuelVol = 0.5) annotation(
        Placement(transformation(origin = {-92.4, 330.2}, extent = {{29.6, -59.2}, {-44.4, 14.8}}, rotation = -90)));
      SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(numGroups = 3, DHYG = {2.3751e-03, 8.7763e-05, 1.9596e-06}, DHlamG = {0.09453, 0.004420, 8.6098E-5}) annotation(
        Placement(transformation(origin = {-18.4, 480.4}, extent = {{-29.6, 29.6}, {44.4, -44.4}}, rotation = 90)));
      SMD_MSR_Modelica.Signals.Constants.ConstantNeutronSource constantNeutronSource(S = 0) annotation(
        Placement(transformation(origin = {390, 430}, extent = {{-27, -27}, {27, 27}}, rotation = 90)));
      MSRR.Components.FuelChannel R1(vol_FN1 = volF1[1], vol_FN2 = volF2[1], vol_GN = volG[1], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[1], kFN2 = kFN2[1], kG = kHT1[1] + kHT2[1], hAnom = hA[1], kHT_FN1 = kHT1[1], kHT_FN2 = kHT2[1], TF1_0 = TF1_0[1], TF2_0 = TF2_0[1], TG_0 = TG_0[1], regionFlowFrac = flowFracRegions[1], KF = kFuel, Ac = Ac[1], LF1 = LF1[1], LF2 = LF2[1], OutterRegion = false, ArF1 = ArF1[1], ArF2 = ArF2[1], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {-450.972, -402.111}, extent = {{-75.1389, -75.1389}, {45.0833, 60.1111}})));
      MSRR.Components.FuelChannel R2(vol_FN1 = volF1[2], vol_FN2 = volF2[2], vol_GN = volG[2], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[2], kFN2 = kFN2[2], kG = kHT1[2] + kHT2[2], hAnom = hA[2], kHT_FN1 = kHT1[2], kHT_FN2 = kHT2[2], TF1_0 = TF1_0[2], TF2_0 = TF2_0[2], TG_0 = TG_0[2], regionFlowFrac = flowFracRegions[2], KF = kFuel, Ac = Ac[2], LF1 = LF1[2], LF2 = LF2[2], OutterRegion = false, ArF1 = ArF1[2], ArF2 = ArF2[2], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {-40.25, -600.333}, extent = {{-75.4167, -75.4167}, {45.25, 60.3333}})));
      MSRR.Components.FuelChannel R3(vol_FN1 = volF1[3], vol_FN2 = volF2[3], vol_GN = volG[3], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[3], kFN2 = kFN2[3], kG = kHT1[3] + kHT2[3], hAnom = hA[3], kHT_FN1 = kHT1[3], kHT_FN2 = kHT2[3], TF1_0 = TF1_0[3], TF2_0 = TF2_0[3], TG_0 = TG_0[3], regionFlowFrac = flowFracRegions[2], KF = kFuel, Ac = Ac[2], LF1 = LF1[3], LF2 = LF2[3], OutterRegion = false, ArF1 = ArF1[3], ArF2 = ArF2[3], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {-40.3333, -400.944}, extent = {{-75.5556, -75.5556}, {45.3333, 60.4444}})));
      MSRR.Components.FuelChannel R4(vol_FN1 = volF1[4], vol_FN2 = volF2[4], vol_GN = volG[4], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[4], kFN2 = kFN2[4], kG = kHT1[4] + kHT2[4], hAnom = hA[4], kHT_FN1 = kHT1[4], kHT_FN2 = kHT2[4], TF1_0 = TF1_0[4], TF2_0 = TF2_0[4], TG_0 = TG_0[4], regionFlowFrac = flowFracRegions[2], KF = kFuel, Ac = Ac[2], LF1 = LF1[4], LF2 = LF2[4], OutterRegion = false, ArF1 = ArF1[4], ArF2 = ArF2[4], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {-42.2778, -200.278}, extent = {{-74.7222, -74.7222}, {44.8333, 59.7778}})));
      MSRR.Components.FuelChannel R5(vol_FN1 = volF1[5], vol_FN2 = volF2[5], vol_GN = volG[5], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[5], kFN2 = kFN2[5], kG = kHT1[5] + kHT2[5], hAnom = hA[5], kHT_FN1 = kHT1[5], kHT_FN2 = kHT2[5], TF1_0 = TF1_0[5], TF2_0 = TF2_0[5], TG_0 = TG_0[5], regionFlowFrac = flowFracRegions[3], KF = kFuel, Ac = Ac[3], LF1 = LF1[5], LF2 = LF2[5], OutterRegion = false, ArF1 = ArF1[5], ArF2 = ArF2[5], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {330.167, -600.833}, extent = {{-74.1667, -74.1667}, {44.5, 59.3333}})));
      MSRR.Components.FuelChannel R6(vol_FN1 = volF1[6], vol_FN2 = volF2[6], vol_GN = volG[6], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[6], kFN2 = kFN2[6], kG = kHT1[6] + kHT2[6], hAnom = hA[6], kHT_FN1 = kHT1[6], kHT_FN2 = kHT2[6], TF1_0 = TF1_0[6], TF2_0 = TF2_0[6], TG_0 = TG_0[6], regionFlowFrac = flowFracRegions[3], KF = kFuel, Ac = Ac[3], LF1 = LF1[6], LF2 = LF2[6], OutterRegion = false, ArF1 = ArF1[6], ArF2 = ArF2[6], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {330.139, -400.861}, extent = {{-75.1389, -75.1389}, {45.0833, 60.1111}})));
      MSRR.Components.FuelChannel R7(vol_FN1 = volF1[7], vol_FN2 = volF2[7], vol_GN = volG[7], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[7], kFN2 = kFN2[7], kG = kHT1[7] + kHT2[7], hAnom = hA[7], kHT_FN1 = kHT1[7], kHT_FN2 = kHT2[7], TF1_0 = TF1_0[7], TF2_0 = TF2_0[7], TG_0 = TG_0[7], regionFlowFrac = flowFracRegions[3], KF = kFuel, Ac = Ac[3], LF1 = LF1[7], LF2 = LF2[7], OutterRegion = false, ArF1 = ArF1[7], ArF2 = ArF2[7], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {330.611, -200.389}, extent = {{-73.6111, -73.6111}, {44.1667, 58.8889}})));
      MSRR.Components.FuelChannel R8(vol_FN1 = volF1[8], vol_FN2 = volF2[8], vol_GN = volG[8], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[8], kFN2 = kFN2[8], kG = kHT1[8] + kHT2[8], hAnom = hA[8], kHT_FN1 = kHT1[8], kHT_FN2 = kHT2[8], TF1_0 = TF1_0[8], TF2_0 = TF2_0[8], TG_0 = TG_0[8], regionFlowFrac = flowFracRegions[4], KF = kFuel, Ac = Ac[4], LF1 = LF1[8], LF2 = LF2[8], OutterRegion = EnableRad, ArF1 = ArF1[8], ArF2 = ArF2[8], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {732.861, -514.889}, extent = {{-74.8611, -74.8611}, {44.9167, 59.8889}})));
      MSRR.Components.FuelChannel R9(vol_FN1 = volF1[9], vol_FN2 = volF2[9], vol_GN = volG[9], rho_fuel = rho_fuel, rho_grap = rho_grap, cP_fuel = cP_fuel, cP_grap = cP_grap, Vdot_fuelNom = volDotFuel, kFN1 = kFN1[9], kFN2 = kFN2[9], kG = kHT1[9] + kHT2[9], hAnom = hA[9], kHT_FN1 = kHT1[9], kHT_FN2 = kHT2[9], TF1_0 = TF1_0[9], TF2_0 = TF2_0[9], TG_0 = TG_0[9], regionFlowFrac = flowFracRegions[4], KF = kFuel, Ac = Ac[4], LF1 = LF1[9], LF2 = LF2[9], OutterRegion = EnableRad, ArF1 = ArF1[9], ArF2 = ArF2[9], e = e, Tinf = T_inf) annotation(
        Placement(transformation(origin = {710.833, -320.889}, extent = {{-73.6111, -73.6111}, {44.1667, 58.8889}})));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF1(a_F = aF, a_G = aG, IF1 = IF1[1], IF2 = IF2[1], IG = IG[1], FuelTempSetPointNode1 = TF1_0[1], FuelTempSetPointNode2 = TF2_0[1], GrapTempSetPoint = TG_0[1]) annotation(
        Placement(transformation(origin = {-280.4, -401.6}, extent = {{-69.6, -46.4}, {46.4, 69.6}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF2(a_F = aF, a_G = aG, IF1 = IF1[2], IF2 = IF2[2], IG = IG[2], FuelTempSetPointNode1 = TF1_0[2], FuelTempSetPointNode2 = TF2_0[2], GrapTempSetPoint = TG_0[2]) annotation(
        Placement(transformation(origin = {137.6, -596.4}, extent = {{-69.6, -46.4}, {46.4, 69.6}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF3(a_F = aF, a_G = aG, IF1 = IF1[3], IF2 = IF2[3], IG = IG[3], FuelTempSetPointNode1 = TF1_0[3], FuelTempSetPointNode2 = TF2_0[3], GrapTempSetPoint = TG_0[3]) annotation(
        Placement(transformation(origin = {137.6, -400}, extent = {{-69.6, -46.4}, {46.4, 69.6}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF4(a_F = aF, a_G = aG, IF1 = IF1[4], IF2 = IF2[4], IG = IG[4], FuelTempSetPointNode1 = TF1_0[4], FuelTempSetPointNode2 = TF2_0[4], GrapTempSetPoint = TG_0[4]) annotation(
        Placement(transformation(origin = {135.6, -200}, extent = {{-69.6, -46.4}, {46.4, 69.6}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF5(a_F = aF, a_G = aG, IF1 = IF1[5], IF2 = IF2[5], IG = IG[5], FuelTempSetPointNode1 = TF1_0[5], FuelTempSetPointNode2 = TF2_0[5], GrapTempSetPoint = TG_0[5]) annotation(
        Placement(transformation(origin = {510.6, -606}, extent = {{-69.6, -46.4}, {46.4, 69.6}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF6(a_F = aF, a_G = aG, IF1 = IF1[6], IF2 = IF2[6], IG = IG[6], FuelTempSetPointNode1 = TF1_0[6], FuelTempSetPointNode2 = TF2_0[6], GrapTempSetPoint = TG_0[6]) annotation(
        Placement(transformation(origin = {510.6, -400}, extent = {{-69.6, -46.4}, {46.4, 69.6}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF7(a_F = aF, a_G = aG, IF1 = IF1[7], IF2 = IF2[7], IG = IG[7], FuelTempSetPointNode1 = TF1_0[7], FuelTempSetPointNode2 = TF2_0[7], GrapTempSetPoint = TG_0[7]) annotation(
        Placement(transformation(origin = {510.6, -200}, extent = {{-69.6, -46.4}, {46.4, 69.6}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF8(a_F = aF, a_G = aG, IF1 = IF1[8], IF2 = IF2[8], IG = IG[8], FuelTempSetPointNode1 = TF1_0[8], FuelTempSetPointNode2 = TF2_0[8], GrapTempSetPoint = TG_0[8]) annotation(
        Placement(transformation(origin = {900.6, -510}, extent = {{-69.6, -46.4}, {46.4, 69.6}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback RF9(a_F = aF, a_G = aG, IF1 = IF1[9], IF2 = IF2[9], IG = IG[9], FuelTempSetPointNode1 = TF1_0[9], FuelTempSetPointNode2 = TF2_0[9], GrapTempSetPoint = TG_0[9]) annotation(
        Placement(transformation(origin = {904.6, -320}, extent = {{-69.6, -46.4}, {46.4, 69.6}}, rotation = 90)));
      SMD_MSR_Modelica.Nuclear.SumReactivity sumFB(numInput = 9) annotation(
        Placement(transformation(origin = {516.544, -5.0561}, extent = {{-37.5439, 56.3158}, {56.3158, -37.5439}})));
      SMD_MSR_Modelica.HeatTransport.FlowDistributor flowDistributor(numOutput = 4, freeConvectionFF = freeConvFF, coastDownK = regionCoastDownK, regionTripTime = regionTripTime) annotation(
        Placement(transformation(origin = {395.602, -954.4}, extent = {{30.3983, 91.195}, {-91.195, -30.3983}})));
      SMD_MSR_Modelica.HeatTransport.MixingPot upperPlenum(numStreams = 4, vol = volUP, VdotNom = volDotFuel, flowFractionsNom = flowFracRegions, rho = rho_fuel, Cp = cP_fuel, K = kFuel, Ac = 1.6118, L = 0.093874, Ar = 0.4225, e = e, T_0 = Tmix_0, Tinf = T_inf, EnableRad = EnableRad) annotation(
        Placement(transformation(origin = {450.182, 131.988}, extent = {{-51.443, 38.5821}, {25.7215, -38.5821}}, rotation = 90)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFractionIn flowFracIn annotation(
        Placement(transformation(origin = {329, -1071}, extent = {{-37, -37}, {37, 37}}), iconTransformation(origin = {-200, -80}, extent = {{-18, -18}, {18, 18}})));
      input SMD_MSR_Modelica.PortsConnectors.TempIn tempIn annotation(
        Placement(transformation(origin = {217, -827}, extent = {{-21, -21}, {21, 21}}), iconTransformation(origin = {-200, -158}, extent = {{-18, -18}, {18, 18}})));
      output SMD_MSR_Modelica.PortsConnectors.TempOut tempOut annotation(
        Placement(transformation(origin = {1117, 235}, extent = {{-28, -28}, {28, 28}}), iconTransformation(origin = {-120, -80}, extent = {{-18, -18}, {18, 18}})));
      output SMD_MSR_Modelica.PortsConnectors.VolumetircPowerOut volumetircPowerOut annotation(
        Placement(transformation(origin = {-144, 522}, extent = {{13, 13}, {-13, -13}}), iconTransformation(origin = {-119, -159}, extent = {{-17, -17}, {17, 17}})));
      input SMD_MSR_Modelica.PortsConnectors.RealIn realIn annotation(
        Placement(transformation(origin = {290, 450}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-160, -80}, extent = {{-18, -18}, {18, 18}})));
    equation
      connect(R1.fuelNode1, RF1.fuelNode1) annotation(
        Line(points = {{-420.916, -447.194}, {-371.416, -447.194}, {-371.416, -448.194}, {-321.916, -448.194}}, color = {204, 0, 0}, thickness = 1));
      connect(R1.fuelNode2, RF1.fuelNode2) annotation(
        Line(points = {{-420.916, -370.553}, {-371.416, -370.553}, {-371.416, -377.553}, {-321.916, -377.553}}, color = {204, 0, 0}, thickness = 1));
      connect(R1.grapNode, RF1.grapNode) annotation(
        Line(points = {{-420.916, -409.625}, {-392.916, -409.625}, {-392.916, -412.625}, {-321.916, -412.625}}, color = {204, 0, 0}, thickness = 1));
      connect(R2.fuelNode1, RF2.fuelNode1) annotation(
        Line(points = {{-10, -646}, {36, -646}, {36, -643}, {91, -643}}, color = {204, 0, 0}, thickness = 1));
      connect(R2.grapNode, RF2.grapNode) annotation(
        Line(points = {{-10, -608}, {91, -608}}, color = {204, 0, 0}, thickness = 1));
      connect(R2.fuelNode2, RF2.fuelNode2) annotation(
        Line(points = {{-10, -568}, {40, -568}, {40, -573}, {91, -573}}, color = {204, 0, 0}, thickness = 1));
      connect(R3.fuelNode1, RF3.fuelNode1) annotation(
        Line(points = {{-10, -446}, {91, -446}}, color = {204, 0, 0}, thickness = 1));
      connect(R3.grapNode, RF3.grapNode) annotation(
        Line(points = {{-10, -408}, {40.5, -408}, {40.5, -412}, {91, -412}}, color = {204, 0, 0}, thickness = 1));
      connect(R3.fuelNode2, RF3.fuelNode2) annotation(
        Line(points = {{-10, -370}, {44, -370}, {44, -377}, {91, -377}}, color = {204, 0, 0}, thickness = 1));
      connect(R4.fuelNode1, RF4.fuelNode1) annotation(
        Line(points = {{-27, -253}, {-27, -258}, {78, -258}}, color = {204, 0, 0}, thickness = 1));
      connect(R4.grapNode, RF4.grapNode) annotation(
        Line(points = {{-27, -215}, {88, -215}, {88, -223}, {78, -223}}, color = {204, 0, 0}, thickness = 1));
      connect(R4.fuelNode2, RF4.fuelNode2) annotation(
        Line(points = {{-27, -176}, {-27, -188}, {78, -188}}, color = {204, 0, 0}, thickness = 1));
      connect(R5.fuelNode1, RF5.fuelNode1) annotation(
        Line(points = {{359.834, -645.333}, {399.334, -645.333}, {399.334, -651.333}, {466.834, -651.333}}, color = {204, 0, 0}, thickness = 1));
      connect(R5.grapNode, RF5.grapNode) annotation(
        Line(points = {{359.834, -608.25}, {466.834, -608.25}, {466.834, -618.25}}, color = {204, 0, 0}, thickness = 1));
      connect(R5.fuelNode2, RF5.fuelNode2) annotation(
        Line(points = {{359.834, -569.683}, {466.834, -569.683}, {466.834, -582.683}}, color = {204, 0, 0}, thickness = 1));
      connect(R6.fuelNode1, RF6.fuelNode1) annotation(
        Line(points = {{360.195, -445.944}, {471.195, -445.944}}, color = {204, 0, 0}, thickness = 1));
      connect(R6.grapNode, RF6.grapNode) annotation(
        Line(points = {{360.195, -408.375}, {471.195, -408.375}, {471.195, -412.375}}, color = {204, 0, 0}, thickness = 1));
      connect(R6.fuelNode2, RF6.fuelNode2) annotation(
        Line(points = {{360.195, -369.303}, {471.195, -369.303}, {471.195, -376.303}}, color = {204, 0, 0}, thickness = 1));
      connect(R7.fuelNode1, RF7.fuelNode1) annotation(
        Line(points = {{345, -252}, {345, -246.556}, {472.055, -246.556}}, color = {204, 0, 0}, thickness = 1));
      connect(R7.grapNode, RF7.grapNode) annotation(
        Line(points = {{345, -215}, {345, -211.75}, {472.055, -211.75}}, color = {204, 0, 0}, thickness = 1));
      connect(R7.fuelNode2, RF7.fuelNode2) annotation(
        Line(points = {{345, -177}, {345, -175.472}, {472.055, -175.472}}, color = {204, 0, 0}, thickness = 1));
      connect(R8.fuelNode1, RF8.fuelNode1) annotation(
        Line(points = {{748, -567}, {803, -567}, {803, -568}, {843, -568}}, color = {204, 0, 0}, thickness = 1));
      connect(R8.grapNode, RF8.grapNode) annotation(
        Line(points = {{748, -530}, {802, -530}, {802, -533}, {843, -533}}, color = {204, 0, 0}, thickness = 1));
      connect(R8.fuelNode2, RF8.fuelNode2) annotation(
        Line(points = {{748, -491}, {805, -491}, {805, -498}, {843, -498}}, color = {204, 0, 0}, thickness = 1));
      connect(R9.fuelNode1, RF9.fuelNode1) annotation(
        Line(points = {{740.277, -365.056}, {798.277, -365.056}, {798.277, -378}, {847, -378}}, color = {204, 0, 0}, thickness = 1));
      connect(R9.grapNode, RF9.grapNode) annotation(
        Line(points = {{740.277, -328.25}, {740.277, -343}, {847, -343}}, color = {204, 0, 0}, thickness = 1));
      connect(R9.fuelNode2, RF9.fuelNode2) annotation(
        Line(points = {{740.277, -289.972}, {740.277, -308}, {847, -308}}, color = {204, 0, 0}, thickness = 1));
      connect(R2.fuelNode2, R3.temp_In) annotation(
        Line(points = {{-10, -568}, {-12, -568}, {-12, -508}, {-94, -508}, {-94, -456}}, color = {204, 0, 0}, thickness = 1));
      connect(R3.fuelNode2, R4.temp_In) annotation(
        Line(points = {{-10, -370}, {-10, -306}, {-111, -306}, {-111, -262}}, color = {204, 0, 0}, thickness = 1));
      connect(R6.fuelNode2, R7.temp_In) annotation(
        Line(points = {{360.195, -369.303}, {360.195, -305.303}, {263, -305.303}, {263, -261}}, color = {204, 0, 0}, thickness = 1));
      connect(R5.fuelNode2, R6.temp_In) annotation(
        Line(points = {{359.834, -569.683}, {357.834, -569.683}, {357.834, -509.683}, {275.834, -509.683}, {275.834, -453.683}}, color = {204, 0, 0}, thickness = 1));
      connect(R8.fuelNode2, R9.temp_In) annotation(
        Line(points = {{748, -491}, {748, -431.447}, {658.805, -431.447}, {658.805, -373.447}}, color = {204, 0, 0}, thickness = 1));
      connect(RF1.feedback, sumFB.reactivityIn[1]) annotation(
        Line(points = {{-257, -413}, {-174, -413}, {-174, -14}, {533, -14}}, color = {78, 154, 6}, thickness = 1));
      connect(RF4.feedback, sumFB.reactivityIn[4]) annotation(
        Line(points = {{147, -223}, {162, -223}, {162, -14}, {533, -14}}, color = {78, 154, 6}, thickness = 1));
      connect(RF3.feedback, sumFB.reactivityIn[3]) annotation(
        Line(points = {{160, -412}, {164, -412}, {164, -14}, {533, -14}}, color = {78, 154, 6}, thickness = 1));
      connect(RF2.feedback, sumFB.reactivityIn[2]) annotation(
        Line(points = {{160, -608}, {160, -14}, {533, -14}}, color = {78, 154, 6}, thickness = 1));
      connect(RF7.feedback, sumFB.reactivityIn[7]) annotation(
        Line(points = {{534, -212}, {534, -111}, {533, -111}, {533, -14}}, color = {78, 154, 6}, thickness = 1));
      connect(RF6.feedback, sumFB.reactivityIn[6]) annotation(
        Line(points = {{534, -412}, {534, -211}, {533, -211}, {533, -14}}, color = {78, 154, 6}, thickness = 1));
      connect(RF5.feedback, sumFB.reactivityIn[5]) annotation(
        Line(points = {{534, -618}, {534, -314}, {533, -314}, {533, -14}}, color = {78, 154, 6}, thickness = 1));
      connect(RF9.feedback, sumFB.reactivityIn[9]) annotation(
        Line(points = {{928, -338}, {916, -9}, {533, -14}}, color = {78, 154, 6}, thickness = 1));
      connect(RF8.feedback, sumFB.reactivityIn[8]) annotation(
        Line(points = {{924, -550}, {912, -31}, {533, -14}}, color = {78, 154, 6}, thickness = 1));
      connect(sumFB.reactivityOut, mpke.feedback) annotation(
        Line(points = {{535, 42}, {535, 333}, {185, 333}}, color = {78, 154, 6}, thickness = 1));
      connect(constantNeutronSource.neutronEmsRateOut, mpke.S) annotation(
        Line(points = {{380.82, 429.46}, {380.82, 354}, {185, 354}}, color = {191, 64, 188}, thickness = 1));
      connect(mpke.n_population, decayHeat.nPop) annotation(
        Line(points = {{122, 343}, {98.4, 343}, {98.4, 494}, {19, 494}}, color = {20, 36, 248}, thickness = 1));
      connect(decayHeat.decayHeat_Out, powerblock.decayNP) annotation(
        Line(points = {{-26, 495}, {-26, 367}, {-115, 367}}, color = {0, 225, 255}, thickness = 1));
      connect(mpke.n_population, powerblock.nPop) annotation(
        Line(points = {{122, 343}, {8.5, 343}, {8.5, 323}, {-115, 323}}, color = {20, 36, 248}, thickness = 1));
      connect(realIn, mpke.ReactivityIn) annotation(
        Line(points = {{290, 450}, {290, 442}, {185, 442}, {185, 375}}, thickness = 1));
      connect(tempIn, R2.temp_In) annotation(
        Line(points = {{217, -827}, {-94, -827}, {-94, -654}}, color = {204, 0, 0}, thickness = 1));
      connect(tempIn, R1.temp_In) annotation(
        Line(points = {{217, -827}, {-505, -827}, {-505, -456}}, color = {204, 0, 0}, thickness = 1));
      connect(tempIn, R5.temp_In) annotation(
        Line(points = {{217, -827}, {277, -827}, {277, -654}}, color = {204, 0, 0}, thickness = 1));
      connect(powerblock.fissionPower, R1.fissionPower) annotation(
        Line(points = {{-159, 323}, {-159, 181}, {-505, 181}, {-505, -360}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R4.fissionPower) annotation(
        Line(points = {{-159, 323}, {-142, 323}, {-142, -166}, {-111, -166}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R3.fissionPower) annotation(
        Line(points = {{-159, 323}, {-142, 323}, {-142, -358}, {-94, -358}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R2.fissionPower) annotation(
        Line(points = {{-159, 323}, {-159, -558}, {-94, -558}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R7.fissionPower) annotation(
        Line(points = {{-159, 323}, {-140, 323}, {-140, -96}, {208, -96}, {208, -167}, {263, -167}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R6.fissionPower) annotation(
        Line(points = {{-159, 323}, {-138, 323}, {-138, -98}, {206, -98}, {206, -359}, {276, -359}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R5.fissionPower) annotation(
        Line(points = {{-159, 323}, {-140, 323}, {-140, -100}, {204, -100}, {204, -559}, {277, -559}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R9.fissionPower) annotation(
        Line(points = {{-159, 323}, {-140, 323}, {-140, -96}, {610, -96}, {610, -280}, {658, -280}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.fissionPower, R8.fissionPower) annotation(
        Line(points = {{-159, 323}, {-140, 323}, {-140, -98}, {610, -98}, {610, -477}, {609.5, -477}, {609.5, -480}, {664, -480}}, color = {129, 61, 156}, thickness = 1));
      connect(powerblock.decayPowerM, R1.decayHeat) annotation(
        Line(points = {{-159, 367}, {-546, 367}, {-546, -393}, {-505, -393}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R4.decayHeat) annotation(
        Line(points = {{-159, 367}, {-198, 367}, {-198, -199}, {-111, -199}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R3.decayHeat) annotation(
        Line(points = {{-159, 367}, {-198, 367}, {-198, -392}, {-94, -392}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R2.decayHeat) annotation(
        Line(points = {{-159, 367}, {-198, 367}, {-198, -592}, {-94, -592}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R7.decayHeat) annotation(
        Line(points = {{-159, 367}, {-198, 367}, {-198, -50}, {226, -50}, {226, -199}, {263, -199}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R6.decayHeat) annotation(
        Line(points = {{-159, 367}, {-198, 367}, {-198, -50}, {224, -50}, {224, -392}, {276, -392}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R5.decayHeat) annotation(
        Line(points = {{-159, 367}, {-196, 367}, {-196, -50}, {224, -50}, {224, -592}, {276, -592}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R9.decayHeat) annotation(
        Line(points = {{-159, 367}, {-196, 367}, {-196, -52}, {590, -52}, {590, -312}, {658, -312}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, R8.decayHeat) annotation(
        Line(points = {{-159, 367}, {-198, 367}, {-198, -54}, {590, -54}, {590, -513}, {664, -513}}, color = {220, 138, 221}, thickness = 1));
      connect(powerblock.decayPowerM, volumetircPowerOut) annotation(
        Line(points = {{-159, 367}, {-159, 489.6}, {-144, 489.6}, {-144, 522}}, color = {220, 138, 221}, thickness = 1));
      connect(flowFracIn, flowDistributor.flowFracIn) annotation(
        Line(points = {{329, -1071}, {329, -929}, {336, -929}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[1], R1.fuelFlowFraction) annotation(
        Line(points = {{336, -865}, {-586, -865}, {-586, -424}, {-504, -424}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[2], R2.fuelFlowFraction) annotation(
        Line(points = {{336, -865}, {-154, -865}, {-154, -622}, {-94, -622}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[2], R3.fuelFlowFraction) annotation(
        Line(points = {{336, -865}, {-156, -865}, {-156, -424}, {-94, -424}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[2], R4.fuelFlowFraction) annotation(
        Line(points = {{336, -865}, {-156, -865}, {-156, -230}, {-110, -230}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[3], R5.fuelFlowFraction) annotation(
        Line(points = {{336, -865}, {242, -865}, {242, -624}, {278, -624}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[3], R6.fuelFlowFraction) annotation(
        Line(points = {{336, -865}, {242, -865}, {242, -424}, {278, -424}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[3], R7.fuelFlowFraction) annotation(
        Line(points = {{336, -865}, {242, -865}, {242, -230}, {264, -230}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[4], R8.fuelFlowFraction) annotation(
        Line(points = {{336, -865}, {336, -546}, {665, -546}, {665, -545}}, color = {255, 120, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut[4], R9.fuelFlowFraction) annotation(
        Line(points = {{336, -865}, {626, -865}, {626, -342}, {660, -342}}, color = {255, 120, 0}, thickness = 1));
      connect(flowFracIn, mpke.fuelFlowFrac) annotation(
        Line(points = {{329, -1071}, {1050, -1071}, {1050, 312}, {185, 312}}, color = {255, 120, 0}, thickness = 1));
      connect(R1.fuelNode2, upperPlenum.temp_In[1]) annotation(
        Line(points = {{-420, -370}, {-414, -370}, {-414, 34}, {445, 34}, {445, 81}, {449, 81}}, color = {204, 0, 0}, thickness = 1));
      connect(R4.fuelNode2, upperPlenum.temp_In[2]) annotation(
        Line(points = {{-27, -176}, {-8, -176}, {-8, 81}, {449, 81}}, color = {204, 0, 0}, thickness = 1));
      connect(R7.fuelNode2, upperPlenum.temp_In[3]) annotation(
        Line(points = {{345, -177}, {345, 36}, {449, 36}, {449, 81}}, color = {204, 0, 0}, thickness = 1));
      connect(R9.fuelNode2, upperPlenum.temp_In[4]) annotation(
        Line(points = {{740, -290}, {740, -76}, {449, -76}, {449, 81}}, color = {204, 0, 0}, thickness = 1));
      connect(flowDistributor.flowFracOut, upperPlenum.flowFraction) annotation(
        Line(points = {{336, -865}, {-594, -865}, {-594, 81}, {423, 81}}, color = {245, 121, 0}, thickness = 1));
      connect(upperPlenum.temp_Out, tempOut) annotation(
        Line(points = {{448, 146}, {448, 235}, {1117, 235}}, color = {204, 0, 0}, thickness = 1));
      connect(powerblock.decayPowerM, upperPlenum.decayHeat) annotation(
        Line(points = {{-159, 367}, {-196, 367}, {-196, 196}, {606, 196}, {606, 81}, {476, 81}}, color = {220, 138, 221}, thickness = 1));
      connect(tempIn, R8.temp_In) annotation(
        Line(points = {{218, -826}, {218, -576}, {664, -576}}, color = {204, 0, 0}, thickness = 1));
      annotation(
        Diagram(coordinateSystem(extent = {{-600, 560}, {1160, -1120}})),
        Icon(coordinateSystem(extent = {{-220, -60}, {-100, -180}}), graphics = {Rectangle(origin = {-160, -120}, lineThickness = 1.5, extent = {{-60, 60}, {60, -60}}), Text(origin = {-157, -119}, extent = {{-53, 31}, {53, -31}}, textString = "Reactor")}));
    end MSRR9R;

    model MSRR1R
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
      
      SMD_MSR_Modelica.Nuclear.mPKE mpke(
        numGroups = numGroups, 
        lambda = lambda, 
        beta = beta, 
        LAMBDA = LAMBDA, 
        n_0 = n_0, 
        addRho0 = addRho0, 
        nomTauLoop = 8.7006, 
        nomTauCore = 34.8025) annotation(
        Placement(transformation(origin = {-205.6, 84.4}, extent = {{-20.4, -20.4}, {13.6, 13.6}})));
      
      MSRR.Components.FuelChannel fuelchannel(
        rho_fuel = rho_fuel, 
        rho_grap = rho_grap, 
        cP_fuel = cP_fuel, 
        cP_grap = cP_grap, 
        Vdot_fuelNom = volDotFuel, 
        kFN1 = 0.4650, 
        kFN2 = 0.4650, 
        kG = 0.07, 
        kHT_FN1 = 0.5, 
        kHT_FN2 = 0.5, 
        TF1_0 = TF1_0, 
        TF2_0 = TF2_0, 
        TG_0 = TG_0, 
        regionFlowFrac = 1, 
        KF = kFuel, 
        Ac = 1.5882, 
        LF1 = 0.7875, 
        LF2 = 0.7875, 
        OutterRegion = EnableRad, 
        ArF1 = 3.5172, 
        ArF2 = 3.5172, 
        e = 0.080000, 
        Tinf = T_inf, 
        vol_FN1 = 0.2, 
        vol_FN2 = 0.2, 
        vol_GN = 1.758, 
        hAnom = 2.4916e+04) annotation(
        Placement(transformation(origin = {-90.3333, -104.278}, extent = {{-67.2222, -67.2222}, {40.3333, 53.7778}})));
      
      SMD_MSR_Modelica.Nuclear.ReactivityFeedback react(
        FuelTempSetPointNode1 = TF1_0, 
        FuelTempSetPointNode2 = TF2_0, 
        GrapTempSetPoint = TG_0, 
        IF1 = 0.5, 
        IF2 = 0.5, 
        IG = 1, 
        a_F = a_F, 
        a_G = a_G) annotation(
        Placement(transformation(origin = {1.69334, 61.8548}, extent = {{-22.2896, 14.8598}, {14.8598, -22.2896}})));
      
      SMD_MSR_Modelica.Nuclear.PowerBlock powerblock(
        P(displayUnit = "W") = 1E6, 
        TotalFuelVol = 0.5) annotation(
        Placement(transformation(origin = {-201.2, 5.6}, extent = {{-12.8, -25.6}, {19.2, 6.4}})));
      
      SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(
        numGroups = 3, 
        DHYG = {2.3751e-03, 8.7763e-05, 1.9596e-06}, 
        DHlamG = {0.09453, 0.004420, 8.6098E-5}) annotation(
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
        Line(points = {{-6, 66}, {-6, 108}, {-206, 108}, {-206, 92}}, color = {78, 154, 6}, thickness = 1));
      connect(flowFracIn, fuelchannel.fuelFlowFraction) annotation(
        Line(points = {{-370, -124}, {-137, -124}}, color = {255, 120, 0}, thickness = 1));
      connect(fuelchannel.fuelNode2, react.fuelNode2) annotation(
        Line(points = {{-64, -76}, {0, -76}, {0, 43}, {5, 43}}, color = {204, 0, 0}, thickness = 1));
      connect(fuelchannel.grapNode, react.grapNode) annotation(
        Line(points = {{-64, -110}, {-6, -110}, {-6, 43}}, color = {204, 0, 0}, thickness = 1));
      connect(fuelchannel.fuelNode1, react.fuelNode1) annotation(
        Line(points = {{-64, -144}, {-12, -144}, {-12, 43}, {-17, 43}}, color = {204, 0, 0}, thickness = 1));
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
    end MSRR1R;
  end Components;
  
  model R1MSRRuhx
  
  parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
  parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
  parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
  parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
  
  parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
  parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
  parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
  parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
  
  parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
  parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
  parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
  parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
  
    MSRR.Components.HeatExchanger heatExchanger(
      vol_P = 0.045379, 
      vol_T = 0.023343, 
      vol_S = 0.055947, 
      rhoP = rhoFuel, 
      rhoT = rhoHXtube, 
      rhoS = rhoCoolant, 
      cP_P = scpFuel, 
      cP_T = scpHXtube, 
      cP_S = scpCoolant, 
      VdotPnom = volDotFuel, 
      VdotSnom = volDotCoolant, 
      hApNom = 7.7966e+04, 
      hAsNom = 1.8150e+04, 
      EnableRad = false, 
      AcShell = 0.1298, 
      AcTube = 0.020276, 
      ArShell = 3.1165, 
      Kp = kFuel, 
      Ks = kCoolant, 
      L_shell = 6.2592, 
      L_tube = 12.5185, 
      e = 0.01, 
      Tinf = 500, 
      TpIn_0 = 580.41, 
      TpOut_0 = 560, 
      TsIn_0 = 500, 
      TsOut_0 = 507.84) annotation(
      Placement(transformation(origin = {-28.2, 32.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
  
    SMD_MSR_Modelica.HeatTransport.UHX uhx(
      Tp_0 = 500, 
      vDot = volDotCoolant, 
      cP = scpCoolant, 
      rho = rhoCoolant, 
      vol = 0.204183209633634)  annotation(
      Placement(transformation(origin = {77.462, -70.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
  
    SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
      vol = 0.2526, 
      volFracNode = 0.1, 
      vDotNom = volDotCoolant, 
      rho = rhoCoolant, 
      cP = scpCoolant, 
      K = kCoolant, 
      Ac = 0.012673, 
      L = 19.932, 
      EnableRad = false, 
      Ar = 7.9559, 
      e = 0.08, 
      Tinf = 644.4, 
      T_0 = 507.84)  annotation(
      Placement(transformation(origin = {-28.8, -81.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
  
    SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
      vol = 0.4419, 
      volFracNode = 0.1, 
      vDotNom = volDotCoolant, 
      rho = rhoCoolant, 
      cP = scpCoolant, 
      K = kCoolant, 
      Ac = 0.012673, 
      L = 34.870, 
      EnableRad = false, 
      Ar = 13.918, 
      e = 0.08, 
      Tinf = 644.4, 
      T_0 = 500) annotation(
      Placement(transformation(origin = {93.2, 7.73333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
  
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
      vol = 0.04, 
      rho = rhoFuel, 
      cP = scpFuel, 
      vDotNom = volDotFuel, 
      DHRS_tK = 10, 
      DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
      DHRS_P_Bleed = 0, 
      DHRS_time = 1000000, 
      K = kFuel, 
      Ac = 0.282743339, 
      L = 0.141471061, 
      Ar = 0.266666667, 
      EnableRad = false, 
      e = 0.08, 
      Tinf = 644.5, 
      T_0 = 580.41)  annotation(
      Placement(transformation(origin = {-267, 32.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
      vol = 4.8738e-03, 
      volFracNode = 0.01, 
      vDotNom = volDotFuel, 
      rho = rhoFuel, 
      cP = scpFuel, 
      K = kFuel, 
      Ac = 2.7321e-03, 
      L = 1.7839, 
      EnableRad = false, 
      Ar = 0.3305, 
      e = 0.08, 
      Tinf = 644.4, 
      T_0 = 580.41)  annotation(
      Placement(transformation(origin = {-135.6, 23.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
      vol = 7.3107e-03, 
      volFracNode = 0.1, 
      vDotNom = volDotFuel, 
      rho = rhoFuel, 
      cP = scpFuel, 
      K = kFuel, 
      Ac = 2.7321e-03, 
      L = 2.6758, 
      EnableRad = false, 
      Ar = 0.4958, 
      e = 0.08, 
      Tinf = 644.4, 
      T_0 = 560) annotation(
      Placement(transformation(origin = {-176.4, -61.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
      vol = 2.4369e-03, 
      volFracNode = 0.1, 
      vDotNom = volDotFuel, 
      rho = rhoFuel, 
      cP = scpFuel, 
      K = kFuel, 
      Ac = 2.7321e-03, 
      L = 0.8919, 
      EnableRad = false, 
      Ar = 0.1653, 
      e = 0.08, 
      Tinf = 644.4, 
      T_0 = 580.41)  annotation(
      Placement(transformation(origin = {-382.4, -32.8}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
    
    MSRR.Components.MSRR1R core1R(
      numGroups = 6, 
      addRho0 = true, 
      EnableRad = false, 
      lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
      beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
      LAMBDA = 2.400E-04, 
      a_F = -6.26E-5,
      a_G = -5.16E-5,
      n_0 = 1, 
      rho_fuel = rhoFuel, 
      rho_grap = rhoGrap, 
      cP_fuel = scpFuel, 
      cP_grap = scpGrap,
      volDotFuel = volDotFuel, T_inf = 570, kFuel = kFuel, TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 570.0000028)  annotation(
      Placement(transformation(origin = {-519.334, -66.0002}, extent = {{-110, -89.9998}, {-49.9999, -29.9999}})));
    
    SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
      Placement(transformation(origin = {-27, -161}, extent = {{-27, -27}, {27, 27}})));
    
    SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
      numRampUp = 1, 
      rampUpK = {1}, 
      rampUpTo = {1}, 
      rampUpTime = {0}, 
      tripK = 50, 
      tripTime = 10000000, 
      freeConvFF = 0)  annotation(
      Placement(transformation(origin = {-554, 158}, extent = {{-24, -32}, {24, 16}})));
    
    SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
      numRampUp = 1, 
      rampUpK = {100}, 
      rampUpTo = {1}, 
      rampUpTime = {0}, 
      tripK = 50, 
      tripTime = 10000000, 
      freeConvFF = 0)  annotation(
      Placement(transformation(origin = {175, 126.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
    
    SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
      Placement(transformation(origin = {43, -209}, extent = {{-27, -27}, {27, 27}})));
    
    SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
      numSteps = 2, 
      stepTime = {0, 4000}, 
      amplitude = {0, 0})  annotation(
      Placement(transformation(origin = {-708.198, 63.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
  
  equation
  
  connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
      Line(points = {{-58, 20}, {-94, 20}, {-94, -58}, {-48, -58}}, color = {204, 0, 0}, thickness = 1));
  connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
      Line(points = {{-18, -58}, {23, -58}, {23, -70}}, color = {204, 0, 0}, thickness = 1));
  connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
      Line(points = {{77, -70}, {140, -70}, {140, 30}, {113, 30}}, color = {204, 0, 0}, thickness = 1));
  connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
      Line(points = {{83, 30}, {60.5, 30}, {60.5, 20}, {52, 20}}, color = {204, 0, 0}, thickness = 1));
  connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
      Line(points = {{-247, 45}, {-205.5, 45}, {-205.5, 55}, {-163, 55}}, color = {204, 0, 0}, thickness = 1));
  connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
      Line(points = {{-121, 55}, {-95, 55}, {-95, 45}, {-58, 45}}, color = {204, 0, 0}, thickness = 1));
  connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
      Line(points = {{52, 45}, {156, 45}, {156, -93}, {-150, -93}}, color = {204, 0, 0}, thickness = 1));
  connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
      Line(points = {{-367, 0}, {-328.5, 0}, {-328.5, 45}, {-300, 45}}, color = {204, 0, 0}, thickness = 1));
  connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
      Line(points = {{-659, -166}, {-409.5, -166}, {-409.5, 0}, {-411, 0}}, color = {204, 0, 0}, thickness = 1));
  connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
      Line(points = {{-191, -93}, {-390, -93}, {-390, -205}, {-699, -205}}, color = {204, 0, 0}, thickness = 1));
  connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
      Line(points = {{-659, -205}, {-448, -205}, {-448, 10}, {-411, 10}}, color = {220, 138, 221}, thickness = 1));
  connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
      Line(points = {{-659, -205}, {-448, -205}, {-448, 70}, {-300, 70}}, color = {220, 138, 221}, thickness = 1));
  connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
      Line(points = {{-659, -205}, {-659, 98}, {-163, 98}, {-163, 65}}, color = {220, 138, 221}, thickness = 1));
  connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
      Line(points = {{-659, -205}, {-659, 114}, {-3, 114}, {-3, 54}}, color = {220, 138, 221}, thickness = 1));
  connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
      Line(points = {{-659, -205}, {-353.5, -205}, {-353.5, -102}, {-150, -102}}, color = {220, 138, 221}, thickness = 1));
  connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
      Line(points = {{-29, -145}, {-74, -145}, {-74, -52}, {-48, -52}}, color = {220, 138, 221}, thickness = 1));
  connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
      Line(points = {{-29, -145}, {148, -145}, {148, 37}, {113, 37}}, color = {220, 138, 221}, thickness = 1));
  connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
      Line(points = {{-554, 126}, {-554, 20}, {-699, 20}, {-699, -166}}, color = {245, 121, 0}, thickness = 1));
  connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
      Line(points = {{-554, 126}, {-468, 126}, {-468, -10}, {-411, -10}}, color = {245, 121, 0}, thickness = 1));
  connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
      Line(points = {{-554, 126}, {-328, 126}, {-328, 20}, {-300, 20}}, color = {245, 121, 0}, thickness = 1));
  connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
      Line(points = {{-554, 126}, {-206, 126}, {-206, 46}, {-163, 46}}, color = {245, 121, 0}, thickness = 1));
  connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
      Line(points = {{-554, 126}, {-554, 92}, {-45, 92}, {-45, 54}}, color = {245, 121, 0}, thickness = 1));
  connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
      Line(points = {{-554, 126}, {-108, 126}, {-108, -83}, {-150, -83}}, color = {245, 121, 0}, thickness = 1));
  connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
      Line(points = {{175, 96}, {174, 96}, {174, 24}, {113, 24}}, color = {245, 121, 0}, thickness = 1));
  connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
      Line(points = {{175, 96}, {175, 98}, {174, 98}, {174, -6}, {40, -6}, {40, 11}}, color = {245, 121, 0}, thickness = 1));
  connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
      Line(points = {{175, 96}, {174, 96}, {174, -43}, {23, -43}}, color = {245, 121, 0}, thickness = 1));
  connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
      Line(points = {{175, 96}, {175, 98}, {176, 98}, {176, -82}, {-48, -82}, {-48, -65}}, color = {245, 121, 0}, thickness = 1));
  connect(externalReact.step, core1R.realIn) annotation(
      Line(points = {{-706, 74}, {-600, 74}, {-600, -106}}, thickness = 1));
  connect(uhxDemand.realOut, uhx.powDemand) annotation(
      Line(points = {{44, -196}, {36, -196}, {36, -98}}, thickness = 1));
  annotation(
      Diagram(coordinateSystem(extent = {{-760, 180}, {220, -280}}), graphics = {Polygon(points = {{108, 32}, {108, 32}})}));
  end R1MSRRuhx;
  
  model R9MSRRuhx
  
  parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
  parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
  parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
  parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
  
  parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
  parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
  parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
  parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
  
  parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
  parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
  parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
  parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
  
    MSRR.Components.HeatExchanger heatExchanger(
      vol_P = 0.045379, 
      vol_T = 0.023343, 
      vol_S = 0.055947, 
      rhoP = rhoFuel, 
      rhoT = rhoHXtube, 
      rhoS = rhoCoolant, 
      cP_P = scpFuel, 
      cP_T = scpHXtube, 
      cP_S = scpCoolant, 
      VdotPnom = volDotFuel, 
      VdotSnom = volDotCoolant, 
      hApNom = 7.7966e+04, 
      hAsNom = 1.8150e+04, 
      EnableRad = false, 
      AcShell = 0.1298, 
      AcTube = 0.020276, 
      ArShell = 3.1165, 
      Kp = kFuel, 
      Ks = kCoolant, 
      L_shell = 6.2592, 
      L_tube = 12.5185, 
      e = 0.01, 
      Tinf = 500, 
      TpIn_0 = 580.41, 
      TpOut_0 = 560, 
      TsIn_0 = 500, 
      TsOut_0 = 507.84) annotation(
      Placement(transformation(origin = {-12.2, 90.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
  
    SMD_MSR_Modelica.HeatTransport.UHX uhx(
      Tp_0 = 500, 
      vDot = volDotCoolant, 
      cP = scpCoolant, 
      rho = rhoCoolant, 
      vol = 0.204183209633634)  annotation(
      Placement(transformation(origin = {103.462, -66.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
  
    SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
      vol = 0.2526, 
      volFracNode = 0.1, 
      vDotNom = volDotCoolant, 
      rho = rhoCoolant, 
      cP = scpCoolant, 
      K = kCoolant, 
      Ac = 0.012673, 
      L = 19.932, 
      EnableRad = false, 
      Ar = 7.9559, 
      e = 0.08, 
      Tinf = 644.4, 
      T_0 = 507.84)  annotation(
      Placement(transformation(origin = {-28.8, -55.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
  
    SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
      vol = 0.4419, 
      volFracNode = 0.1, 
      vDotNom = volDotCoolant, 
      rho = rhoCoolant, 
      cP = scpCoolant, 
      K = kCoolant, 
      Ac = 0.012673, 
      L = 34.870, 
      EnableRad = false, 
      Ar = 13.918, 
      e = 0.08, 
      Tinf = 644.4, 
      T_0 = 500) annotation(
      Placement(transformation(origin = {95.2, 27.7333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
  
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
      vol = 0.04, 
      rho = rhoFuel, 
      cP = scpFuel, 
      vDotNom = volDotFuel, 
      DHRS_tK = 10, 
      DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
      DHRS_P_Bleed = 0, 
      DHRS_time = 1000000, 
      K = kFuel, 
      Ac = 0.282743339, 
      L = 0.141471061, 
      Ar = 0.266666667, 
      EnableRad = false, 
      e = 0.08, 
      Tinf = 644.5, 
      T_0 = 580.41)  annotation(
      Placement(transformation(origin = {-249, 30.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
      vol = 4.8738e-03, 
      volFracNode = 0.01, 
      vDotNom = volDotFuel, 
      rho = rhoFuel, 
      cP = scpFuel, 
      K = kFuel, 
      Ac = 2.7321e-03, 
      L = 1.7839, 
      EnableRad = false, 
      Ar = 0.3305, 
      e = 0.08, 
      Tinf = 644.4, 
      T_0 = 580.41)  annotation(
      Placement(transformation(origin = {-141.6, 11.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
      vol = 7.3107e-03, 
      volFracNode = 0.1, 
      vDotNom = volDotFuel, 
      rho = rhoFuel, 
      cP = scpFuel, 
      K = kFuel, 
      Ac = 2.7321e-03, 
      L = 2.6758, 
      EnableRad = false, 
      Ar = 0.4958, 
      e = 0.08, 
      Tinf = 644.4, 
      T_0 = 560) annotation(
      Placement(transformation(origin = {-158.4, -105.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
      vol = 2.4369e-03, 
      volFracNode = 0.1, 
      vDotNom = volDotFuel, 
      rho = rhoFuel, 
      cP = scpFuel, 
      K = kFuel, 
      Ac = 2.7321e-03, 
      L = 0.8919, 
      EnableRad = false, 
      Ar = 0.1653, 
      e = 0.08, 
      Tinf = 644.4, 
      T_0 = 580.41)  annotation(
      Placement(transformation( origin = {-374, 6}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
    
    SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
      Placement(transformation(origin = {-37, -175}, extent = {{-27, -27}, {27, 27}})));
    
    SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
      numRampUp = 1, 
      rampUpK = {1}, 
      rampUpTo = {1}, 
      rampUpTime = {0}, 
      tripK = 50, 
      tripTime = 1000000, 
      freeConvFF = 0)  annotation(
      Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
    
    SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
      numRampUp = 1, 
      rampUpK = {100}, 
      rampUpTo = {1}, 
      rampUpTime = {0}, 
      tripK = 50, 
      tripTime = 10000000, 
      freeConvFF = 0)  annotation(
      Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
    
    SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
      Placement(transformation(origin = {45, -251}, extent = {{-27, -27}, {27, 27}})));
    
    SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
      numSteps = 2, 
      stepTime = {0, 4000}, 
      amplitude = {0, 0})  annotation(
      Placement(transformation(origin = {-596.198, -38.1983}, extent = {{-23.7983, 23.7983}, {23.7983, -23.7983}}, rotation = -0)));
    
    MSRR.Components.MSRR9R msre9r(
      numGroups = 6, 
      lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
      beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
      LAMBDA = 2.400E-04, 
      n_0 = 1, 
      aF = -6.26E-5, 
      aG = -5.16E-5, 
      addRho0 = true, 
      volF1 = {0.003795391373, 0.012869720861, 0.007038081469, 0.008797601837, 0.021767608195, 0.011889109660, 0.014880331986, 0.059823315476, 0.040594513827}, 
      volF2 = {0.003971456514, 0.008772341956, 0.007038081469, 0.017142787894, 0.014880331986, 0.011889109660, 0.028956494924, 0.034788134316, 0.068118358784}, 
      volG = {0.03488047800, 0.10533760200, 0.08002416000, 0.10244745000, 0.17818736400, 0.13543456200, 0.17330364000, 0.47895127800, 0.46943522400}, 
      volUP = 2*volDotFuel, 
      hA = {513.29, 1430.28, 930.26, 1714.35, 2421.98, 1571.45, 2897.08, 6252.67, 7184.60}, 
      kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, 
      kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, 
      kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, 
      kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, 
      TF1_0 = {566.94, 562.89, 584.55, 609.63, 561.59, 574.98, 590.48, 560.60, 565.19}, 
      TF2_0 = {576.02, 572.69, 597.23, 615.83, 567.64, 582.82, 594.30, 562.63, 566.46}, 
      TG_0 = {570.89, 566.21, 591.18, 613.04, 564.19, 580.21, 593.17, 562.06, 566.66}, 
      Tmix_0 = 580.41, 
      IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, 
      IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, 
      IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, 
      flowFracRegions = {0.061410, 0.138550, 0.234231, 0.565809}, 
      rho_fuel = rhoFuel, 
      cP_fuel = scpFuel, 
      kFuel = kFuel, 
      volDotFuel = volDotFuel,
      rho_grap = rhoGrap, 
      cP_grap = scpGrap, 
      regionTripTime = {200000, 200000, 200000, 200000}, 
      regionCoastDownK = 0.02, 
      freeConvFF = 0.01, 
      EnableRad = false, 
      T_inf = 500, LF1 = {0.7412, 0.3166, 0.1731, 0.2164, 0.3167, 0.1730, 0.2165, 0.4463, 0.3028}, LF2 = {0.7756, 0.2158, 0.1731, 0.4217, 0.2165, 0.1730, 0.4213, 0.2595, 0.5082}, Ac = {0.016189, 0.036524, 0.061748, 0.149159}, ArF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, ArF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, e = 0.08)  annotation(
      Placement(transformation(origin = {-258.333, -35}, extent = {{-179.667, -147}, {-81.6667, -49}})));
  
  equation
    connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
      Line(points = {{-42, 78}, {-94, 78}, {-94, -32}, {-51, -32}}, color = {204, 0, 0}, thickness = 1));
    connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
      Line(points = {{-15, -32}, {28, -32}, {28, -66}, {49, -66}}, color = {204, 0, 0}, thickness = 1));
    connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
      Line(points = {{103, -66}, {140, -66}, {140, 50}, {118, 50}}, color = {204, 0, 0}, thickness = 1));
    connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
      Line(points = {{82, 50}, {57.5, 50}, {57.5, 78}, {68, 78}}, color = {204, 0, 0}, thickness = 1));
    connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
      Line(points = {{-229, 43}, {-174, 43}}, color = {204, 0, 0}, thickness = 1));
    connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
      Line(points = {{-122, 43}, {-76, 43}, {-76, 103}, {-42, 103}}, color = {204, 0, 0}, thickness = 1));
    connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
      Line(points = {{68, 103}, {156, 103}, {156, -137}, {-132, -137}}, color = {204, 0, 0}, thickness = 1));
    connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
      Line(points = {{-359, 39}, {-322.5, 39}, {-322.5, 43}, {-282, 43}}, color = {204, 0, 0}, thickness = 1));
    connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
      Line(points = {{-39, -159}, {-74, -159}, {-74, -23}, {-51, -23}}, color = {220, 138, 221}, thickness = 1));
    connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
      Line(points = {{-39, -159}, {148, -159}, {148, 59}, {118, 59}}, color = {220, 138, 221}, thickness = 1));
    connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
      Line(points = {{-550, 92}, {-468, 92}, {-468, 29}, {-403, 29}}, color = {245, 121, 0}, thickness = 1));
    connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
      Line(points = {{-550, 92}, {-328, 92}, {-328, 18}, {-282, 18}}, color = {245, 121, 0}, thickness = 1));
    connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
      Line(points = {{-550, 92}, {-206, 92}, {-206, 30}, {-174, 30}}, color = {245, 121, 0}, thickness = 1));
    connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
      Line(points = {{-550, 92}, {-299.5, 92}, {-299.5, 112}, {-29, 112}}, color = {245, 121, 0}, thickness = 1));
    connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
      Line(points = {{-550, 92}, {-108, 92}, {-108, -127}, {-132, -127}}, color = {245, 121, 0}, thickness = 1));
    connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
      Line(points = {{175, 92}, {174, 92}, {174, 41}, {118, 41}}, color = {245, 121, 0}, thickness = 1));
    connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
      Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {56, -6}, {56, 69}}, color = {245, 121, 0}, thickness = 1));
    connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
      Line(points = {{175, 92}, {174, 92}, {174, -39}, {49, -39}}, color = {245, 121, 0}, thickness = 1));
    connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
      Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-51, -82}, {-51, -41}}, color = {245, 121, 0}, thickness = 1));
    connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
      Line(points = {{-173, -137}, {-175, -262}, {-552, -262}}, color = {204, 0, 0}, thickness = 1));
    connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
      Line(points = {{-486, -263}, {-319, -263}, {-319, -146}, {-132, -146}}, color = {220, 138, 221}, thickness = 1));
    connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
      Line(points = {{-486, -263}, {-434, -263}, {-434, 49}, {-403, 49}}, color = {220, 138, 221}, thickness = 1));
    connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
      Line(points = {{-486, -263}, {-318, -263}, {-318, 68}, {-282, 68}}, color = {220, 138, 221}, thickness = 1));
    connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
      Line(points = {{-486, -263}, {-200, -263}, {-200, 56}, {-174, 56}}, color = {220, 138, 221}, thickness = 1));
    connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
      Line(points = {{-486, -263}, {-434, -263}, {-434, 112}, {13, 112}}, color = {220, 138, 221}, thickness = 1));
    connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
      Line(points = {{-550, 92}, {-550, -53}, {-552, -53}, {-552, -198}}, color = {245, 121, 0}, thickness = 1));
    connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
      Line(points = {{-487, -198}, {-486.5, -198}, {-486.5, 39}, {-403, 39}}, color = {204, 0, 0}, thickness = 1));
  connect(externalReact.step, msre9r.realIn) annotation(
      Line(points = {{-594, -46}, {-520, -46}, {-520, -198}}, thickness = 1));
  connect(uhxDemand.realOut, uhx.powDemand) annotation(
      Line(points = {{46, -238}, {62, -238}, {62, -94}}, thickness = 1));
    annotation(
      Diagram(coordinateSystem(extent = {{-680, 160}, {220, -320}})));
  end R9MSRRuhx;

  package Transients
   
package R1fullSteps
      model R1MSRR2dol
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
      parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
      parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
      parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
      
      parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
      parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
      
        MSRR.Components.HeatExchanger heatExchanger(
          vol_P = 0.045379, 
          vol_T = 0.023343, 
          vol_S = 0.055947, 
          rhoP = rhoFuel, 
          rhoT = rhoHXtube, 
          rhoS = rhoCoolant, 
          cP_P = scpFuel, 
          cP_T = scpHXtube, 
          cP_S = scpCoolant, 
          VdotPnom = volDotFuel, 
          VdotSnom = volDotCoolant, 
          hApNom = 7.7966e+04, 
          hAsNom = 1.8150e+04, 
          EnableRad = false, 
          AcShell = 0.1298, 
          AcTube = 0.020276, 
          ArShell = 3.1165, 
          Kp = kFuel, 
          Ks = kCoolant, 
          L_shell = 6.2592, 
          L_tube = 12.5185, 
          e = 0.01, 
          Tinf = 500, 
          TpIn_0 = 580.41, 
          TpOut_0 = 560, 
          TsIn_0 = 500, 
          TsOut_0 = 507.84) annotation(
          Placement(transformation(origin = {-28.2, 32.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
      
        SMD_MSR_Modelica.HeatTransport.UHX uhx(
          Tp_0 = 500, 
          vDot = volDotCoolant, 
          cP = scpCoolant, 
          rho = rhoCoolant, 
          vol = 0.204183209633634)  annotation(
          Placement(transformation(origin = {77.462, -70.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
          vol = 0.2526, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 19.932, 
          EnableRad = false, 
          Ar = 7.9559, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 507.84)  annotation(
          Placement(transformation(origin = {-28.8, -81.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
          vol = 0.4419, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 34.870, 
          EnableRad = false, 
          Ar = 13.918, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 500) annotation(
          Placement(transformation(origin = {93.2, 7.73333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
          vol = 0.04, 
          rho = rhoFuel, 
          cP = scpFuel, 
          vDotNom = volDotFuel, 
          DHRS_tK = 10, 
          DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
          DHRS_P_Bleed = 0, 
          DHRS_time = 1000000, 
          K = kFuel, 
          Ac = 0.282743339, 
          L = 0.141471061, 
          Ar = 0.266666667, 
          EnableRad = false, 
          e = 0.08, 
          Tinf = 644.5, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-267, 32.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
          vol = 4.8738e-03, 
          volFracNode = 0.01, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 1.7839, 
          EnableRad = false, 
          Ar = 0.3305, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-135.6, 23.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
          vol = 7.3107e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 2.6758, 
          EnableRad = false, 
          Ar = 0.4958, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 560) annotation(
          Placement(transformation(origin = {-176.4, -61.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
          vol = 2.4369e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 0.8919, 
          EnableRad = false, 
          Ar = 0.1653, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-382.4, -32.8}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        
        MSRR.Components.MSRR1R core1R(
          numGroups = 6, 
          addRho0 = true, 
          EnableRad = false, 
          lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
          beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
          LAMBDA = 2.400E-04, 
          a_F = -6.26E-5,
          a_G = -5.16E-5,
          n_0 = 1, 
          rho_fuel = rhoFuel, 
          rho_grap = rhoGrap, 
          cP_fuel = scpFuel, 
          cP_grap = scpGrap,
          volDotFuel = volDotFuel, T_inf = 570, kFuel = kFuel, TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 570.0000028)  annotation(
          Placement(transformation(origin = {-519.334, -66.0002}, extent = {{-110, -89.9998}, {-49.9999, -29.9999}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
          Placement(transformation(origin = {-27, -161}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
          numRampUp = 1, 
          rampUpK = {1}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {-554, 158}, extent = {{-24, -32}, {24, 16}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
          numRampUp = 1, 
          rampUpK = {100}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {175, 126.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
          Placement(transformation(origin = {43, -209}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
          numSteps = 2, 
          stepTime = {0, 4000}, 
          amplitude = {0, 1178.138})  annotation(
          Placement(transformation(origin = {-708.198, 63.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
      
      equation
      
      connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-58, 20}, {-94, 20}, {-94, -58}, {-48, -58}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-18, -58}, {23, -58}, {23, -70}}, color = {204, 0, 0}, thickness = 1));
      connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{77, -70}, {140, -70}, {140, 30}, {113, 30}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{83, 30}, {60.5, 30}, {60.5, 20}, {52, 20}}, color = {204, 0, 0}, thickness = 1));
      connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-247, 45}, {-205.5, 45}, {-205.5, 55}, {-163, 55}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-121, 55}, {-95, 55}, {-95, 45}, {-58, 45}}, color = {204, 0, 0}, thickness = 1));
      connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{52, 45}, {156, 45}, {156, -93}, {-150, -93}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-367, 0}, {-328.5, 0}, {-328.5, 45}, {-300, 45}}, color = {204, 0, 0}, thickness = 1));
      connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-659, -166}, {-409.5, -166}, {-409.5, 0}, {-411, 0}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
          Line(points = {{-191, -93}, {-390, -93}, {-390, -205}, {-699, -205}}, color = {204, 0, 0}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-448, -205}, {-448, 10}, {-411, 10}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-659, -205}, {-448, -205}, {-448, 70}, {-300, 70}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-659, 98}, {-163, 98}, {-163, 65}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-659, -205}, {-659, 114}, {-3, 114}, {-3, 54}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-353.5, -205}, {-353.5, -102}, {-150, -102}}, color = {220, 138, 221}, thickness = 1));
      connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-29, -145}, {-74, -145}, {-74, -52}, {-48, -52}}, color = {220, 138, 221}, thickness = 1));
      connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-29, -145}, {148, -145}, {148, 37}, {113, 37}}, color = {220, 138, 221}, thickness = 1));
      connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
          Line(points = {{-554, 126}, {-554, 20}, {-699, 20}, {-699, -166}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-554, 126}, {-468, 126}, {-468, -10}, {-411, -10}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-554, 126}, {-328, 126}, {-328, 20}, {-300, 20}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-554, 126}, {-206, 126}, {-206, 46}, {-163, 46}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-554, 126}, {-554, 92}, {-45, 92}, {-45, 54}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-554, 126}, {-108, 126}, {-108, -83}, {-150, -83}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 96}, {174, 96}, {174, 24}, {113, 24}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 96}, {175, 98}, {174, 98}, {174, -6}, {40, -6}, {40, 11}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 96}, {174, 96}, {174, -43}, {23, -43}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 96}, {175, 98}, {176, 98}, {176, -82}, {-48, -82}, {-48, -65}}, color = {245, 121, 0}, thickness = 1));
      connect(externalReact.step, core1R.realIn) annotation(
          Line(points = {{-706, 74}, {-600, 74}, {-600, -106}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{44, -196}, {36, -196}, {36, -98}}, thickness = 1));
      annotation(
          Diagram(coordinateSystem(extent = {{-760, 180}, {220, -280}}), graphics = {Polygon(points = {{108, 32}, {108, 32}})}));
      end R1MSRR2dol;
      
      model R1MSRR1dol
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
      parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
      parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
      parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
      
      parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
      parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
      
        MSRR.Components.HeatExchanger heatExchanger(
          vol_P = 0.045379, 
          vol_T = 0.023343, 
          vol_S = 0.055947, 
          rhoP = rhoFuel, 
          rhoT = rhoHXtube, 
          rhoS = rhoCoolant, 
          cP_P = scpFuel, 
          cP_T = scpHXtube, 
          cP_S = scpCoolant, 
          VdotPnom = volDotFuel, 
          VdotSnom = volDotCoolant, 
          hApNom = 7.7966e+04, 
          hAsNom = 1.8150e+04, 
          EnableRad = false, 
          AcShell = 0.1298, 
          AcTube = 0.020276, 
          ArShell = 3.1165, 
          Kp = kFuel, 
          Ks = kCoolant, 
          L_shell = 6.2592, 
          L_tube = 12.5185, 
          e = 0.01, 
          Tinf = 500, 
          TpIn_0 = 580.41, 
          TpOut_0 = 560, 
          TsIn_0 = 500, 
          TsOut_0 = 507.84) annotation(
          Placement(transformation(origin = {-28.2, 32.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
      
        SMD_MSR_Modelica.HeatTransport.UHX uhx(
          Tp_0 = 500, 
          vDot = volDotCoolant, 
          cP = scpCoolant, 
          rho = rhoCoolant, 
          vol = 0.204183209633634)  annotation(
          Placement(transformation(origin = {77.462, -70.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
          vol = 0.2526, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 19.932, 
          EnableRad = false, 
          Ar = 7.9559, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 507.84)  annotation(
          Placement(transformation(origin = {-28.8, -81.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
          vol = 0.4419, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 34.870, 
          EnableRad = false, 
          Ar = 13.918, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 500) annotation(
          Placement(transformation(origin = {93.2, 7.73333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
          vol = 0.04, 
          rho = rhoFuel, 
          cP = scpFuel, 
          vDotNom = volDotFuel, 
          DHRS_tK = 10, 
          DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
          DHRS_P_Bleed = 0, 
          DHRS_time = 1000000, 
          K = kFuel, 
          Ac = 0.282743339, 
          L = 0.141471061, 
          Ar = 0.266666667, 
          EnableRad = false, 
          e = 0.08, 
          Tinf = 644.5, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-267, 32.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
          vol = 4.8738e-03, 
          volFracNode = 0.01, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 1.7839, 
          EnableRad = false, 
          Ar = 0.3305, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-135.6, 23.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
          vol = 7.3107e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 2.6758, 
          EnableRad = false, 
          Ar = 0.4958, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 560) annotation(
          Placement(transformation(origin = {-176.4, -61.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
          vol = 2.4369e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 0.8919, 
          EnableRad = false, 
          Ar = 0.1653, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-382.4, -32.8}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        
        MSRR.Components.MSRR1R core1R(
          numGroups = 6, 
          addRho0 = true, 
          EnableRad = false, 
          lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
          beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
          LAMBDA = 2.400E-04, 
          a_F = -6.26E-5,
          a_G = -5.16E-5,
          n_0 = 1, 
          rho_fuel = rhoFuel, 
          rho_grap = rhoGrap, 
          cP_fuel = scpFuel, 
          cP_grap = scpGrap,
          volDotFuel = volDotFuel, T_inf = 570, kFuel = kFuel, TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 570.0000028)  annotation(
          Placement(transformation(origin = {-519.334, -66.0002}, extent = {{-110, -89.9998}, {-49.9999, -29.9999}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
          Placement(transformation(origin = {-27, -161}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
          numRampUp = 1, 
          rampUpK = {1}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {-554, 158}, extent = {{-24, -32}, {24, 16}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
          numRampUp = 1, 
          rampUpK = {100}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {175, 126.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
          Placement(transformation(origin = {43, -209}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
          numSteps = 2, 
          stepTime = {0, 4000}, 
          amplitude = {0, 589.069})  annotation(
          Placement(transformation(origin = {-708.198, 63.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
      
      equation
      
      connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-58, 20}, {-94, 20}, {-94, -58}, {-48, -58}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-18, -58}, {23, -58}, {23, -70}}, color = {204, 0, 0}, thickness = 1));
      connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{77, -70}, {140, -70}, {140, 30}, {113, 30}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{83, 30}, {60.5, 30}, {60.5, 20}, {52, 20}}, color = {204, 0, 0}, thickness = 1));
      connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-247, 45}, {-205.5, 45}, {-205.5, 55}, {-163, 55}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-121, 55}, {-95, 55}, {-95, 45}, {-58, 45}}, color = {204, 0, 0}, thickness = 1));
      connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{52, 45}, {156, 45}, {156, -93}, {-150, -93}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-367, 0}, {-328.5, 0}, {-328.5, 45}, {-300, 45}}, color = {204, 0, 0}, thickness = 1));
      connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-659, -166}, {-409.5, -166}, {-409.5, 0}, {-411, 0}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
          Line(points = {{-191, -93}, {-390, -93}, {-390, -205}, {-699, -205}}, color = {204, 0, 0}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-448, -205}, {-448, 10}, {-411, 10}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-659, -205}, {-448, -205}, {-448, 70}, {-300, 70}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-659, 98}, {-163, 98}, {-163, 65}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-659, -205}, {-659, 114}, {-3, 114}, {-3, 54}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-353.5, -205}, {-353.5, -102}, {-150, -102}}, color = {220, 138, 221}, thickness = 1));
      connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-29, -145}, {-74, -145}, {-74, -52}, {-48, -52}}, color = {220, 138, 221}, thickness = 1));
      connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-29, -145}, {148, -145}, {148, 37}, {113, 37}}, color = {220, 138, 221}, thickness = 1));
      connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
          Line(points = {{-554, 126}, {-554, 20}, {-699, 20}, {-699, -166}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-554, 126}, {-468, 126}, {-468, -10}, {-411, -10}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-554, 126}, {-328, 126}, {-328, 20}, {-300, 20}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-554, 126}, {-206, 126}, {-206, 46}, {-163, 46}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-554, 126}, {-554, 92}, {-45, 92}, {-45, 54}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-554, 126}, {-108, 126}, {-108, -83}, {-150, -83}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 96}, {174, 96}, {174, 24}, {113, 24}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 96}, {175, 98}, {174, 98}, {174, -6}, {40, -6}, {40, 11}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 96}, {174, 96}, {174, -43}, {23, -43}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 96}, {175, 98}, {176, 98}, {176, -82}, {-48, -82}, {-48, -65}}, color = {245, 121, 0}, thickness = 1));
      connect(externalReact.step, core1R.realIn) annotation(
          Line(points = {{-706, 74}, {-600, 74}, {-600, -106}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{44, -196}, {36, -196}, {36, -98}}, thickness = 1));
      annotation(
          Diagram(coordinateSystem(extent = {{-760, 180}, {220, -280}}), graphics = {Polygon(points = {{108, 32}, {108, 32}})}));
      end R1MSRR1dol;
      
      model R1MSRRhalfDol
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
      parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
      parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
      parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
      
      parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
      parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
      
        MSRR.Components.HeatExchanger heatExchanger(
          vol_P = 0.045379, 
          vol_T = 0.023343, 
          vol_S = 0.055947, 
          rhoP = rhoFuel, 
          rhoT = rhoHXtube, 
          rhoS = rhoCoolant, 
          cP_P = scpFuel, 
          cP_T = scpHXtube, 
          cP_S = scpCoolant, 
          VdotPnom = volDotFuel, 
          VdotSnom = volDotCoolant, 
          hApNom = 7.7966e+04, 
          hAsNom = 1.8150e+04, 
          EnableRad = false, 
          AcShell = 0.1298, 
          AcTube = 0.020276, 
          ArShell = 3.1165, 
          Kp = kFuel, 
          Ks = kCoolant, 
          L_shell = 6.2592, 
          L_tube = 12.5185, 
          e = 0.01, 
          Tinf = 500, 
          TpIn_0 = 580.41, 
          TpOut_0 = 560, 
          TsIn_0 = 500, 
          TsOut_0 = 507.84) annotation(
          Placement(transformation(origin = {-28.2, 32.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
      
        SMD_MSR_Modelica.HeatTransport.UHX uhx(
          Tp_0 = 500, 
          vDot = volDotCoolant, 
          cP = scpCoolant, 
          rho = rhoCoolant, 
          vol = 0.204183209633634)  annotation(
          Placement(transformation(origin = {77.462, -70.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
          vol = 0.2526, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 19.932, 
          EnableRad = false, 
          Ar = 7.9559, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 507.84)  annotation(
          Placement(transformation(origin = {-28.8, -81.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
          vol = 0.4419, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 34.870, 
          EnableRad = false, 
          Ar = 13.918, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 500) annotation(
          Placement(transformation(origin = {93.2, 7.73333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
          vol = 0.04, 
          rho = rhoFuel, 
          cP = scpFuel, 
          vDotNom = volDotFuel, 
          DHRS_tK = 10, 
          DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
          DHRS_P_Bleed = 0, 
          DHRS_time = 1000000, 
          K = kFuel, 
          Ac = 0.282743339, 
          L = 0.141471061, 
          Ar = 0.266666667, 
          EnableRad = false, 
          e = 0.08, 
          Tinf = 644.5, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-267, 32.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
          vol = 4.8738e-03, 
          volFracNode = 0.01, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 1.7839, 
          EnableRad = false, 
          Ar = 0.3305, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-135.6, 23.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
          vol = 7.3107e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 2.6758, 
          EnableRad = false, 
          Ar = 0.4958, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 560) annotation(
          Placement(transformation(origin = {-176.4, -61.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
          vol = 2.4369e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 0.8919, 
          EnableRad = false, 
          Ar = 0.1653, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-382.4, -32.8}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        
        MSRR.Components.MSRR1R core1R(
          numGroups = 6, 
          addRho0 = true, 
          EnableRad = false, 
          lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
          beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
          LAMBDA = 2.400E-04, 
          a_F = -6.26E-5,
          a_G = -5.16E-5,
          n_0 = 1, 
          rho_fuel = rhoFuel, 
          rho_grap = rhoGrap, 
          cP_fuel = scpFuel, 
          cP_grap = scpGrap,
          volDotFuel = volDotFuel, T_inf = 570, kFuel = kFuel, TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 570.0000028)  annotation(
          Placement(transformation(origin = {-519.334, -66.0002}, extent = {{-110, -89.9998}, {-49.9999, -29.9999}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
          Placement(transformation(origin = {-27, -161}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
          numRampUp = 1, 
          rampUpK = {1}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {-554, 158}, extent = {{-24, -32}, {24, 16}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
          numRampUp = 1, 
          rampUpK = {100}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {175, 126.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
          Placement(transformation(origin = {43, -209}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
          numSteps = 2, 
          stepTime = {0, 4000}, 
          amplitude = {0, 294.535})  annotation(
          Placement(transformation(origin = {-708.198, 63.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
      
      equation
      
      connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-58, 20}, {-94, 20}, {-94, -58}, {-48, -58}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-18, -58}, {23, -58}, {23, -70}}, color = {204, 0, 0}, thickness = 1));
      connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{77, -70}, {140, -70}, {140, 30}, {113, 30}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{83, 30}, {60.5, 30}, {60.5, 20}, {52, 20}}, color = {204, 0, 0}, thickness = 1));
      connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-247, 45}, {-205.5, 45}, {-205.5, 55}, {-163, 55}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-121, 55}, {-95, 55}, {-95, 45}, {-58, 45}}, color = {204, 0, 0}, thickness = 1));
      connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{52, 45}, {156, 45}, {156, -93}, {-150, -93}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-367, 0}, {-328.5, 0}, {-328.5, 45}, {-300, 45}}, color = {204, 0, 0}, thickness = 1));
      connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-659, -166}, {-409.5, -166}, {-409.5, 0}, {-411, 0}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
          Line(points = {{-191, -93}, {-390, -93}, {-390, -205}, {-699, -205}}, color = {204, 0, 0}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-448, -205}, {-448, 10}, {-411, 10}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-659, -205}, {-448, -205}, {-448, 70}, {-300, 70}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-659, 98}, {-163, 98}, {-163, 65}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-659, -205}, {-659, 114}, {-3, 114}, {-3, 54}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-353.5, -205}, {-353.5, -102}, {-150, -102}}, color = {220, 138, 221}, thickness = 1));
      connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-29, -145}, {-74, -145}, {-74, -52}, {-48, -52}}, color = {220, 138, 221}, thickness = 1));
      connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-29, -145}, {148, -145}, {148, 37}, {113, 37}}, color = {220, 138, 221}, thickness = 1));
      connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
          Line(points = {{-554, 126}, {-554, 20}, {-699, 20}, {-699, -166}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-554, 126}, {-468, 126}, {-468, -10}, {-411, -10}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-554, 126}, {-328, 126}, {-328, 20}, {-300, 20}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-554, 126}, {-206, 126}, {-206, 46}, {-163, 46}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-554, 126}, {-554, 92}, {-45, 92}, {-45, 54}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-554, 126}, {-108, 126}, {-108, -83}, {-150, -83}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 96}, {174, 96}, {174, 24}, {113, 24}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 96}, {175, 98}, {174, 98}, {174, -6}, {40, -6}, {40, 11}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 96}, {174, 96}, {174, -43}, {23, -43}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 96}, {175, 98}, {176, 98}, {176, -82}, {-48, -82}, {-48, -65}}, color = {245, 121, 0}, thickness = 1));
      connect(externalReact.step, core1R.realIn) annotation(
          Line(points = {{-706, 74}, {-600, 74}, {-600, -106}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{44, -196}, {36, -196}, {36, -98}}, thickness = 1));
      annotation(
          Diagram(coordinateSystem(extent = {{-760, 180}, {220, -280}}), graphics = {Polygon(points = {{108, 32}, {108, 32}})}));
      end R1MSRRhalfDol;
      
      model R1MSRRpOneDol
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
      parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
      parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
      parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
      
      parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
      parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
      
        MSRR.Components.HeatExchanger heatExchanger(
          vol_P = 0.045379, 
          vol_T = 0.023343, 
          vol_S = 0.055947, 
          rhoP = rhoFuel, 
          rhoT = rhoHXtube, 
          rhoS = rhoCoolant, 
          cP_P = scpFuel, 
          cP_T = scpHXtube, 
          cP_S = scpCoolant, 
          VdotPnom = volDotFuel, 
          VdotSnom = volDotCoolant, 
          hApNom = 7.7966e+04, 
          hAsNom = 1.8150e+04, 
          EnableRad = false, 
          AcShell = 0.1298, 
          AcTube = 0.020276, 
          ArShell = 3.1165, 
          Kp = kFuel, 
          Ks = kCoolant, 
          L_shell = 6.2592, 
          L_tube = 12.5185, 
          e = 0.01, 
          Tinf = 500, 
          TpIn_0 = 580.41, 
          TpOut_0 = 560, 
          TsIn_0 = 500, 
          TsOut_0 = 507.84) annotation(
          Placement(transformation(origin = {-28.2, 32.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
      
        SMD_MSR_Modelica.HeatTransport.UHX uhx(
          Tp_0 = 500, 
          vDot = volDotCoolant, 
          cP = scpCoolant, 
          rho = rhoCoolant, 
          vol = 0.204183209633634)  annotation(
          Placement(transformation(origin = {77.462, -70.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
          vol = 0.2526, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 19.932, 
          EnableRad = false, 
          Ar = 7.9559, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 507.84)  annotation(
          Placement(transformation(origin = {-28.8, -81.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
          vol = 0.4419, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 34.870, 
          EnableRad = false, 
          Ar = 13.918, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 500) annotation(
          Placement(transformation(origin = {93.2, 7.73333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
          vol = 0.04, 
          rho = rhoFuel, 
          cP = scpFuel, 
          vDotNom = volDotFuel, 
          DHRS_tK = 10, 
          DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
          DHRS_P_Bleed = 0, 
          DHRS_time = 1000000, 
          K = kFuel, 
          Ac = 0.282743339, 
          L = 0.141471061, 
          Ar = 0.266666667, 
          EnableRad = false, 
          e = 0.08, 
          Tinf = 644.5, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-267, 32.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
          vol = 4.8738e-03, 
          volFracNode = 0.01, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 1.7839, 
          EnableRad = false, 
          Ar = 0.3305, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-135.6, 23.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
          vol = 7.3107e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 2.6758, 
          EnableRad = false, 
          Ar = 0.4958, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 560) annotation(
          Placement(transformation(origin = {-176.4, -61.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
          vol = 2.4369e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 0.8919, 
          EnableRad = false, 
          Ar = 0.1653, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-382.4, -32.8}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        
        MSRR.Components.MSRR1R core1R(
          numGroups = 6, 
          addRho0 = true, 
          EnableRad = false, 
          lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
          beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
          LAMBDA = 2.400E-04, 
          a_F = -6.26E-5,
          a_G = -5.16E-5,
          n_0 = 1, 
          rho_fuel = rhoFuel, 
          rho_grap = rhoGrap, 
          cP_fuel = scpFuel, 
          cP_grap = scpGrap,
          volDotFuel = volDotFuel, T_inf = 570, kFuel = kFuel, TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 570.0000028)  annotation(
          Placement(transformation(origin = {-519.334, -66.0002}, extent = {{-110, -89.9998}, {-49.9999, -29.9999}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
          Placement(transformation(origin = {-27, -161}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
          numRampUp = 1, 
          rampUpK = {1}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {-554, 158}, extent = {{-24, -32}, {24, 16}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
          numRampUp = 1, 
          rampUpK = {100}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {175, 126.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
          Placement(transformation(origin = {43, -209}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
          numSteps = 2, 
          stepTime = {0, 4000}, 
          amplitude = {0, 58.907})  annotation(
          Placement(transformation(origin = {-708.198, 63.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
      
      equation
      
      connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-58, 20}, {-94, 20}, {-94, -58}, {-48, -58}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-18, -58}, {23, -58}, {23, -70}}, color = {204, 0, 0}, thickness = 1));
      connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{77, -70}, {140, -70}, {140, 30}, {113, 30}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{83, 30}, {60.5, 30}, {60.5, 20}, {52, 20}}, color = {204, 0, 0}, thickness = 1));
      connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-247, 45}, {-205.5, 45}, {-205.5, 55}, {-163, 55}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-121, 55}, {-95, 55}, {-95, 45}, {-58, 45}}, color = {204, 0, 0}, thickness = 1));
      connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{52, 45}, {156, 45}, {156, -93}, {-150, -93}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-367, 0}, {-328.5, 0}, {-328.5, 45}, {-300, 45}}, color = {204, 0, 0}, thickness = 1));
      connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-659, -166}, {-409.5, -166}, {-409.5, 0}, {-411, 0}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
          Line(points = {{-191, -93}, {-390, -93}, {-390, -205}, {-699, -205}}, color = {204, 0, 0}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-448, -205}, {-448, 10}, {-411, 10}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-659, -205}, {-448, -205}, {-448, 70}, {-300, 70}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-659, 98}, {-163, 98}, {-163, 65}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-659, -205}, {-659, 114}, {-3, 114}, {-3, 54}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-353.5, -205}, {-353.5, -102}, {-150, -102}}, color = {220, 138, 221}, thickness = 1));
      connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-29, -145}, {-74, -145}, {-74, -52}, {-48, -52}}, color = {220, 138, 221}, thickness = 1));
      connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-29, -145}, {148, -145}, {148, 37}, {113, 37}}, color = {220, 138, 221}, thickness = 1));
      connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
          Line(points = {{-554, 126}, {-554, 20}, {-699, 20}, {-699, -166}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-554, 126}, {-468, 126}, {-468, -10}, {-411, -10}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-554, 126}, {-328, 126}, {-328, 20}, {-300, 20}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-554, 126}, {-206, 126}, {-206, 46}, {-163, 46}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-554, 126}, {-554, 92}, {-45, 92}, {-45, 54}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-554, 126}, {-108, 126}, {-108, -83}, {-150, -83}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 96}, {174, 96}, {174, 24}, {113, 24}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 96}, {175, 98}, {174, 98}, {174, -6}, {40, -6}, {40, 11}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 96}, {174, 96}, {174, -43}, {23, -43}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 96}, {175, 98}, {176, 98}, {176, -82}, {-48, -82}, {-48, -65}}, color = {245, 121, 0}, thickness = 1));
      connect(externalReact.step, core1R.realIn) annotation(
          Line(points = {{-706, 74}, {-600, 74}, {-600, -106}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{44, -196}, {36, -196}, {36, -98}}, thickness = 1));
      annotation(
          Diagram(coordinateSystem(extent = {{-760, 180}, {220, -280}}), graphics = {Polygon(points = {{108, 32}, {108, 32}})}));
      end R1MSRRpOneDol;
      
      model R1MSRR100pcm
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
      parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
      parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
      parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
      
      parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
      parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
      
        MSRR.Components.HeatExchanger heatExchanger(
          vol_P = 0.045379, 
          vol_T = 0.023343, 
          vol_S = 0.055947, 
          rhoP = rhoFuel, 
          rhoT = rhoHXtube, 
          rhoS = rhoCoolant, 
          cP_P = scpFuel, 
          cP_T = scpHXtube, 
          cP_S = scpCoolant, 
          VdotPnom = volDotFuel, 
          VdotSnom = volDotCoolant, 
          hApNom = 7.7966e+04, 
          hAsNom = 1.8150e+04, 
          EnableRad = false, 
          AcShell = 0.1298, 
          AcTube = 0.020276, 
          ArShell = 3.1165, 
          Kp = kFuel, 
          Ks = kCoolant, 
          L_shell = 6.2592, 
          L_tube = 12.5185, 
          e = 0.01, 
          Tinf = 500, 
          TpIn_0 = 580.41, 
          TpOut_0 = 560, 
          TsIn_0 = 500, 
          TsOut_0 = 507.84) annotation(
          Placement(transformation(origin = {-28.2, 32.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
      
        SMD_MSR_Modelica.HeatTransport.UHX uhx(
          Tp_0 = 500, 
          vDot = volDotCoolant, 
          cP = scpCoolant, 
          rho = rhoCoolant, 
          vol = 0.204183209633634)  annotation(
          Placement(transformation(origin = {77.462, -70.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
          vol = 0.2526, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 19.932, 
          EnableRad = false, 
          Ar = 7.9559, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 507.84)  annotation(
          Placement(transformation(origin = {-28.8, -81.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
          vol = 0.4419, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 34.870, 
          EnableRad = false, 
          Ar = 13.918, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 500) annotation(
          Placement(transformation(origin = {93.2, 7.73333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
          vol = 0.04, 
          rho = rhoFuel, 
          cP = scpFuel, 
          vDotNom = volDotFuel, 
          DHRS_tK = 10, 
          DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
          DHRS_P_Bleed = 0, 
          DHRS_time = 1000000, 
          K = kFuel, 
          Ac = 0.282743339, 
          L = 0.141471061, 
          Ar = 0.266666667, 
          EnableRad = false, 
          e = 0.08, 
          Tinf = 644.5, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-267, 32.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
          vol = 4.8738e-03, 
          volFracNode = 0.01, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 1.7839, 
          EnableRad = false, 
          Ar = 0.3305, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-135.6, 23.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
          vol = 7.3107e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 2.6758, 
          EnableRad = false, 
          Ar = 0.4958, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 560) annotation(
          Placement(transformation(origin = {-176.4, -61.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
          vol = 2.4369e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 0.8919, 
          EnableRad = false, 
          Ar = 0.1653, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-382.4, -32.8}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        
        MSRR.Components.MSRR1R core1R(
          numGroups = 6, 
          addRho0 = true, 
          EnableRad = false, 
          lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
          beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
          LAMBDA = 2.400E-04, 
          a_F = -6.26E-5,
          a_G = -5.16E-5,
          n_0 = 1, 
          rho_fuel = rhoFuel, 
          rho_grap = rhoGrap, 
          cP_fuel = scpFuel, 
          cP_grap = scpGrap,
          volDotFuel = volDotFuel, T_inf = 570, kFuel = kFuel, TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 570.0000028)  annotation(
          Placement(transformation(origin = {-519.334, -66.0002}, extent = {{-110, -89.9998}, {-49.9999, -29.9999}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
          Placement(transformation(origin = {-27, -161}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
          numRampUp = 1, 
          rampUpK = {1}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {-554, 158}, extent = {{-24, -32}, {24, 16}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
          numRampUp = 1, 
          rampUpK = {100}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {175, 126.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
          Placement(transformation(origin = {43, -209}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
          numSteps = 2, 
          stepTime = {0, 4000}, 
          amplitude = {0, 100})  annotation(
          Placement(transformation(origin = {-708.198, 63.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
      
      equation
      
      connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-58, 20}, {-94, 20}, {-94, -58}, {-48, -58}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-18, -58}, {23, -58}, {23, -70}}, color = {204, 0, 0}, thickness = 1));
      connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{77, -70}, {140, -70}, {140, 30}, {113, 30}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{83, 30}, {60.5, 30}, {60.5, 20}, {52, 20}}, color = {204, 0, 0}, thickness = 1));
      connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-247, 45}, {-205.5, 45}, {-205.5, 55}, {-163, 55}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-121, 55}, {-95, 55}, {-95, 45}, {-58, 45}}, color = {204, 0, 0}, thickness = 1));
      connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{52, 45}, {156, 45}, {156, -93}, {-150, -93}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-367, 0}, {-328.5, 0}, {-328.5, 45}, {-300, 45}}, color = {204, 0, 0}, thickness = 1));
      connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-659, -166}, {-409.5, -166}, {-409.5, 0}, {-411, 0}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
          Line(points = {{-191, -93}, {-390, -93}, {-390, -205}, {-699, -205}}, color = {204, 0, 0}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-448, -205}, {-448, 10}, {-411, 10}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-659, -205}, {-448, -205}, {-448, 70}, {-300, 70}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-659, 98}, {-163, 98}, {-163, 65}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-659, -205}, {-659, 114}, {-3, 114}, {-3, 54}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-353.5, -205}, {-353.5, -102}, {-150, -102}}, color = {220, 138, 221}, thickness = 1));
      connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-29, -145}, {-74, -145}, {-74, -52}, {-48, -52}}, color = {220, 138, 221}, thickness = 1));
      connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-29, -145}, {148, -145}, {148, 37}, {113, 37}}, color = {220, 138, 221}, thickness = 1));
      connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
          Line(points = {{-554, 126}, {-554, 20}, {-699, 20}, {-699, -166}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-554, 126}, {-468, 126}, {-468, -10}, {-411, -10}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-554, 126}, {-328, 126}, {-328, 20}, {-300, 20}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-554, 126}, {-206, 126}, {-206, 46}, {-163, 46}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-554, 126}, {-554, 92}, {-45, 92}, {-45, 54}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-554, 126}, {-108, 126}, {-108, -83}, {-150, -83}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 96}, {174, 96}, {174, 24}, {113, 24}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 96}, {175, 98}, {174, 98}, {174, -6}, {40, -6}, {40, 11}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 96}, {174, 96}, {174, -43}, {23, -43}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 96}, {175, 98}, {176, 98}, {176, -82}, {-48, -82}, {-48, -65}}, color = {245, 121, 0}, thickness = 1));
      connect(externalReact.step, core1R.realIn) annotation(
          Line(points = {{-706, 74}, {-600, 74}, {-600, -106}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{44, -196}, {36, -196}, {36, -98}}, thickness = 1));
      annotation(
          Diagram(coordinateSystem(extent = {{-760, 180}, {220, -280}}), graphics = {Polygon(points = {{108, 32}, {108, 32}})}));
      end R1MSRR100pcm;
      
      model R1MSRR10pcm
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
      parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
      parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
      parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
      
      parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
      parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
      
        MSRR.Components.HeatExchanger heatExchanger(
          vol_P = 0.045379, 
          vol_T = 0.023343, 
          vol_S = 0.055947, 
          rhoP = rhoFuel, 
          rhoT = rhoHXtube, 
          rhoS = rhoCoolant, 
          cP_P = scpFuel, 
          cP_T = scpHXtube, 
          cP_S = scpCoolant, 
          VdotPnom = volDotFuel, 
          VdotSnom = volDotCoolant, 
          hApNom = 7.7966e+04, 
          hAsNom = 1.8150e+04, 
          EnableRad = false, 
          AcShell = 0.1298, 
          AcTube = 0.020276, 
          ArShell = 3.1165, 
          Kp = kFuel, 
          Ks = kCoolant, 
          L_shell = 6.2592, 
          L_tube = 12.5185, 
          e = 0.01, 
          Tinf = 500, 
          TpIn_0 = 580.41, 
          TpOut_0 = 560, 
          TsIn_0 = 500, 
          TsOut_0 = 507.84) annotation(
          Placement(transformation(origin = {-28.2, 32.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
      
        SMD_MSR_Modelica.HeatTransport.UHX uhx(
          Tp_0 = 500, 
          vDot = volDotCoolant, 
          cP = scpCoolant, 
          rho = rhoCoolant, 
          vol = 0.204183209633634)  annotation(
          Placement(transformation(origin = {77.462, -70.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
          vol = 0.2526, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 19.932, 
          EnableRad = false, 
          Ar = 7.9559, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 507.84)  annotation(
          Placement(transformation(origin = {-28.8, -81.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
          vol = 0.4419, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 34.870, 
          EnableRad = false, 
          Ar = 13.918, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 500) annotation(
          Placement(transformation(origin = {93.2, 7.73333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
          vol = 0.04, 
          rho = rhoFuel, 
          cP = scpFuel, 
          vDotNom = volDotFuel, 
          DHRS_tK = 10, 
          DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
          DHRS_P_Bleed = 0, 
          DHRS_time = 1000000, 
          K = kFuel, 
          Ac = 0.282743339, 
          L = 0.141471061, 
          Ar = 0.266666667, 
          EnableRad = false, 
          e = 0.08, 
          Tinf = 644.5, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-267, 32.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
          vol = 4.8738e-03, 
          volFracNode = 0.01, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 1.7839, 
          EnableRad = false, 
          Ar = 0.3305, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-135.6, 23.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
          vol = 7.3107e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 2.6758, 
          EnableRad = false, 
          Ar = 0.4958, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 560) annotation(
          Placement(transformation(origin = {-176.4, -61.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
          vol = 2.4369e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 0.8919, 
          EnableRad = false, 
          Ar = 0.1653, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-382.4, -32.8}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        
        MSRR.Components.MSRR1R core1R(
          numGroups = 6, 
          addRho0 = true, 
          EnableRad = false, 
          lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
          beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
          LAMBDA = 2.400E-04, 
          a_F = -6.26E-5,
          a_G = -5.16E-5,
          n_0 = 1, 
          rho_fuel = rhoFuel, 
          rho_grap = rhoGrap, 
          cP_fuel = scpFuel, 
          cP_grap = scpGrap,
          volDotFuel = volDotFuel, T_inf = 570, kFuel = kFuel, TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 570.0000028)  annotation(
          Placement(transformation(origin = {-519.334, -66.0002}, extent = {{-110, -89.9998}, {-49.9999, -29.9999}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
          Placement(transformation(origin = {-27, -161}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
          numRampUp = 1, 
          rampUpK = {1}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {-554, 158}, extent = {{-24, -32}, {24, 16}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
          numRampUp = 1, 
          rampUpK = {100}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {175, 126.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
          Placement(transformation(origin = {43, -209}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
          numSteps = 2, 
          stepTime = {0, 4000}, 
          amplitude = {0, 10})  annotation(
          Placement(transformation(origin = {-708.198, 63.8017}, extent = {{-29.7983, -29.7983}, {29.7983, 29.7983}})));
      
      equation
      
      connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-58, 20}, {-94, 20}, {-94, -58}, {-48, -58}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-18, -58}, {23, -58}, {23, -70}}, color = {204, 0, 0}, thickness = 1));
      connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{77, -70}, {140, -70}, {140, 30}, {113, 30}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{83, 30}, {60.5, 30}, {60.5, 20}, {52, 20}}, color = {204, 0, 0}, thickness = 1));
      connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-247, 45}, {-205.5, 45}, {-205.5, 55}, {-163, 55}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-121, 55}, {-95, 55}, {-95, 45}, {-58, 45}}, color = {204, 0, 0}, thickness = 1));
      connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{52, 45}, {156, 45}, {156, -93}, {-150, -93}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-367, 0}, {-328.5, 0}, {-328.5, 45}, {-300, 45}}, color = {204, 0, 0}, thickness = 1));
      connect(core1R.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-659, -166}, {-409.5, -166}, {-409.5, 0}, {-411, 0}}, color = {204, 0, 0}, thickness = 1));
      connect(pipeHXtoCore.PiTempOut, core1R.tempIn) annotation(
          Line(points = {{-191, -93}, {-390, -93}, {-390, -205}, {-699, -205}}, color = {204, 0, 0}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-448, -205}, {-448, 10}, {-411, 10}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-659, -205}, {-448, -205}, {-448, 70}, {-300, 70}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-659, 98}, {-163, 98}, {-163, 65}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-659, -205}, {-659, 114}, {-3, 114}, {-3, 54}}, color = {220, 138, 221}, thickness = 1));
      connect(core1R.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-659, -205}, {-353.5, -205}, {-353.5, -102}, {-150, -102}}, color = {220, 138, 221}, thickness = 1));
      connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-29, -145}, {-74, -145}, {-74, -52}, {-48, -52}}, color = {220, 138, 221}, thickness = 1));
      connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-29, -145}, {148, -145}, {148, 37}, {113, 37}}, color = {220, 138, 221}, thickness = 1));
      connect(primaryPump.flowFrac, core1R.flowFracIn) annotation(
          Line(points = {{-554, 126}, {-554, 20}, {-699, 20}, {-699, -166}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-554, 126}, {-468, 126}, {-468, -10}, {-411, -10}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-554, 126}, {-328, 126}, {-328, 20}, {-300, 20}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-554, 126}, {-206, 126}, {-206, 46}, {-163, 46}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-554, 126}, {-554, 92}, {-45, 92}, {-45, 54}}, color = {245, 121, 0}, thickness = 1));
      connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-554, 126}, {-108, 126}, {-108, -83}, {-150, -83}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 96}, {174, 96}, {174, 24}, {113, 24}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 96}, {175, 98}, {174, 98}, {174, -6}, {40, -6}, {40, 11}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 96}, {174, 96}, {174, -43}, {23, -43}}, color = {245, 121, 0}, thickness = 1));
      connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 96}, {175, 98}, {176, 98}, {176, -82}, {-48, -82}, {-48, -65}}, color = {245, 121, 0}, thickness = 1));
      connect(externalReact.step, core1R.realIn) annotation(
          Line(points = {{-706, 74}, {-600, 74}, {-600, -106}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{44, -196}, {36, -196}, {36, -98}}, thickness = 1));
      annotation(
          Diagram(coordinateSystem(extent = {{-760, 180}, {220, -280}}), graphics = {Polygon(points = {{108, 32}, {108, 32}})}));
      end R1MSRR10pcm;
    end R1fullSteps;

    package R9fullSteps
      model R9MSRR2dol
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
      parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
      parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
      parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
      
      parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
      parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
      
        MSRR.Components.HeatExchanger heatExchanger(
          vol_P = 0.045379, 
          vol_T = 0.023343, 
          vol_S = 0.055947, 
          rhoP = rhoFuel, 
          rhoT = rhoHXtube, 
          rhoS = rhoCoolant, 
          cP_P = scpFuel, 
          cP_T = scpHXtube, 
          cP_S = scpCoolant, 
          VdotPnom = volDotFuel, 
          VdotSnom = volDotCoolant, 
          hApNom = 7.7966e+04, 
          hAsNom = 1.8150e+04, 
          EnableRad = false, 
          AcShell = 0.1298, 
          AcTube = 0.020276, 
          ArShell = 3.1165, 
          Kp = kFuel, 
          Ks = kCoolant, 
          L_shell = 6.2592, 
          L_tube = 12.5185, 
          e = 0.01, 
          Tinf = 500, 
          TpIn_0 = 580.41, 
          TpOut_0 = 560, 
          TsIn_0 = 500, 
          TsOut_0 = 507.84) annotation(
          Placement(transformation(origin = {-12.2, 90.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
      
        SMD_MSR_Modelica.HeatTransport.UHX uhx(
          Tp_0 = 500, 
          vDot = volDotCoolant, 
          cP = scpCoolant, 
          rho = rhoCoolant, 
          vol = 0.204183209633634)  annotation(
          Placement(transformation(origin = {103.462, -66.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
          vol = 0.2526, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 19.932, 
          EnableRad = false, 
          Ar = 7.9559, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 507.84)  annotation(
          Placement(transformation(origin = {-28.8, -55.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
          vol = 0.4419, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 34.870, 
          EnableRad = false, 
          Ar = 13.918, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 500) annotation(
          Placement(transformation(origin = {95.2, 27.7333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
          vol = 0.04, 
          rho = rhoFuel, 
          cP = scpFuel, 
          vDotNom = volDotFuel, 
          DHRS_tK = 10, 
          DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
          DHRS_P_Bleed = 0, 
          DHRS_time = 1000000, 
          K = kFuel, 
          Ac = 0.282743339, 
          L = 0.141471061, 
          Ar = 0.266666667, 
          EnableRad = false, 
          e = 0.08, 
          Tinf = 644.5, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-249, 30.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
          vol = 4.8738e-03, 
          volFracNode = 0.01, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 1.7839, 
          EnableRad = false, 
          Ar = 0.3305, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-141.6, 11.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
          vol = 7.3107e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 2.6758, 
          EnableRad = false, 
          Ar = 0.4958, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 560) annotation(
          Placement(transformation(origin = {-158.4, -105.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
          vol = 2.4369e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 0.8919, 
          EnableRad = false, 
          Ar = 0.1653, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation( origin = {-374, 6}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
          Placement(transformation(origin = {-37, -175}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
          numRampUp = 1, 
          rampUpK = {1}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 1000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
          numRampUp = 1, 
          rampUpK = {100}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
          Placement(transformation(origin = {45, -251}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
          numSteps = 2, 
          stepTime = {0, 4000}, 
          amplitude = {0, 1178.138})  annotation(
          Placement(transformation(origin = {-596.198, -38.1983}, extent = {{-23.7983, 23.7983}, {23.7983, -23.7983}}, rotation = -0)));
        
        MSRR.Components.MSRR9R msre9r(
          numGroups = 6, 
          lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
          beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
          LAMBDA = 2.400E-04, 
          n_0 = 1, 
          aF = -6.26E-5, 
          aG = -5.16E-5, 
          addRho0 = true, 
          volF1 = {0.003795391373, 0.012869720861, 0.007038081469, 0.008797601837, 0.021767608195, 0.011889109660, 0.014880331986, 0.059823315476, 0.040594513827}, 
          volF2 = {0.003971456514, 0.008772341956, 0.007038081469, 0.017142787894, 0.014880331986, 0.011889109660, 0.028956494924, 0.034788134316, 0.068118358784}, 
          volG = {0.03488047800, 0.10533760200, 0.08002416000, 0.10244745000, 0.17818736400, 0.13543456200, 0.17330364000, 0.47895127800, 0.46943522400}, 
          volUP = 2*volDotFuel, 
          hA = {513.29, 1430.28, 930.26, 1714.35, 2421.98, 1571.45, 2897.08, 6252.67, 7184.60}, 
          kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, 
          kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, 
          kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, 
          kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, 
          TF1_0 = {566.94, 562.89, 584.55, 609.63, 561.59, 574.98, 590.48, 560.60, 565.19}, 
          TF2_0 = {576.02, 572.69, 597.23, 615.83, 567.64, 582.82, 594.30, 562.63, 566.46}, 
          TG_0 = {570.89, 566.21, 591.18, 613.04, 564.19, 580.21, 593.17, 562.06, 566.66}, 
          Tmix_0 = 580.41, 
          IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, 
          IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, 
          IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, 
          flowFracRegions = {0.061410, 0.138550, 0.234231, 0.565809}, 
          rho_fuel = rhoFuel, 
          cP_fuel = scpFuel, 
          kFuel = kFuel, 
          volDotFuel = volDotFuel,
          rho_grap = rhoGrap, 
          cP_grap = scpGrap, 
          regionTripTime = {200000, 200000, 200000, 200000}, 
          regionCoastDownK = 0.02, 
          freeConvFF = 0.01, 
          EnableRad = false, 
          T_inf = 500, LF1 = {0.7412, 0.3166, 0.1731, 0.2164, 0.3167, 0.1730, 0.2165, 0.4463, 0.3028}, LF2 = {0.7756, 0.2158, 0.1731, 0.4217, 0.2165, 0.1730, 0.4213, 0.2595, 0.5082}, Ac = {0.016189, 0.036524, 0.061748, 0.149159}, ArF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, ArF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, e = 0.08)  annotation(
          Placement(transformation(origin = {-258.333, -15}, extent = {{-179.667, -147}, {-81.6667, -49}})));
      
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-42, 78}, {-94, 78}, {-94, -32}, {-51, -32}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-15, -32}, {28, -32}, {28, -66}, {49, -66}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{103, -66}, {140, -66}, {140, 50}, {118, 50}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{82, 50}, {57.5, 50}, {57.5, 78}, {68, 78}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-229, 43}, {-174, 43}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-122, 43}, {-76, 43}, {-76, 103}, {-42, 103}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{68, 103}, {156, 103}, {156, -137}, {-132, -137}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-359, 39}, {-322.5, 39}, {-322.5, 43}, {-282, 43}}, color = {204, 0, 0}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-39, -159}, {-74, -159}, {-74, -23}, {-51, -23}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-39, -159}, {148, -159}, {148, 59}, {118, 59}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-550, 92}, {-468, 92}, {-468, 29}, {-403, 29}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-550, 92}, {-328, 92}, {-328, 18}, {-282, 18}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-550, 92}, {-206, 92}, {-206, 30}, {-174, 30}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-550, 92}, {-299.5, 92}, {-299.5, 112}, {-29, 112}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-550, 92}, {-108, 92}, {-108, -127}, {-132, -127}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, 41}, {118, 41}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {56, -6}, {56, 69}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, -39}, {49, -39}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-51, -82}, {-51, -41}}, color = {245, 121, 0}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
          Line(points = {{-173, -137}, {-173, -242}, {-552, -242}}, color = {204, 0, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-486, -243}, {-319, -243}, {-319, -146}, {-132, -146}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-486, -243}, {-434, -243}, {-434, 49}, {-403, 49}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-486, -243}, {-318, -243}, {-318, 68}, {-282, 68}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-486, -243}, {-200, -243}, {-200, 56}, {-174, 56}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-486, -243}, {-434, -243}, {-434, 112}, {13, 112}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
          Line(points = {{-550, 92}, {-550, -53}, {-552, -53}, {-552, -178}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-487, -178}, {-486.5, -178}, {-486.5, 39}, {-403, 39}}, color = {204, 0, 0}, thickness = 1));
  connect(externalReact.step, msre9r.realIn) annotation(
          Line(points = {{-594, -46}, {-520, -46}, {-520, -178}}, thickness = 1));
  connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{46, -238}, {62, -238}, {62, -94}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-680, 160}, {220, -320}})));
      end R9MSRR2dol;
      
      model R9MSRR1dol
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
      parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
      parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
      parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
      
      parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
      parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
      
        MSRR.Components.HeatExchanger heatExchanger(
          vol_P = 0.045379, 
          vol_T = 0.023343, 
          vol_S = 0.055947, 
          rhoP = rhoFuel, 
          rhoT = rhoHXtube, 
          rhoS = rhoCoolant, 
          cP_P = scpFuel, 
          cP_T = scpHXtube, 
          cP_S = scpCoolant, 
          VdotPnom = volDotFuel, 
          VdotSnom = volDotCoolant, 
          hApNom = 7.7966e+04, 
          hAsNom = 1.8150e+04, 
          EnableRad = false, 
          AcShell = 0.1298, 
          AcTube = 0.020276, 
          ArShell = 3.1165, 
          Kp = kFuel, 
          Ks = kCoolant, 
          L_shell = 6.2592, 
          L_tube = 12.5185, 
          e = 0.01, 
          Tinf = 500, 
          TpIn_0 = 580.41, 
          TpOut_0 = 560, 
          TsIn_0 = 500, 
          TsOut_0 = 507.84) annotation(
          Placement(transformation(origin = {-10.2, 90.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
      
        SMD_MSR_Modelica.HeatTransport.UHX uhx(
          Tp_0 = 500, 
          vDot = volDotCoolant, 
          cP = scpCoolant, 
          rho = rhoCoolant, 
          vol = 0.204183209633634)  annotation(
          Placement(transformation(origin = {101.462, -66.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
          vol = 0.2526, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 19.932, 
          EnableRad = false, 
          Ar = 7.9559, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 507.84)  annotation(
          Placement(transformation(origin = {-28.8, -55.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
          vol = 0.4419, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 34.870, 
          EnableRad = false, 
          Ar = 13.918, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 500) annotation(
          Placement(transformation(origin = {95.2, 27.7333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
          vol = 0.04, 
          rho = rhoFuel, 
          cP = scpFuel, 
          vDotNom = volDotFuel, 
          DHRS_tK = 10, 
          DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
          DHRS_P_Bleed = 0, 
          DHRS_time = 1000000, 
          K = kFuel, 
          Ac = 0.282743339, 
          L = 0.141471061, 
          Ar = 0.266666667, 
          EnableRad = false, 
          e = 0.08, 
          Tinf = 644.5, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-249, 30.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
          vol = 4.8738e-03, 
          volFracNode = 0.01, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 1.7839, 
          EnableRad = false, 
          Ar = 0.3305, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-141.6, 11.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
          vol = 7.3107e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 2.6758, 
          EnableRad = false, 
          Ar = 0.4958, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 560) annotation(
          Placement(transformation(origin = {-158.4, -105.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
          vol = 2.4369e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 0.8919, 
          EnableRad = false, 
          Ar = 0.1653, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation( origin = {-374, 6}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
          Placement(transformation(origin = {-37, -175}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
          numRampUp = 1, 
          rampUpK = {1}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 1000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
          numRampUp = 1, 
          rampUpK = {100}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
          Placement(transformation(origin = {43, -251}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
          numSteps = 2, 
          stepTime = {0, 4000}, 
          amplitude = {0, 589.069})  annotation(
          Placement(transformation(origin = {-596.198, -38.1983}, extent = {{-23.7983, 23.7983}, {23.7983, -23.7983}}, rotation = -0)));
        
        MSRR.Components.MSRR9R msre9r(
          numGroups = 6, 
          lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
          beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
          LAMBDA = 2.400E-04, 
          n_0 = 1, 
          aF = -6.26E-5, 
          aG = -5.16E-5, 
          addRho0 = true, 
          volF1 = {0.003795391373, 0.012869720861, 0.007038081469, 0.008797601837, 0.021767608195, 0.011889109660, 0.014880331986, 0.059823315476, 0.040594513827}, 
          volF2 = {0.003971456514, 0.008772341956, 0.007038081469, 0.017142787894, 0.014880331986, 0.011889109660, 0.028956494924, 0.034788134316, 0.068118358784}, 
          volG = {0.03488047800, 0.10533760200, 0.08002416000, 0.10244745000, 0.17818736400, 0.13543456200, 0.17330364000, 0.47895127800, 0.46943522400}, 
          volUP = 2*volDotFuel, 
          hA = {513.29, 1430.28, 930.26, 1714.35, 2421.98, 1571.45, 2897.08, 6252.67, 7184.60}, 
          kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, 
          kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, 
          kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, 
          kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, 
          TF1_0 = {566.94, 562.89, 584.55, 609.63, 561.59, 574.98, 590.48, 560.60, 565.19}, 
          TF2_0 = {576.02, 572.69, 597.23, 615.83, 567.64, 582.82, 594.30, 562.63, 566.46}, 
          TG_0 = {570.89, 566.21, 591.18, 613.04, 564.19, 580.21, 593.17, 562.06, 566.66}, 
          Tmix_0 = 580.41, 
          IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, 
          IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, 
          IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, 
          flowFracRegions = {0.061410, 0.138550, 0.234231, 0.565809}, 
          rho_fuel = rhoFuel, 
          cP_fuel = scpFuel, 
          kFuel = kFuel, 
          volDotFuel = volDotFuel,
          rho_grap = rhoGrap, 
          cP_grap = scpGrap, 
          regionTripTime = {200000, 200000, 200000, 200000}, 
          regionCoastDownK = 0.02, 
          freeConvFF = 0.01, 
          EnableRad = false, 
          T_inf = 500, LF1 = {0.7412, 0.3166, 0.1731, 0.2164, 0.3167, 0.1730, 0.2165, 0.4463, 0.3028}, LF2 = {0.7756, 0.2158, 0.1731, 0.4217, 0.2165, 0.1730, 0.4213, 0.2595, 0.5082}, Ac = {0.016189, 0.036524, 0.061748, 0.149159}, ArF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, ArF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, e = 0.08)  annotation(
          Placement(transformation(origin = {-260.333, -15}, extent = {{-179.667, -147}, {-81.6667, -49}})));
      
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-40, 78}, {-94, 78}, {-94, -32}, {-51, -32}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-15, -32}, {28, -32}, {28, -66}, {47, -66}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{101, -66}, {140, -66}, {140, 50}, {118, 50}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{82, 50}, {57.5, 50}, {57.5, 78}, {70, 78}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-229, 43}, {-174, 43}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-122, 43}, {-76, 43}, {-76, 103}, {-40, 103}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{70, 103}, {156, 103}, {156, -137}, {-132, -137}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-359, 39}, {-322.5, 39}, {-322.5, 43}, {-282, 43}}, color = {204, 0, 0}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-39, -159}, {-74, -159}, {-74, -23}, {-51, -23}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-39, -159}, {148, -159}, {148, 59}, {118, 59}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-550, 92}, {-468, 92}, {-468, 29}, {-403, 29}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-550, 92}, {-328, 92}, {-328, 18}, {-282, 18}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-550, 92}, {-206, 92}, {-206, 30}, {-174, 30}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-550, 92}, {-299.5, 92}, {-299.5, 112}, {-27, 112}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-550, 92}, {-108, 92}, {-108, -127}, {-132, -127}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, 41}, {118, 41}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {58, -6}, {58, 69}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, -39}, {47, -39}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-51, -82}, {-51, -41}}, color = {245, 121, 0}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
          Line(points = {{-173, -137}, {-175, -242}, {-554, -242}}, color = {204, 0, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-488, -243}, {-319, -243}, {-319, -146}, {-132, -146}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-488, -243}, {-434, -243}, {-434, 49}, {-403, 49}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-488, -243}, {-318, -243}, {-318, 68}, {-282, 68}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-488, -243}, {-200, -243}, {-200, 56}, {-174, 56}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-488, -243}, {-434, -243}, {-434, 112}, {15, 112}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
          Line(points = {{-550, 92}, {-550, -53}, {-554, -53}, {-554, -178}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-489, -178}, {-486.5, -178}, {-486.5, 39}, {-403, 39}}, color = {204, 0, 0}, thickness = 1));
      connect(externalReact.step, msre9r.realIn) annotation(
          Line(points = {{-594, -46}, {-520, -46}, {-522, -178}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{44, -237.5}, {62, -237.5}, {62, -94}, {47, -94}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-680, 160}, {220, -320}})));
      end R9MSRR1dol;
      
      model R9MSRRhalfDol
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
      parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
      parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
      parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
      
      parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
      parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
      
        MSRR.Components.HeatExchanger heatExchanger(
          vol_P = 0.045379, 
          vol_T = 0.023343, 
          vol_S = 0.055947, 
          rhoP = rhoFuel, 
          rhoT = rhoHXtube, 
          rhoS = rhoCoolant, 
          cP_P = scpFuel, 
          cP_T = scpHXtube, 
          cP_S = scpCoolant, 
          VdotPnom = volDotFuel, 
          VdotSnom = volDotCoolant, 
          hApNom = 7.7966e+04, 
          hAsNom = 1.8150e+04, 
          EnableRad = false, 
          AcShell = 0.1298, 
          AcTube = 0.020276, 
          ArShell = 3.1165, 
          Kp = kFuel, 
          Ks = kCoolant, 
          L_shell = 6.2592, 
          L_tube = 12.5185, 
          e = 0.01, 
          Tinf = 500, 
          TpIn_0 = 580.41, 
          TpOut_0 = 560, 
          TsIn_0 = 500, 
          TsOut_0 = 507.84) annotation(
          Placement(transformation(origin = {-10.2, 90.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
      
        SMD_MSR_Modelica.HeatTransport.UHX uhx(
          Tp_0 = 500, 
          vDot = volDotCoolant, 
          cP = scpCoolant, 
          rho = rhoCoolant, 
          vol = 0.204183209633634)  annotation(
          Placement(transformation(origin = {101.462, -66.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
          vol = 0.2526, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 19.932, 
          EnableRad = false, 
          Ar = 7.9559, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 507.84)  annotation(
          Placement(transformation(origin = {-28.8, -55.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
          vol = 0.4419, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 34.870, 
          EnableRad = false, 
          Ar = 13.918, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 500) annotation(
          Placement(transformation(origin = {95.2, 27.7333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
          vol = 0.04, 
          rho = rhoFuel, 
          cP = scpFuel, 
          vDotNom = volDotFuel, 
          DHRS_tK = 10, 
          DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
          DHRS_P_Bleed = 0, 
          DHRS_time = 1000000, 
          K = kFuel, 
          Ac = 0.282743339, 
          L = 0.141471061, 
          Ar = 0.266666667, 
          EnableRad = false, 
          e = 0.08, 
          Tinf = 644.5, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-249, 30.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
          vol = 4.8738e-03, 
          volFracNode = 0.01, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 1.7839, 
          EnableRad = false, 
          Ar = 0.3305, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-141.6, 11.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
          vol = 7.3107e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 2.6758, 
          EnableRad = false, 
          Ar = 0.4958, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 560) annotation(
          Placement(transformation(origin = {-158.4, -105.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
          vol = 2.4369e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 0.8919, 
          EnableRad = false, 
          Ar = 0.1653, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation( origin = {-374, 6}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
          Placement(transformation(origin = {-37, -175}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
          numRampUp = 1, 
          rampUpK = {1}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 1000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
          numRampUp = 1, 
          rampUpK = {100}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
          Placement(transformation(origin = {43, -251}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
          numSteps = 2, 
          stepTime = {0, 4000}, 
          amplitude = {0, 294.535})  annotation(
          Placement(transformation(origin = {-596.198, -38.1983}, extent = {{-23.7983, 23.7983}, {23.7983, -23.7983}}, rotation = -0)));
        
        MSRR.Components.MSRR9R msre9r(
          numGroups = 6, 
          lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
          beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
          LAMBDA = 2.400E-04, 
          n_0 = 1, 
          aF = -6.26E-5, 
          aG = -5.16E-5, 
          addRho0 = true, 
          volF1 = {0.003795391373, 0.012869720861, 0.007038081469, 0.008797601837, 0.021767608195, 0.011889109660, 0.014880331986, 0.059823315476, 0.040594513827}, 
          volF2 = {0.003971456514, 0.008772341956, 0.007038081469, 0.017142787894, 0.014880331986, 0.011889109660, 0.028956494924, 0.034788134316, 0.068118358784}, 
          volG = {0.03488047800, 0.10533760200, 0.08002416000, 0.10244745000, 0.17818736400, 0.13543456200, 0.17330364000, 0.47895127800, 0.46943522400}, 
          volUP = 2*volDotFuel, 
          hA = {513.29, 1430.28, 930.26, 1714.35, 2421.98, 1571.45, 2897.08, 6252.67, 7184.60}, 
          kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, 
          kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, 
          kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, 
          kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, 
          TF1_0 = {566.94, 562.89, 584.55, 609.63, 561.59, 574.98, 590.48, 560.60, 565.19}, 
          TF2_0 = {576.02, 572.69, 597.23, 615.83, 567.64, 582.82, 594.30, 562.63, 566.46}, 
          TG_0 = {570.89, 566.21, 591.18, 613.04, 564.19, 580.21, 593.17, 562.06, 566.66}, 
          Tmix_0 = 580.41, 
          IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, 
          IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, 
          IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, 
          flowFracRegions = {0.061410, 0.138550, 0.234231, 0.565809}, 
          rho_fuel = rhoFuel, 
          cP_fuel = scpFuel, 
          kFuel = kFuel, 
          volDotFuel = volDotFuel,
          rho_grap = rhoGrap, 
          cP_grap = scpGrap, 
          regionTripTime = {200000, 200000, 200000, 200000}, 
          regionCoastDownK = 0.02, 
          freeConvFF = 0.01, 
          EnableRad = false, 
          T_inf = 500, LF1 = {0.7412, 0.3166, 0.1731, 0.2164, 0.3167, 0.1730, 0.2165, 0.4463, 0.3028}, LF2 = {0.7756, 0.2158, 0.1731, 0.4217, 0.2165, 0.1730, 0.4213, 0.2595, 0.5082}, Ac = {0.016189, 0.036524, 0.061748, 0.149159}, ArF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, ArF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, e = 0.08)  annotation(
          Placement(transformation(origin = {-260.333, -15}, extent = {{-179.667, -147}, {-81.6667, -49}})));
      
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-40, 78}, {-94, 78}, {-94, -32}, {-51, -32}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-15, -32}, {28, -32}, {28, -66}, {47, -66}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{101, -66}, {140, -66}, {140, 50}, {118, 50}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{82, 50}, {57.5, 50}, {57.5, 78}, {70, 78}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-229, 43}, {-174, 43}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-122, 43}, {-76, 43}, {-76, 103}, {-40, 103}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{70, 103}, {156, 103}, {156, -137}, {-132, -137}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-359, 39}, {-322.5, 39}, {-322.5, 43}, {-282, 43}}, color = {204, 0, 0}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-39, -159}, {-74, -159}, {-74, -23}, {-51, -23}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-39, -159}, {148, -159}, {148, 59}, {118, 59}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-550, 92}, {-468, 92}, {-468, 29}, {-403, 29}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-550, 92}, {-328, 92}, {-328, 18}, {-282, 18}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-550, 92}, {-206, 92}, {-206, 30}, {-174, 30}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-550, 92}, {-299.5, 92}, {-299.5, 112}, {-27, 112}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-550, 92}, {-108, 92}, {-108, -127}, {-132, -127}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, 41}, {118, 41}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {58, -6}, {58, 69}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, -39}, {47, -39}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-51, -82}, {-51, -41}}, color = {245, 121, 0}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
          Line(points = {{-173, -137}, {-175, -242}, {-554, -242}}, color = {204, 0, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-488, -243}, {-319, -243}, {-319, -146}, {-132, -146}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-488, -243}, {-434, -243}, {-434, 49}, {-403, 49}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-488, -243}, {-318, -243}, {-318, 68}, {-282, 68}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-488, -243}, {-200, -243}, {-200, 56}, {-174, 56}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-488, -243}, {-434, -243}, {-434, 112}, {15, 112}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
          Line(points = {{-550, 92}, {-550, -53}, {-554, -53}, {-554, -178}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-489, -178}, {-486.5, -178}, {-486.5, 39}, {-403, 39}}, color = {204, 0, 0}, thickness = 1));
      connect(externalReact.step, msre9r.realIn) annotation(
          Line(points = {{-594, -46}, {-520, -46}, {-522, -178}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{44, -237.5}, {62, -237.5}, {62, -94}, {47, -94}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-680, 160}, {220, -320}})));
      end R9MSRRhalfDol;
      
      model R9MSRRpOneDol
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
      parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
      parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
      parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
      
      parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
      parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
      
        MSRR.Components.HeatExchanger heatExchanger(
          vol_P = 0.045379, 
          vol_T = 0.023343, 
          vol_S = 0.055947, 
          rhoP = rhoFuel, 
          rhoT = rhoHXtube, 
          rhoS = rhoCoolant, 
          cP_P = scpFuel, 
          cP_T = scpHXtube, 
          cP_S = scpCoolant, 
          VdotPnom = volDotFuel, 
          VdotSnom = volDotCoolant, 
          hApNom = 7.7966e+04, 
          hAsNom = 1.8150e+04, 
          EnableRad = false, 
          AcShell = 0.1298, 
          AcTube = 0.020276, 
          ArShell = 3.1165, 
          Kp = kFuel, 
          Ks = kCoolant, 
          L_shell = 6.2592, 
          L_tube = 12.5185, 
          e = 0.01, 
          Tinf = 500, 
          TpIn_0 = 580.41, 
          TpOut_0 = 560, 
          TsIn_0 = 500, 
          TsOut_0 = 507.84) annotation(
          Placement(transformation(origin = {-10.2, 90.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
      
        SMD_MSR_Modelica.HeatTransport.UHX uhx(
          Tp_0 = 500, 
          vDot = volDotCoolant, 
          cP = scpCoolant, 
          rho = rhoCoolant, 
          vol = 0.204183209633634)  annotation(
          Placement(transformation(origin = {101.462, -66.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
          vol = 0.2526, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 19.932, 
          EnableRad = false, 
          Ar = 7.9559, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 507.84)  annotation(
          Placement(transformation(origin = {-28.8, -55.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
          vol = 0.4419, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 34.870, 
          EnableRad = false, 
          Ar = 13.918, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 500) annotation(
          Placement(transformation(origin = {95.2, 27.7333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
          vol = 0.04, 
          rho = rhoFuel, 
          cP = scpFuel, 
          vDotNom = volDotFuel, 
          DHRS_tK = 10, 
          DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
          DHRS_P_Bleed = 0, 
          DHRS_time = 1000000, 
          K = kFuel, 
          Ac = 0.282743339, 
          L = 0.141471061, 
          Ar = 0.266666667, 
          EnableRad = false, 
          e = 0.08, 
          Tinf = 644.5, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-249, 30.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
          vol = 4.8738e-03, 
          volFracNode = 0.01, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 1.7839, 
          EnableRad = false, 
          Ar = 0.3305, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-141.6, 11.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
          vol = 7.3107e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 2.6758, 
          EnableRad = false, 
          Ar = 0.4958, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 560) annotation(
          Placement(transformation(origin = {-158.4, -105.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
          vol = 2.4369e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 0.8919, 
          EnableRad = false, 
          Ar = 0.1653, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation( origin = {-374, 6}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
          Placement(transformation(origin = {-37, -175}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
          numRampUp = 1, 
          rampUpK = {1}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 1000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
          numRampUp = 1, 
          rampUpK = {100}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
          Placement(transformation(origin = {43, -251}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
          numSteps = 2, 
          stepTime = {0, 4000}, 
          amplitude = {0, 58.907})  annotation(
          Placement(transformation(origin = {-596.198, -38.1983}, extent = {{-23.7983, 23.7983}, {23.7983, -23.7983}}, rotation = -0)));
        
        MSRR.Components.MSRR9R msre9r(
          numGroups = 6, 
          lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
          beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
          LAMBDA = 2.400E-04, 
          n_0 = 1, 
          aF = -6.26E-5, 
          aG = -5.16E-5, 
          addRho0 = true, 
          volF1 = {0.003795391373, 0.012869720861, 0.007038081469, 0.008797601837, 0.021767608195, 0.011889109660, 0.014880331986, 0.059823315476, 0.040594513827}, 
          volF2 = {0.003971456514, 0.008772341956, 0.007038081469, 0.017142787894, 0.014880331986, 0.011889109660, 0.028956494924, 0.034788134316, 0.068118358784}, 
          volG = {0.03488047800, 0.10533760200, 0.08002416000, 0.10244745000, 0.17818736400, 0.13543456200, 0.17330364000, 0.47895127800, 0.46943522400}, 
          volUP = 2*volDotFuel, 
          hA = {513.29, 1430.28, 930.26, 1714.35, 2421.98, 1571.45, 2897.08, 6252.67, 7184.60}, 
          kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, 
          kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, 
          kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, 
          kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, 
          TF1_0 = {566.94, 562.89, 584.55, 609.63, 561.59, 574.98, 590.48, 560.60, 565.19}, 
          TF2_0 = {576.02, 572.69, 597.23, 615.83, 567.64, 582.82, 594.30, 562.63, 566.46}, 
          TG_0 = {570.89, 566.21, 591.18, 613.04, 564.19, 580.21, 593.17, 562.06, 566.66}, 
          Tmix_0 = 580.41, 
          IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, 
          IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, 
          IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, 
          flowFracRegions = {0.061410, 0.138550, 0.234231, 0.565809}, 
          rho_fuel = rhoFuel, 
          cP_fuel = scpFuel, 
          kFuel = kFuel, 
          volDotFuel = volDotFuel,
          rho_grap = rhoGrap, 
          cP_grap = scpGrap, 
          regionTripTime = {200000, 200000, 200000, 200000}, 
          regionCoastDownK = 0.02, 
          freeConvFF = 0.01, 
          EnableRad = false, 
          T_inf = 500, LF1 = {0.7412, 0.3166, 0.1731, 0.2164, 0.3167, 0.1730, 0.2165, 0.4463, 0.3028}, LF2 = {0.7756, 0.2158, 0.1731, 0.4217, 0.2165, 0.1730, 0.4213, 0.2595, 0.5082}, Ac = {0.016189, 0.036524, 0.061748, 0.149159}, ArF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, ArF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, e = 0.08)  annotation(
          Placement(transformation(origin = {-260.333, -15}, extent = {{-179.667, -147}, {-81.6667, -49}})));
      
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-40, 78}, {-94, 78}, {-94, -32}, {-51, -32}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-15, -32}, {28, -32}, {28, -66}, {47, -66}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{101, -66}, {140, -66}, {140, 50}, {118, 50}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{82, 50}, {57.5, 50}, {57.5, 78}, {70, 78}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-229, 43}, {-174, 43}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-122, 43}, {-76, 43}, {-76, 103}, {-40, 103}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{70, 103}, {156, 103}, {156, -137}, {-132, -137}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-359, 39}, {-322.5, 39}, {-322.5, 43}, {-282, 43}}, color = {204, 0, 0}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-39, -159}, {-74, -159}, {-74, -23}, {-51, -23}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-39, -159}, {148, -159}, {148, 59}, {118, 59}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-550, 92}, {-468, 92}, {-468, 29}, {-403, 29}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-550, 92}, {-328, 92}, {-328, 18}, {-282, 18}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-550, 92}, {-206, 92}, {-206, 30}, {-174, 30}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-550, 92}, {-299.5, 92}, {-299.5, 112}, {-27, 112}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-550, 92}, {-108, 92}, {-108, -127}, {-132, -127}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, 41}, {118, 41}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {58, -6}, {58, 69}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, -39}, {47, -39}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-51, -82}, {-51, -41}}, color = {245, 121, 0}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
          Line(points = {{-173, -137}, {-175, -242}, {-554, -242}}, color = {204, 0, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-488, -243}, {-319, -243}, {-319, -146}, {-132, -146}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-488, -243}, {-434, -243}, {-434, 49}, {-403, 49}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-488, -243}, {-318, -243}, {-318, 68}, {-282, 68}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-488, -243}, {-200, -243}, {-200, 56}, {-174, 56}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-488, -243}, {-434, -243}, {-434, 112}, {15, 112}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
          Line(points = {{-550, 92}, {-550, -53}, {-554, -53}, {-554, -178}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-489, -178}, {-486.5, -178}, {-486.5, 39}, {-403, 39}}, color = {204, 0, 0}, thickness = 1));
      connect(externalReact.step, msre9r.realIn) annotation(
          Line(points = {{-594, -46}, {-520, -46}, {-522, -178}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{44, -237.5}, {62, -237.5}, {62, -94}, {47, -94}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-680, 160}, {220, -320}})));
      end R9MSRRpOneDol;
      
      model R9MSRR100pcm
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
      parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
      parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
      parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
      
      parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
      parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
      
        MSRR.Components.HeatExchanger heatExchanger(
          vol_P = 0.045379, 
          vol_T = 0.023343, 
          vol_S = 0.055947, 
          rhoP = rhoFuel, 
          rhoT = rhoHXtube, 
          rhoS = rhoCoolant, 
          cP_P = scpFuel, 
          cP_T = scpHXtube, 
          cP_S = scpCoolant, 
          VdotPnom = volDotFuel, 
          VdotSnom = volDotCoolant, 
          hApNom = 7.7966e+04, 
          hAsNom = 1.8150e+04, 
          EnableRad = false, 
          AcShell = 0.1298, 
          AcTube = 0.020276, 
          ArShell = 3.1165, 
          Kp = kFuel, 
          Ks = kCoolant, 
          L_shell = 6.2592, 
          L_tube = 12.5185, 
          e = 0.01, 
          Tinf = 500, 
          TpIn_0 = 580.41, 
          TpOut_0 = 560, 
          TsIn_0 = 500, 
          TsOut_0 = 507.84) annotation(
          Placement(transformation(origin = {-10.2, 90.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
      
        SMD_MSR_Modelica.HeatTransport.UHX uhx(
          Tp_0 = 500, 
          vDot = volDotCoolant, 
          cP = scpCoolant, 
          rho = rhoCoolant, 
          vol = 0.204183209633634)  annotation(
          Placement(transformation(origin = {101.462, -66.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
          vol = 0.2526, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 19.932, 
          EnableRad = false, 
          Ar = 7.9559, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 507.84)  annotation(
          Placement(transformation(origin = {-28.8, -55.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
          vol = 0.4419, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 34.870, 
          EnableRad = false, 
          Ar = 13.918, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 500) annotation(
          Placement(transformation(origin = {95.2, 27.7333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
          vol = 0.04, 
          rho = rhoFuel, 
          cP = scpFuel, 
          vDotNom = volDotFuel, 
          DHRS_tK = 10, 
          DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
          DHRS_P_Bleed = 0, 
          DHRS_time = 1000000, 
          K = kFuel, 
          Ac = 0.282743339, 
          L = 0.141471061, 
          Ar = 0.266666667, 
          EnableRad = false, 
          e = 0.08, 
          Tinf = 644.5, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-249, 30.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
          vol = 4.8738e-03, 
          volFracNode = 0.01, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 1.7839, 
          EnableRad = false, 
          Ar = 0.3305, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-141.6, 11.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
          vol = 7.3107e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 2.6758, 
          EnableRad = false, 
          Ar = 0.4958, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 560) annotation(
          Placement(transformation(origin = {-158.4, -105.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
          vol = 2.4369e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 0.8919, 
          EnableRad = false, 
          Ar = 0.1653, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation( origin = {-374, 6}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
          Placement(transformation(origin = {-37, -175}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
          numRampUp = 1, 
          rampUpK = {1}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 1000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
          numRampUp = 1, 
          rampUpK = {100}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
          Placement(transformation(origin = {43, -251}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
          numSteps = 2, 
          stepTime = {0, 4000}, 
          amplitude = {0, 100})  annotation(
          Placement(transformation(origin = {-596.198, -38.1983}, extent = {{-23.7983, 23.7983}, {23.7983, -23.7983}}, rotation = -0)));
        
        MSRR.Components.MSRR9R msre9r(
          numGroups = 6, 
          lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
          beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
          LAMBDA = 2.400E-04, 
          n_0 = 1, 
          aF = -6.26E-5, 
          aG = -5.16E-5, 
          addRho0 = true, 
          volF1 = {0.003795391373, 0.012869720861, 0.007038081469, 0.008797601837, 0.021767608195, 0.011889109660, 0.014880331986, 0.059823315476, 0.040594513827}, 
          volF2 = {0.003971456514, 0.008772341956, 0.007038081469, 0.017142787894, 0.014880331986, 0.011889109660, 0.028956494924, 0.034788134316, 0.068118358784}, 
          volG = {0.03488047800, 0.10533760200, 0.08002416000, 0.10244745000, 0.17818736400, 0.13543456200, 0.17330364000, 0.47895127800, 0.46943522400}, 
          volUP = 2*volDotFuel, 
          hA = {513.29, 1430.28, 930.26, 1714.35, 2421.98, 1571.45, 2897.08, 6252.67, 7184.60}, 
          kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, 
          kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, 
          kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, 
          kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, 
          TF1_0 = {566.94, 562.89, 584.55, 609.63, 561.59, 574.98, 590.48, 560.60, 565.19}, 
          TF2_0 = {576.02, 572.69, 597.23, 615.83, 567.64, 582.82, 594.30, 562.63, 566.46}, 
          TG_0 = {570.89, 566.21, 591.18, 613.04, 564.19, 580.21, 593.17, 562.06, 566.66}, 
          Tmix_0 = 580.41, 
          IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, 
          IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, 
          IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, 
          flowFracRegions = {0.061410, 0.138550, 0.234231, 0.565809}, 
          rho_fuel = rhoFuel, 
          cP_fuel = scpFuel, 
          kFuel = kFuel, 
          volDotFuel = volDotFuel,
          rho_grap = rhoGrap, 
          cP_grap = scpGrap, 
          regionTripTime = {200000, 200000, 200000, 200000}, 
          regionCoastDownK = 0.02, 
          freeConvFF = 0.01, 
          EnableRad = false, 
          T_inf = 500, LF1 = {0.7412, 0.3166, 0.1731, 0.2164, 0.3167, 0.1730, 0.2165, 0.4463, 0.3028}, LF2 = {0.7756, 0.2158, 0.1731, 0.4217, 0.2165, 0.1730, 0.4213, 0.2595, 0.5082}, Ac = {0.016189, 0.036524, 0.061748, 0.149159}, ArF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, ArF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, e = 0.08)  annotation(
          Placement(transformation(origin = {-260.333, -15}, extent = {{-179.667, -147}, {-81.6667, -49}})));
      
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-40, 78}, {-94, 78}, {-94, -32}, {-51, -32}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-15, -32}, {28, -32}, {28, -66}, {47, -66}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{101, -66}, {140, -66}, {140, 50}, {118, 50}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{82, 50}, {57.5, 50}, {57.5, 78}, {70, 78}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-229, 43}, {-174, 43}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-122, 43}, {-76, 43}, {-76, 103}, {-40, 103}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{70, 103}, {156, 103}, {156, -137}, {-132, -137}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-359, 39}, {-322.5, 39}, {-322.5, 43}, {-282, 43}}, color = {204, 0, 0}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-39, -159}, {-74, -159}, {-74, -23}, {-51, -23}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-39, -159}, {148, -159}, {148, 59}, {118, 59}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-550, 92}, {-468, 92}, {-468, 29}, {-403, 29}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-550, 92}, {-328, 92}, {-328, 18}, {-282, 18}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-550, 92}, {-206, 92}, {-206, 30}, {-174, 30}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-550, 92}, {-299.5, 92}, {-299.5, 112}, {-27, 112}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-550, 92}, {-108, 92}, {-108, -127}, {-132, -127}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, 41}, {118, 41}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {58, -6}, {58, 69}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, -39}, {47, -39}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-51, -82}, {-51, -41}}, color = {245, 121, 0}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
          Line(points = {{-173, -137}, {-175, -242}, {-554, -242}}, color = {204, 0, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-488, -243}, {-319, -243}, {-319, -146}, {-132, -146}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-488, -243}, {-434, -243}, {-434, 49}, {-403, 49}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-488, -243}, {-318, -243}, {-318, 68}, {-282, 68}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-488, -243}, {-200, -243}, {-200, 56}, {-174, 56}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-488, -243}, {-434, -243}, {-434, 112}, {15, 112}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
          Line(points = {{-550, 92}, {-550, -53}, {-554, -53}, {-554, -178}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-489, -178}, {-486.5, -178}, {-486.5, 39}, {-403, 39}}, color = {204, 0, 0}, thickness = 1));
      connect(externalReact.step, msre9r.realIn) annotation(
          Line(points = {{-594, -46}, {-520, -46}, {-522, -178}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{44, -237.5}, {62, -237.5}, {62, -94}, {47, -94}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-680, 160}, {220, -320}})));
      end R9MSRR100pcm;
      
      model R9MSRR10pcm
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotFuel = 0.011493;
      parameter SMD_MSR_Modelica.Units.Density rhoFuel = 2079.449700;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpFuel = 2009.66;
      parameter SMD_MSR_Modelica.Units.Conductivity kFuel = 0;
      
      parameter SMD_MSR_Modelica.Units.VolumeticFlowRate volDotCoolant = 0.030219;
      parameter SMD_MSR_Modelica.Units.Density rhoCoolant = 1767.1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpCoolant = 2390;
      parameter SMD_MSR_Modelica.Units.Conductivity kCoolant = 0;
      
      parameter SMD_MSR_Modelica.Units.Density rhoGrap = 1800;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpGrap = 1773;
      parameter SMD_MSR_Modelica.Units.Density rhoHXtube = 8774.5;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity scpHXtube = 450;
      
        MSRR.Components.HeatExchanger heatExchanger(
          vol_P = 0.045379, 
          vol_T = 0.023343, 
          vol_S = 0.055947, 
          rhoP = rhoFuel, 
          rhoT = rhoHXtube, 
          rhoS = rhoCoolant, 
          cP_P = scpFuel, 
          cP_T = scpHXtube, 
          cP_S = scpCoolant, 
          VdotPnom = volDotFuel, 
          VdotSnom = volDotCoolant, 
          hApNom = 7.7966e+04, 
          hAsNom = 1.8150e+04, 
          EnableRad = false, 
          AcShell = 0.1298, 
          AcTube = 0.020276, 
          ArShell = 3.1165, 
          Kp = kFuel, 
          Ks = kCoolant, 
          L_shell = 6.2592, 
          L_tube = 12.5185, 
          e = 0.01, 
          Tinf = 500, 
          TpIn_0 = 580.41, 
          TpOut_0 = 560, 
          TsIn_0 = 500, 
          TsOut_0 = 507.84) annotation(
          Placement(transformation(origin = {-10.2, 90.6}, extent = {{-50.8, -25.4}, {76.2, 25.4}})));
      
        SMD_MSR_Modelica.HeatTransport.UHX uhx(
          Tp_0 = 500, 
          vDot = volDotCoolant, 
          cP = scpCoolant, 
          rho = rhoCoolant, 
          vol = 0.204183209633634)  annotation(
          Placement(transformation(origin = {101.462, -66.3463}, extent = {{-54.8618, -41.1463}, {27.4309, 41.1463}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoUHX(
          vol = 0.2526, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 19.932, 
          EnableRad = false, 
          Ar = 7.9559, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 507.84)  annotation(
          Placement(transformation(origin = {-28.8, -55.0667}, extent = {{-27.2, 9.06667}, {18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.Pipe pipeUHXtoHX(
          vol = 0.4419, 
          volFracNode = 0.1, 
          vDotNom = volDotCoolant, 
          rho = rhoCoolant, 
          cP = scpCoolant, 
          K = kCoolant, 
          Ac = 0.012673, 
          L = 34.870, 
          EnableRad = false, 
          Ar = 13.918, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 500) annotation(
          Placement(transformation(origin = {95.2, 27.7333}, extent = {{27.2, 9.06667}, {-18.1333, 36.2667}})));
      
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(
          vol = 0.04, 
          rho = rhoFuel, 
          cP = scpFuel, 
          vDotNom = volDotFuel, 
          DHRS_tK = 10, 
          DHRS_MaxP_Rm(displayUnit = "MW") = 320000, 
          DHRS_P_Bleed = 0, 
          DHRS_time = 1000000, 
          K = kFuel, 
          Ac = 0.282743339, 
          L = 0.141471061, 
          Ar = 0.266666667, 
          EnableRad = false, 
          e = 0.08, 
          Tinf = 644.5, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-249, 30.6667}, extent = {{-49.3333, -24.6667}, {37, 49.3333}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeDHRStoHX(
          vol = 4.8738e-03, 
          volFracNode = 0.01, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 1.7839, 
          EnableRad = false, 
          Ar = 0.3305, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation(origin = {-141.6, 11.2}, extent = {{-38.4, 12.8}, {25.6, 51.2}})));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeHXtoCore(
          vol = 7.3107e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 2.6758, 
          EnableRad = false, 
          Ar = 0.4958, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 560) annotation(
          Placement(transformation(origin = {-158.4, -105.6}, extent = {{-37.2, 12.4}, {24.8, 49.6}}, rotation = 180)));
        
        SMD_MSR_Modelica.HeatTransport.Pipe pipeCoreToDHRS(
          vol = 2.4369e-03, 
          volFracNode = 0.1, 
          vDotNom = volDotFuel, 
          rho = rhoFuel, 
          cP = scpFuel, 
          K = kFuel, 
          Ac = 2.7321e-03, 
          L = 0.8919, 
          EnableRad = false, 
          Ar = 0.1653, 
          e = 0.08, 
          Tinf = 644.4, 
          T_0 = 580.41)  annotation(
          Placement(transformation( origin = {-374, 6}, extent = {{-39.6, 13.2}, {26.4, 52.8}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantVolumetricPower constantVolumetricPower(Q_volumetric = 0)  annotation(
          Placement(transformation(origin = {-37, -175}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump primaryPump(
          numRampUp = 1, 
          rampUpK = {1}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 1000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {-550, 116}, extent = {{-24, -32}, {24, 16}})));
        
        SMD_MSR_Modelica.HeatTransport.Pump secondaryPump(
          numRampUp = 1, 
          rampUpK = {100}, 
          rampUpTo = {1}, 
          rampUpTime = {0}, 
          tripK = 50, 
          tripTime = 10000000, 
          freeConvFF = 0)  annotation(
          Placement(transformation(origin = {175, 114.667}, extent = {{-23, -30.6667}, {23, 15.3333}})));
        
        SMD_MSR_Modelica.Signals.Constants.ConstantReal uhxDemand(magnitude = 1E6)  annotation(
          Placement(transformation(origin = {43, -251}, extent = {{-27, -27}, {27, 27}})));
        
        SMD_MSR_Modelica.Signals.TimeDependent.Stepper externalReact(
          numSteps = 2, 
          stepTime = {0, 4000}, 
          amplitude = {0, 10})  annotation(
          Placement(transformation(origin = {-596.198, -38.1983}, extent = {{-23.7983, 23.7983}, {23.7983, -23.7983}}, rotation = -0)));
        
        MSRR.Components.MSRR9R msre9r(
          numGroups = 6, 
          lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, 
          beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, 
          LAMBDA = 2.400E-04, 
          n_0 = 1, 
          aF = -6.26E-5, 
          aG = -5.16E-5, 
          addRho0 = true, 
          volF1 = {0.003795391373, 0.012869720861, 0.007038081469, 0.008797601837, 0.021767608195, 0.011889109660, 0.014880331986, 0.059823315476, 0.040594513827}, 
          volF2 = {0.003971456514, 0.008772341956, 0.007038081469, 0.017142787894, 0.014880331986, 0.011889109660, 0.028956494924, 0.034788134316, 0.068118358784}, 
          volG = {0.03488047800, 0.10533760200, 0.08002416000, 0.10244745000, 0.17818736400, 0.13543456200, 0.17330364000, 0.47895127800, 0.46943522400}, 
          volUP = 2*volDotFuel, 
          hA = {513.29, 1430.28, 930.26, 1714.35, 2421.98, 1571.45, 2897.08, 6252.67, 7184.60}, 
          kFN1 = {0.014930, 0.027360, 0.045040, 0.051260, 0.036010, 0.060140, 0.068450, 0.061790, 0.093330}, 
          kFN2 = {0.017210, 0.045500, 0.046560, 0.042610, 0.060690, 0.062180, 0.056640, 0.077070, 0.073110}, 
          kHT1 = {9.4600e-04, 1.6850e-03, 3.0290e-03, 3.4470e-03, 2.2160e-03, 4.0440e-03, 4.6030e-03, 3.9200e-03, 6.2770e-03}, 
          kHT2 = {1.0810e-03, 3.0600e-03, 3.1310e-03, 2.3950e-03, 4.0810e-03, 4.1820e-03, 3.1840e-03, 5.1830e-03, 4.3050e-03}, 
          TF1_0 = {566.94, 562.89, 584.55, 609.63, 561.59, 574.98, 590.48, 560.60, 565.19}, 
          TF2_0 = {576.02, 572.69, 597.23, 615.83, 567.64, 582.82, 594.30, 562.63, 566.46}, 
          TG_0 = {570.89, 566.21, 591.18, 613.04, 564.19, 580.21, 593.17, 562.06, 566.66}, 
          Tmix_0 = 580.41, 
          IF1 = {0.021680, 0.021970, 0.078970, 0.082490, 0.022540, 0.082550, 0.086230, 0.027450, 0.069360}, 
          IF2 = {0.026780, 0.065190, 0.084380, 0.041240, 0.068010, 0.088230, 0.042900, 0.055290, 0.034730}, 
          IG = {0.044430, 0.088350, 0.166710, 0.120770, 0.091810, 0.174290, 0.126120, 0.084080, 0.103430}, 
          flowFracRegions = {0.061410, 0.138550, 0.234231, 0.565809}, 
          rho_fuel = rhoFuel, 
          cP_fuel = scpFuel, 
          kFuel = kFuel, 
          volDotFuel = volDotFuel,
          rho_grap = rhoGrap, 
          cP_grap = scpGrap, 
          regionTripTime = {200000, 200000, 200000, 200000}, 
          regionCoastDownK = 0.02, 
          freeConvFF = 0.01, 
          EnableRad = false, 
          T_inf = 500, LF1 = {0.7412, 0.3166, 0.1731, 0.2164, 0.3167, 0.1730, 0.2165, 0.4463, 0.3028}, LF2 = {0.7756, 0.2158, 0.1731, 0.4217, 0.2165, 0.1730, 0.4213, 0.2595, 0.5082}, Ac = {0.016189, 0.036524, 0.061748, 0.149159}, ArF1 = {0.7946, 0.3287, 0.1798, 0.2247, 0.3288, 0.1796, 0.2248, 0.4634, 0.3145}, ArF2 = {0.8314, 0.2241, 0.1798, 0.4379, 0.2248, 0.1796, 0.4374, 0.2695, 0.5277}, e = 0.08)  annotation(
          Placement(transformation(origin = {-368.333, -141}, extent = {{-179.667, -147}, {-81.6667, -49}})));
      
      equation
        connect(heatExchanger.T_out_sFluid, pipeHXtoUHX.PiTemp_IN) annotation(
          Line(points = {{-40, 78}, {-94, 78}, {-94, -32}, {-51, -32}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeHXtoUHX.PiTempOut, uhx.tempIn) annotation(
          Line(points = {{-15, -32}, {28, -32}, {28, -66}, {47, -66}}, color = {204, 0, 0}, thickness = 1));
        connect(uhx.tempOut, pipeUHXtoHX.PiTemp_IN) annotation(
          Line(points = {{101, -66}, {140, -66}, {140, 50}, {118, 50}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeUHXtoHX.PiTempOut, heatExchanger.T_in_sFluid) annotation(
          Line(points = {{82, 50}, {57.5, 50}, {57.5, 78}, {70, 78}}, color = {204, 0, 0}, thickness = 1));
        connect(dhrs.tempOut, pipeDHRStoHX.PiTemp_IN) annotation(
          Line(points = {{-229, 43}, {-174, 43}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeDHRStoHX.PiTempOut, heatExchanger.T_in_pFluid) annotation(
          Line(points = {{-122, 43}, {-76, 43}, {-76, 103}, {-40, 103}}, color = {204, 0, 0}, thickness = 1));
        connect(heatExchanger.T_out_pFluid, pipeHXtoCore.PiTemp_IN) annotation(
          Line(points = {{70, 103}, {156, 103}, {156, -137}, {-132, -137}}, color = {204, 0, 0}, thickness = 1));
        connect(pipeCoreToDHRS.PiTempOut, dhrs.tempIn) annotation(
          Line(points = {{-359, 39}, {-322.5, 39}, {-322.5, 43}, {-282, 43}}, color = {204, 0, 0}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeHXtoUHX.PiDecay_Heat) annotation(
          Line(points = {{-39, -159}, {-74, -159}, {-74, -23}, {-51, -23}}, color = {220, 138, 221}, thickness = 1));
        connect(constantVolumetricPower.volPow, pipeUHXtoHX.PiDecay_Heat) annotation(
          Line(points = {{-39, -159}, {148, -159}, {148, 59}, {118, 59}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, pipeCoreToDHRS.flowFrac) annotation(
          Line(points = {{-550, 92}, {-468, 92}, {-468, 29}, {-403, 29}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, dhrs.flowFrac) annotation(
          Line(points = {{-550, 92}, {-328, 92}, {-328, 18}, {-282, 18}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeDHRStoHX.flowFrac) annotation(
          Line(points = {{-550, 92}, {-206, 92}, {-206, 30}, {-174, 30}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, heatExchanger.primaryFF) annotation(
          Line(points = {{-550, 92}, {-299.5, 92}, {-299.5, 112}, {-27, 112}}, color = {245, 121, 0}, thickness = 1));
        connect(primaryPump.flowFrac, pipeHXtoCore.flowFrac) annotation(
          Line(points = {{-550, 92}, {-108, 92}, {-108, -127}, {-132, -127}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeUHXtoHX.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, 41}, {118, 41}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, heatExchanger.secondaryFF) annotation(
          Line(points = {{175, 92}, {175, 98}, {174, 98}, {174, -6}, {58, -6}, {58, 69}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, uhx.flowFrac) annotation(
          Line(points = {{175, 92}, {174, 92}, {174, -39}, {47, -39}}, color = {245, 121, 0}, thickness = 1));
        connect(secondaryPump.flowFrac, pipeHXtoUHX.flowFrac) annotation(
          Line(points = {{175, 92}, {175, 98}, {176, 98}, {176, -82}, {-51, -82}, {-51, -41}}, color = {245, 121, 0}, thickness = 1));
        connect(pipeHXtoCore.PiTempOut, msre9r.tempIn) annotation(
          Line(points = {{-173, -137}, {-173, -368}, {-662, -368}}, color = {204, 0, 0}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeHXtoCore.PiDecay_Heat) annotation(
          Line(points = {{-596, -369}, {-319, -369}, {-319, -146}, {-132, -146}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeCoreToDHRS.PiDecay_Heat) annotation(
          Line(points = {{-596, -369}, {-596, 49}, {-403, 49}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, dhrs.pDecay) annotation(
          Line(points = {{-596, -369}, {-318, -369}, {-318, 68}, {-282, 68}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, pipeDHRStoHX.PiDecay_Heat) annotation(
          Line(points = {{-596, -369}, {-200, -369}, {-200, 56}, {-174, 56}}, color = {220, 138, 221}, thickness = 1));
        connect(msre9r.volumetircPowerOut, heatExchanger.P_decay) annotation(
          Line(points = {{-596, -369}, {-596, 112}, {15, 112}}, color = {220, 138, 221}, thickness = 1));
        connect(primaryPump.flowFrac, msre9r.flowFracIn) annotation(
          Line(points = {{-550, 92}, {-550, -53}, {-662, -53}, {-662, -304}}, color = {245, 121, 0}, thickness = 1));
        connect(msre9r.tempOut, pipeCoreToDHRS.PiTemp_IN) annotation(
          Line(points = {{-597, -304}, {-486.5, -304}, {-486.5, 39}, {-403, 39}}, color = {204, 0, 0}, thickness = 1));
      connect(externalReact.step, msre9r.realIn) annotation(
          Line(points = {{-594, -46}, {-520, -46}, {-520, -304}, {-630, -304}}, thickness = 1));
      connect(uhxDemand.realOut, uhx.powDemand) annotation(
          Line(points = {{44, -237.5}, {62, -237.5}, {62, -94}, {47, -94}}, thickness = 1));
        annotation(
          Diagram(coordinateSystem(extent = {{-680, 160}, {220, -320}})));
      end R9MSRR10pcm;
    end R9fullSteps;
  end Transients;
end MSRR;