package MSRR
  
  model MSRRCoreRegion
    parameter SMD_MSR_Modelica.Units.Convection hAnom;
    parameter SMD_MSR_Modelica.Units.Mass m_FN1;
    parameter SMD_MSR_Modelica.Units.Mass m_FN2;
    parameter SMD_MSR_Modelica.Units.Mass m_GN;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_fuel;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_graphite;
    parameter SMD_MSR_Modelica.Units.MassFlowRate mdot_fuelNom;
    parameter SMD_MSR_Modelica.Units.VolumeImportance kFN1;
    parameter SMD_MSR_Modelica.Units.VolumeImportance kFN2;
    parameter SMD_MSR_Modelica.Units.VolumeImportance kG;
    parameter SMD_MSR_Modelica.Units.HeatTransferFraction kHT_FN1;
    parameter SMD_MSR_Modelica.Units.HeatTransferFraction kHT_FN2;
    parameter SMD_MSR_Modelica.Units.Temperature TF1_0;
    parameter SMD_MSR_Modelica.Units.Temperature TF2_0;
    parameter SMD_MSR_Modelica.Units.Temperature TG_0;
    parameter SMD_MSR_Modelica.Units.FlowFraction regionFlowFrac;
    parameter SMD_MSR_Modelica.Units.Conductivity KF;
    parameter SMD_MSR_Modelica.Units.Area Ac;
    parameter SMD_MSR_Modelica.Units.Area ArF1;
    parameter SMD_MSR_Modelica.Units.Area ArF2;
    parameter SMD_MSR_Modelica.Units.Length LF1;
    parameter SMD_MSR_Modelica.Units.Length LF2;
    parameter SMD_MSR_Modelica.Units.Emissivity e;
    parameter SMD_MSR_Modelica.Units.Temperature Tinf;
    parameter Boolean OutterRegion;
    constant SMD_MSR_Modelica.Units.StefanBoltzmannConstant SigSBK = 5.670374419E-14;
  
    SMD_MSR_Modelica.Units.MassFlowRate mdot_fuel;
    SMD_MSR_Modelica.Units.Convection hA;
    SMD_MSR_Modelica.Units.RadPower radF1;
    SMD_MSR_Modelica.Units.RadPower radF2;
    input SMD_MSR_Modelica.PortsConnectors.volSpecificHeat_in decayHeat annotation(
      Placement(visible = true, transformation(origin = {-72, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-72, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    output SMD_MSR_Modelica.PortsConnectors.Temp_Out grapNode annotation(
      Placement(visible = true, transformation(origin = {40, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {40, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    output SMD_MSR_Modelica.PortsConnectors.Temp_Out fuelNode1 annotation(
      Placement(visible = true, transformation(origin = {40, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {40, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    input SMD_MSR_Modelica.PortsConnectors.Temp_In temp_In annotation(
      Placement(visible = true, transformation(origin = {-72, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-72, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    input SMD_MSR_Modelica.PortsConnectors.PowerIn fissionPower annotation(
      Placement(visible = true, transformation(origin = {-72, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-72, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    output SMD_MSR_Modelica.PortsConnectors.Temp_Out fuelNode2 annotation(
      Placement(visible = true, transformation(origin = {40, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {40, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In fuelFlowFraction annotation(
      Placement(visible = true, transformation(origin = {-70, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  initial equation
    fuelNode2.T = TF2_0;
    fuelNode1.T = TF1_0;
    grapNode.T = TG_0;
  equation
    hA = 0.02492*fuelFlowFraction.FF^(0.33);
    
    if OutterRegion == true then
      radF1 = e*SigSBK*ArF1*((Tinf + 273.15)^4 - (fuelNode1.T + 273.15)^4);
      radF2 = e*SigSBK*ArF2*((Tinf + 273.15)^4 - (fuelNode1.T + 273.15)^4);
    else
      radF1 = 0;
      radF2 = 0;
    end if;
    
    mdot_fuel = mdot_fuelNom*fuelFlowFraction.FF*regionFlowFrac;
    if fuelFlowFraction.FF > 0 then
      m_FN1*cP_fuel*der(fuelNode1.T) = mdot_fuel*cP_fuel*(temp_In.T - fuelNode1.T) + kFN1*fissionPower.P - hA*(kHT_FN1/(kHT_FN1 + kHT_FN2))*(fuelNode1.T - grapNode.T) + decayHeat.Q*m_FN1 + radF1;
      m_FN2*cP_fuel*der(fuelNode2.T) = mdot_fuel*cP_fuel*(fuelNode1.T - fuelNode2.T) + kFN2*fissionPower.P - hA*(kHT_FN2/(kHT_FN1 + kHT_FN2))*(fuelNode1.T - grapNode.T) + decayHeat.Q*m_FN2 + radF2;
    else
      m_FN1*cP_fuel*der(fuelNode1.T) = ((KF*Ac)/LF1)*(temp_In.T - fuelNode1.T) - hA*(kHT_FN1/(kHT_FN1 + kHT_FN2))*(fuelNode1.T - grapNode.T) + kFN1*fissionPower.P + decayHeat.Q*m_FN1 + radF1;
      m_FN2*cP_fuel*der(fuelNode2.T) = ((KF*Ac)/LF2)*(fuelNode1.T - fuelNode2.T) - hA*(kHT_FN1/(kHT_FN1 + kHT_FN2))*(fuelNode1.T - grapNode.T) + kFN2*fissionPower.P + decayHeat.Q*m_FN2 + radF2;
    end if;
    m_GN*cP_graphite*der(grapNode.T) = hA*(fuelNode1.T - grapNode.T) + kG*fissionPower.P;
    annotation(
      Diagram(coordinateSystem(extent = {{-100, 80}, {60, -100}}), graphics = {Text(origin = {9, 68}, extent = {{-47, 12}, {47, -12}}, textString = "CoreRegion"), Rectangle(origin = {-15.43, -10.01}, lineThickness = 2, extent = {{-75.81, 90.1}, {75.81, -90.1}}), Text(origin = {44, -76}, extent = {{-13, 13}, {13, -13}}, textString = "T_FN1"), Rectangle(origin = {-34, 20}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}}), Text(origin = {44, 24}, extent = {{-13, 13}, {13, -13}}, textString = "T_FN2"), Rectangle(origin = {-15.43, -10.01}, lineThickness = 2, extent = {{-75.81, 90.1}, {75.81, -90.1}}), Rectangle(origin = {-15.43, -10.01}, lineThickness = 2, extent = {{-75.81, 90.1}, {75.81, -90.1}}), Text(origin = {-71, 36}, extent = {{-6, 8}, {6, -8}}, textString = "FH"), Text(origin = {-71, -45}, extent = {{-6, 8}, {6, -8}}, textString = "FF"), Rectangle(origin = {11, -10}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-15, 50}, {15, -50}}), Rectangle(origin = {-34, -40}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}}), Text(origin = {-71, -5}, extent = {{-6, 8}, {6, -8}}, textString = "DH"), Text(origin = {-70, -86}, extent = {{-9, 10}, {9, -10}}, textString = "T_In"), Text(origin = {42, -25}, extent = {{-11, 12}, {11, -12}}, textString = "T_G1")}),
      Icon(coordinateSystem(extent = {{-100, 80}, {60, -100}}), graphics = {Text(origin = {-73, 36}, extent = {{-6, 8}, {6, -8}}, textString = "FH"), Text(origin = {-70, -86}, extent = {{-9, 10}, {9, -10}}, textString = "T_In"), Rectangle(origin = {-34, 20}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}}), Text(origin = {42, -25}, extent = {{-11, 12}, {11, -12}}, textString = "T_G1"), Text(origin = {44, -76}, extent = {{-13, 13}, {13, -13}}, textString = "T_FN1"), Text(origin = {44, 24}, extent = {{-13, 13}, {13, -13}}, textString = "T_FN2"), Rectangle(origin = {11, -10}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-15, 50}, {15, -50}}), Text(origin = {9, 68}, extent = {{-47, 12}, {47, -12}}, textString = "CoreRegion"), Rectangle(origin = {-15.43, -10.01}, lineThickness = 2, extent = {{-75.81, 90.1}, {75.81, -90.1}}), Text(origin = {-71, -45}, extent = {{-6, 8}, {6, -8}}, textString = "FF"), Text(origin = {-71, -5}, extent = {{-6, 8}, {6, -8}}, textString = "DH"), Rectangle(origin = {-15.43, -10.01}, lineThickness = 2, extent = {{-75.81, 90.1}, {75.81, -90.1}}), Rectangle(origin = {-15.43, -10.01}, lineThickness = 2, extent = {{-75.81, 90.1}, {75.81, -90.1}}), Rectangle(origin = {-34, -40}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}})}));
  end MSRRCoreRegion;
  
  model HeatExchanger
    parameter SMD_MSR_Modelica.Units.Mass m_PN1;
    parameter SMD_MSR_Modelica.Units.Mass m_PN2;
    parameter SMD_MSR_Modelica.Units.Mass m_PN3;
    parameter SMD_MSR_Modelica.Units.Mass m_PN4;
    parameter SMD_MSR_Modelica.Units.Mass m_TN1;
    parameter SMD_MSR_Modelica.Units.Mass m_TN2;
    parameter SMD_MSR_Modelica.Units.Mass m_SN1;
    parameter SMD_MSR_Modelica.Units.Mass m_SN2;
    parameter SMD_MSR_Modelica.Units.Mass m_SN3;
    parameter SMD_MSR_Modelica.Units.Mass m_SN4;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_pFluid;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_Tube;
    parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_sFluid;
    parameter SMD_MSR_Modelica.Units.MassFlowRate m_dot_pFluidNom;
    parameter SMD_MSR_Modelica.Units.MassFlowRate m_dot_sFluidNom;
    parameter SMD_MSR_Modelica.Units.Convection hApNom;
    parameter SMD_MSR_Modelica.Units.Convection hAsNom;
    parameter SMD_MSR_Modelica.Units.Temperature T_PN1_initial;
    parameter SMD_MSR_Modelica.Units.Temperature T_PN2_initial;
    parameter SMD_MSR_Modelica.Units.Temperature T_PN3_initial;
    parameter SMD_MSR_Modelica.Units.Temperature T_PN4_initial;
    parameter SMD_MSR_Modelica.Units.Temperature T_TN1_initial;
    parameter SMD_MSR_Modelica.Units.Temperature T_TN2_initial;
    parameter SMD_MSR_Modelica.Units.Temperature T_SN1_initial;
    parameter SMD_MSR_Modelica.Units.Temperature T_SN2_initial;
    parameter SMD_MSR_Modelica.Units.Temperature T_SN3_initial;
    parameter SMD_MSR_Modelica.Units.Temperature T_SN4_initial;
    parameter SMD_MSR_Modelica.Units.Conductivity KF;
    parameter SMD_MSR_Modelica.Units.Conductivity KS;
    parameter SMD_MSR_Modelica.Units.Area AcShell;
    parameter SMD_MSR_Modelica.Units.Area AcTube;
    parameter SMD_MSR_Modelica.Units.Area ArShell;
    parameter SMD_MSR_Modelica.Units.Length LF;
    parameter SMD_MSR_Modelica.Units.Length LS;
    parameter SMD_MSR_Modelica.Units.Emissivity e;
    parameter SMD_MSR_Modelica.Units.Temperature Tinf;
    constant SMD_MSR_Modelica.Units.StefanBoltzmannConstant SigSBK = 5.670374419E-14;
    SMD_MSR_Modelica.Units.MassFlowRate m_dot_pFluid;
    SMD_MSR_Modelica.Units.MassFlowRate m_dot_sFluid;
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
    input SMD_MSR_Modelica.PortsConnectors.Temp_In T_in_pFluid annotation(
      Placement(visible = true, transformation(origin = {-100, 32}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    input SMD_MSR_Modelica.PortsConnectors.Temp_In T_in_sFluid annotation(
      Placement(visible = true, transformation(origin = {163, -31}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {160, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    output SMD_MSR_Modelica.PortsConnectors.Temp_Out T_out_sFluid annotation(
      Placement(visible = true, transformation(origin = {-100, -38}, extent = {{-14, -14}, {14, 14}}, rotation = 0), iconTransformation(origin = {-100, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    output SMD_MSR_Modelica.PortsConnectors.Temp_Out T_out_pFluid annotation(
      Placement(visible = true, transformation(origin = {158, 28}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {160, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In primaryFF annotation(
      Placement(visible = true, transformation(origin = {-69, 49}, extent = {{-9, -9}, {9, 9}}, rotation = 0), iconTransformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In secondaryFF annotation(
      Placement(visible = true, transformation(origin = {123, -47}, extent = {{-11, -11}, {11, 11}}, rotation = 0), iconTransformation(origin = {130, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    input SMD_MSR_Modelica.PortsConnectors.volSpecificHeat_in P_decay annotation(
      Placement(visible = true, transformation(origin = {30, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {30, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  initial equation
    T_PN1 = T_PN1_initial;
    T_PN2 = T_PN2_initial;
    T_PN3 = T_PN3_initial;
    T_out_pFluid.T = T_PN4_initial;
    T_TN1 = T_TN1_initial;
    T_TN2 = T_TN2_initial;
    T_SN1 = T_SN1_initial;
    T_SN2 = T_SN2_initial;
    T_SN3 = T_SN3_initial;
    T_out_sFluid.T = T_SN4_initial;
  equation
    m_dot_pFluid = m_dot_pFluidNom*primaryFF.FF;
    m_dot_sFluid = m_dot_sFluidNom*secondaryFF.FF;
    hApn = 0.07797*primaryFF.FF^(0.33);
    hAsn = (1-0.01)*hAsNom*secondaryFF.FF + 0.01*hAsNom;
    if primaryFF.FF > 0 then
      m_PN1*cP_pFluid*der(T_PN1) = m_dot_pFluid*cP_pFluid*(T_in_pFluid.T - T_PN1) - hApn*(T_PN1 - T_TN1) + P_decay.Q*m_PN1 + e*SigSBK*ArShell*((Tinf + 273.15)^4 - (T_PN1 + 273.15)^4) + P_decay.Q*m_PN1;
      m_PN2*cP_pFluid*der(T_PN2) = m_dot_pFluid*cP_pFluid*(T_PN1 - T_PN2) - hApn*(T_PN1 - T_TN1) + P_decay.Q*m_PN2 + e*SigSBK*ArShell*((Tinf + 273.15)^4 - (T_PN2 + 273.15)^4) + P_decay.Q*m_PN2;
      m_PN3*cP_pFluid*der(T_PN3) = m_dot_pFluid*cP_pFluid*(T_PN2 - T_PN3) - hApn*(T_PN3 - T_TN2) + P_decay.Q*m_PN3 + e*SigSBK*ArShell*((Tinf + 273.15)^4 - (T_PN3 + 273.15)^4) + P_decay.Q*m_PN3;
      m_PN4*cP_pFluid*der(T_out_pFluid.T) = m_dot_pFluid*cP_pFluid*(T_PN3 - T_out_pFluid.T) - hApn*(T_PN3 - T_TN2) + P_decay.Q*m_PN4 + e*SigSBK*ArShell*((Tinf + 273.15)^4 - (T_out_pFluid.T + 273.15)^4) + P_decay.Q*m_PN4;
    else
      m_PN1*cP_pFluid*der(T_PN1) = ((KF*AcShell)/LF)*(T_in_pFluid.T - T_PN1) + e*SigSBK*ArShell*((Tinf + 273.15)^4 - (T_PN1 + 273.15)^4) + P_decay.Q*m_PN1;
      m_PN2*cP_pFluid*der(T_PN2) = ((KF*AcShell)/LF)*(T_PN1 - T_PN2) + e*SigSBK*ArShell*((Tinf + 273.15)^4 - (T_PN2 + 273.15)^4) + P_decay.Q*m_PN2;
      m_PN3*cP_pFluid*der(T_PN3) = ((KF*AcShell)/LF)*(T_PN2 - T_PN3) + e*SigSBK*ArShell*((Tinf + 273.15)^4 - (T_PN3 + 273.15)^4) + P_decay.Q*m_PN3;
      m_PN4*cP_pFluid*der(T_out_pFluid.T) = ((KF*AcShell)/LF)*(T_PN3 - T_out_pFluid.T) + e*SigSBK*ArShell*((Tinf + 273.15)^4 - (T_out_pFluid.T + 273.15)^4) + P_decay.Q*m_PN4;
    end if;
    if secondaryFF.FF > 0 then
      m_SN1*cP_sFluid*der(T_SN1) = m_dot_sFluid*cP_sFluid*(T_in_sFluid.T - T_SN1) + hAsn*(T_TN2 - T_SN1);
      m_SN2*cP_sFluid*der(T_SN2) = m_dot_sFluid*cP_sFluid*(T_SN1 - T_SN2) + hAsn*(T_TN2 - T_SN1);
      m_SN3*cP_sFluid*der(T_SN3) = m_dot_sFluid*cP_sFluid*(T_SN2 - T_SN3) + hAsn*(T_TN1 - T_SN3);
      m_SN4*cP_sFluid*der(T_out_sFluid.T) = m_dot_sFluid*cP_sFluid*(T_SN3 - T_out_sFluid.T) + hAsn*(T_TN1 - T_SN3);
    else
      m_SN1*cP_sFluid*der(T_SN1) = ((KS*AcTube)/LS)*(T_in_sFluid.T - T_SN1);
      m_SN2*cP_sFluid*der(T_SN2) = ((KS*AcTube)/LS)*(T_SN1 - T_SN2);
      m_SN3*cP_sFluid*der(T_SN3) = ((KS*AcTube)/LS)*(T_SN2 - T_SN3);
      m_SN4*cP_sFluid*der(T_out_sFluid.T) = ((KS*AcTube)/LS)*(T_SN3 - T_out_sFluid.T);
  
    end if;
    m_TN1*cP_Tube*der(T_TN1) = 2*hApn*(T_PN1 - T_TN1) - 2*hAsn*(T_TN1 - T_SN3);
    m_TN2*cP_Tube*der(T_TN2) = 2*hApn*(T_PN3 - T_TN2) - 2*hAsn*(T_TN2 - T_SN1);
    annotation(
      Diagram(graphics = {Rectangle(origin = {30.2864, 0.0340212}, extent = {{-149.242, 59.993}, {149.242, -59.993}}), Rectangle(origin = {-60, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {120, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {0, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-59.95, -29.61}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.63, -9.61}, {19.63, 9.61}}), Rectangle(origin = {-0.12, -29.82}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.88, 9.82}, {19.88, -9.82}}), Rectangle(origin = {59.8, -29.85}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.8, 9.85}, {19.8, -9.85}}), Rectangle(origin = {119.85, -29.89}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.85, 9.89}, {19.85, -9.89}}), Rectangle(origin = {89.92, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {60, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-30.08, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Text(origin = {150, 49}, extent = {{-28, 9}, {28, -9}}, textString = "HX")}, coordinateSystem(extent = {{-120, 60}, {180, -60}})),
      Icon(graphics = {Rectangle(origin = {30.29, 0.03}, lineThickness = 2, extent = {{-149.24, 59.99}, {149.24, -59.99}}), Rectangle(origin = {-30.08, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {60, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {89.92, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {-0.12, -29.82}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.88, 9.82}, {19.88, -9.82}}), Rectangle(origin = {0, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-60, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {120, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {59.8, -29.85}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.8, 9.85}, {19.8, -9.85}}), Rectangle(origin = {119.85, -29.89}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.85, 9.89}, {19.85, -9.89}}), Rectangle(origin = {-59.95, -29.61}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.63, -9.61}, {19.63, 9.61}}), Text(origin = {150, 49}, extent = {{-28, 9}, {28, -9}}, textString = "HX")}, coordinateSystem(extent = {{-120, 60}, {180, -60}})));
  end HeatExchanger;
  
  model MSRRuhx
  MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
      Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
  
  SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
      Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
  
  MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
      Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
  
  SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
      Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
  
  SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
      Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
      Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
      Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
      Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
  
  SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
      Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
  
  SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
      Placement(visible = true, transformation(origin = {-5891, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
      Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
      Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
      Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
  
  SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
      Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
  SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
      Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
  SMD_MSR_Modelica.InputSignals.Stepper stepper2(amplitude = {37}, numSteps = 1, stepTime = {0})  annotation(
      Placement(visible = true, transformation(origin = {-5172, 1754}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
      Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
  SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {100000})  annotation(
      Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
      Placement(visible = true, transformation(origin = {-5094, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0, primaryPumpK = 50, tripPrimaryPump = 100000)  annotation(
      Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
  equation
    connect(pke.n_population, powerBlock.nPop) annotation(
      Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
    connect(pke.n_population, decayHeat.nPop) annotation(
      Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
    connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
      Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
    connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
      Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
    connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
      Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
    connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
      Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
    connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
      Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
    connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
      Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
    connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
      Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
    connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
      Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
    connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
      Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
    connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
      Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
    connect(reactivityFeedback.feedback, pke.feedback) annotation(
      Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
    connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
      Line(points = {{-5892, 2092}, {-5808, 2092}, {-5808, 2052}, {-5810, 2052}}, color = {191, 64, 188}));
    connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
      Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
    connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
      Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
    connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
      Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
    connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
      Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
    connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
      Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
    connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
      Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
    connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
      Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
    connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
      Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
    connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
      Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
    connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
      Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
    connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
      Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
    connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
      Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
    connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
      Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
    connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
      Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
    connect(stepper.step, pke.ReactivityIn) annotation(
      Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
    connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
      Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5084, 1878}}, color = {204, 0, 0}));
    connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
      Line(points = {{-5125, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
    connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
      Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1899}, {-5125, 1899}}, color = {245, 121, 0}));
    connect(stepper1.step, uhx.PowerDemand) annotation(
      Line(points = {{-5126, 1738}, {-5126, 1857}, {-5125, 1857}}));
    connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
      Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
      Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
      Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
      Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
      Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
      Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
      Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
    annotation(
      Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
  end MSRRuhx;
  
  model MSRRrad
  MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
      Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
  
  SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
      Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
  
  MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 (displayUnit = "km")= 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
      Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
  
  SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
      Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
  
  SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
      Placement(visible = true, transformation(origin = {-5903, 1999}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
      Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
      Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
      Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
  
  SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
      Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
  
  SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
      Placement(visible = true, transformation(origin = {-5893, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
      Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
      Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
  
  SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
      Placement(visible = true, transformation(origin = {-5290, 1892}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
  
  SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
      Placement(visible = true, transformation(origin = {-5176, 1968}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
  SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
      Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
  SMD_MSR_Modelica.InputSignals.Stepper stepper2(amplitude = {37}, numSteps = 1, stepTime = {0})  annotation(
      Placement(visible = true, transformation(origin = {-5190, 1746}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
      Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0, primaryPumpK = 1, tripPrimaryPump = 1000000)  annotation(
      Placement(visible = true, transformation(origin = {-5560, 1726}, extent = {{-36, -36}, {36, 36}}, rotation = 0)));
  SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {100000})  annotation(
      Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.Radiator radiator(Cp_air = 1.0085E-3, Cp_salt = 2.3900E-3, Tair_0 = 79, Tp_0 = 503.9177, hA_radiatorNom = 0.0026, m_airN = 0.3949*1.1237, m_dot_air = 11.8044, m_dot_saltNom = 53.39, m_saltN = (0.03473*1767.1)) annotation(
      Placement(visible = true, transformation(origin = {-5126, 1864}, extent = {{-90, -90}, {90, 90}}, rotation = 0)));
  equation
    connect(pke.n_population, powerBlock.nPop) annotation(
      Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
    connect(pke.n_population, decayHeat.nPop) annotation(
      Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2020}, {-5898, 2020}}, color = {20, 36, 248}));
    connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
      Line(points = {{-5898, 1989}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
    connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
      Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
    connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
      Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
    connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
      Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
    connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
      Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
    connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
      Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
    connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
      Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
    connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
      Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
    connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
      Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
    connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
      Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
    connect(reactivityFeedback.feedback, pke.feedback) annotation(
      Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
    connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
      Line(points = {{-5894, 2092}, {-5810, 2092}, {-5810, 2052}}, color = {191, 64, 188}));
    connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
      Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
    connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
      Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
    connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
      Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
    connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
      Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
    connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
      Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
    connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
      Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
    connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
      Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
    connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
      Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
    connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
      Line(points = {{-5308, 1904}, {-5352, 1904}, {-5352, 2014}}, color = {204, 0, 0}));
    connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
      Line(points = {{-5178, 1980}, {-5176, 1980}, {-5176, 2064}, {-5087, 2064}, {-5087, 2028}}, color = {220, 138, 221}));
    connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
      Line(points = {{-5178, 1980}, {-5308, 1980}, {-5308, 1920}}, color = {220, 138, 221}));
    connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
      Line(points = {{-5268, 1688}, {-5308, 1688}, {-5308, 1890}}, color = {245, 121, 0}));
    connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
      Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
    connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
      Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
      Line(points = {{-5560, 1704}, {-5790, 1704}, {-5790, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
      Line(points = {{-5560, 1704}, {-5513, 1704}, {-5513, 2023}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
      Line(points = {{-5560, 1704}, {-5460, 1704}, {-5460, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
      Line(points = {{-5560, 1704}, {-5494, 1704}, {-5494, 1980}, {-5590, 1980}, {-5590, 2026}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
      Line(points = {{-5560, 1704}, {-5378, 1704}, {-5378, 1770}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
      Line(points = {{-5560, 1704}, {-5456, 1704}, {-5456, 2084}, {-5340, 2084}, {-5340, 2046}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
      Line(points = {{-5560, 1704}, {-5532, 1704}, {-5532, 2078}, {-5794, 2078}, {-5794, 2052}}, color = {245, 121, 0}));
    connect(stepper.step, pke.ReactivityIn) annotation(
      Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
  connect(stepper2.step, radiator.airTempIn) annotation(
      Line(points = {{-5189, 1752}, {-5190, 1752}, {-5190, 1855}, {-5187, 1855}}));
  connect(stepper1.step, radiator.airFlowFrac) annotation(
      Line(points = {{-5126, 1738}, {-5126, 1792}, {-5151, 1792}, {-5151, 1855}}));
  connect(pipe3.PiTemp_IN, radiator.RadTempOut) annotation(
      Line(points = {{-5086, 2014}, {-5034, 2014}, {-5034, 1882}, {-5101, 1882}}, color = {204, 0, 0}));
  connect(radiator.RadTempIN, pipe4.PiTemp_Out) annotation(
      Line(points = {{-5188, 1882}, {-5242, 1882}, {-5242, 1904}, {-5278, 1904}}, color = {204, 0, 0}));
  connect(secondaryPump.secondaryFlowFrac, radiator.secondaryFlowFrac) annotation(
      Line(points = {{-5268, 1688}, {-5230, 1688}, {-5230, 1900}, {-5188, 1900}}, color = {245, 121, 0}));
    annotation(
      Diagram(coordinateSystem(extent = {{-5940, 2120}, {-5080, 1640}})));
  end MSRRrad;

  package MSRRsteps
    model MSRRuhxHalfD
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5893, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper2(amplitude = {37}, numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5172, 1754}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {294.5347}, numSteps = 1, stepTime = {4000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0, primaryPumpK = 50, tripPrimaryPump = 100000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5894, 2092}, {-5810, 2092}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
    connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRuhxHalfD;
    
    model MSRRuhxOneD
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5893, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper2(amplitude = {37}, numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5172, 1754}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {589.0693}, numSteps = 1, stepTime = {4000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0, primaryPumpK = 50, tripPrimaryPump = 100000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5894, 2092}, {-5810, 2092}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
    connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRuhxOneD;
    
    model MSRRuhxTwoD
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5893, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper2(amplitude = {37}, numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5172, 1754}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {1178.1}, numSteps = 1, stepTime = {4000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0, primaryPumpK = 50, tripPrimaryPump = 100000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5894, 2092}, {-5810, 2092}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
    connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRuhxTwoD;
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false)));
  end MSRRsteps;

  package MSRRramp
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false)));
  end MSRRramp;

  package MSRRstepVarFlow
    model MSRR13Flow
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5893, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper2(amplitude = {37}, numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5172, 1754}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {589.0693}, numSteps = 1, stepTime = {6000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 1/3, primaryPumpK = 50, tripPrimaryPump = 1000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5894, 2092}, {-5810, 2092}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRR13Flow;
    
    model MSRR23Flow
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5893, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper2(amplitude = {37}, numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5172, 1754}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {589.0693}, numSteps = 1, stepTime = {6000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 2/3, primaryPumpK = 50, tripPrimaryPump = 1000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5894, 2092}, {-5810, 2092}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRR23Flow;
    
    model MSRR1Flow
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5893, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper2(amplitude = {37}, numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5172, 1754}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {589.0693}, numSteps = 1, stepTime = {6000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0, primaryPumpK = 50, tripPrimaryPump = 100000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5894, 2092}, {-5810, 2092}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRR1Flow;
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false)));
  end MSRRstepVarFlow;

  package MSRRstartUp
    model MSRRstartUpUHX
    MSRR.HeatExchanger heatExchanger(AcShell = 0.0072, AcTube = 0.0045, ArShell = 0.4132*12.5185,KF = 1E-6, KS = 1.1E-6, LF = 12.5185, LS = 12.5185*2, T_PN1_initial = 570, T_PN2_initial = 570, T_PN3_initial = 570, T_PN4_initial = 570, T_SN1_initial = 570, T_SN2_initial = 570, T_SN3_initial = 570, T_SN4_initial = 570, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 570,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0.5, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = false, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 0, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51/2, LF2 = 1.51/2, OutterRegion = true,TF1_0 = 570, TF2_0 = 570, TG_0 = 570, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.5, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 0.04, Ar = 0.7081*1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 1E-6, L = 1, T_0 = 570, Tinf = 570, e = 0.5, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5488.1, 2032.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 0.0027, Ar = 0.1853*1.0703,CpFluid = 2009.66E-6, KF = 1E-6, L = 1.0703, T_0 = 570, Tinf = 570, e = 0.5, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0029)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 0.0027, Ar = 0.1853*2.1407,CpFluid = 2009.66E-6, KF = 1E-6, L = 2.1407, T_0 = 570, Tinf = 570, e = 0.5, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0058)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 0.0027, Ar = 0.1853*3.2110, CpFluid = 2009.66E-6, KF = 1E-6, L = 3.2110, T_0 = 570, Tinf = 570, e = 0.5, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0088)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 0.0027, Ar = 0.1853*2.1407,CpFluid = 2.3900E-3, KF = 1.1E-6, L = 2.1407, T_0 = 570, Tinf = 570, e = 0.5, mDotNom = 53.39, mFracNode = 0.1, mPi = 1767.1*0.0088)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 0.0027, Ar = 0.1853*2.1407,CpFluid = 2.3900E-3, KF = 1.1E-6, L = 2.1407, T_0 = 570, Tinf = 570, e = 0.5, mDotNom = 53.39,mFracNode = 0.1, mPi = 1767.1*0.0088)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {0},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper2(amplitude = {37}, numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5172, 1754}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {-3600, -3600, -1072.9, -1033.6, -887.5173, -745.4337, -602.7790, -459.4196, -314.6013, -166.0349, -98.5330, -54.0801, -47.091, -35.2460, -19.2281, -9.3705, -5.2325, -2.536, 6.6989, 13.826, 24.94}, numSteps = 21, stepTime = {0, 3600, 7200, 10800, 14400, 18000, 21600, 25200, 28800, 32400, 36000, 39600, 46800, 57600, 61200, 68400, 79200, 82800, 86400, 97200, 100800})  annotation(
        Placement(visible = true, transformation(origin = {-6004, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, Tp_0 = 570, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.Pump pump(freeConvFF = 0, numRampUp = 3, rampUpK = {50, 50, 50}, rampUpTime = {50400, 72000, 90000}, rampUpTo = {0.25, 0.50, 1.0}, tripK = 50, tripTime = 75600)  annotation(
        Placement(visible = true, transformation(origin = {-5619, 1731}, extent = {{-49, -49}, {49, 49}}, rotation = 0)));
  MSRR.MSRRstartUp.ConstantNeutronSource constantNeutronSource(S = 1E6)  annotation(
        Placement(visible = true, transformation(origin = {-5989, 2084}, extent = {{-27, -22}, {27, 22}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5512, 2040}, {-5511, 2040}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5473, 2040}, {-5451.5, 2040}, {-5451.5, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2055}, {-5511, 2055}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6003, 2001}, {-5996, 2001}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(pump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5618, 1702}, {-5378, 1702}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(pump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5618, 1702}, {-5428, 1702}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(pump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5618, 1702}, {-5511, 1702}, {-5511, 2025}}, color = {245, 121, 0}));
      connect(pump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5618, 1702}, {-5606, 1702}, {-5606, 2010}, {-5590, 2010}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(pump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5618, 1702}, {-5764, 1702}, {-5764, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(pump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5618, 1702}, {-5764, 1702}, {-5764, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(pump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5618, 1702}, {-5366, 1702}, {-5366, 2078}, {-5340, 2078}, {-5340, 2046}}, color = {245, 121, 0}));
  connect(constantNeutronSource.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5990, 2092}, {-5810, 2092}, {-5810, 2052}}, color = {191, 64, 188}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRstartUpUHX;
    
    model ConstantNeutronSource
      parameter SMD_MSR_Modelica.Units.NeutronEmissionRate S;
      output SMD_MSR_Modelica.PortsConnectors.NeutronEmsRateOut neutronEmsRateOut annotation(
        Placement(visible = true, transformation(origin = {-2, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-2, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      when time >= 3600 then
        neutronEmsRateOut.nDot = S;
      elsewhen time >=  46900 then
        neutronEmsRateOut.nDot = 0;
      elsewhen time >= 54000 then
        neutronEmsRateOut.nDot = S;
      elsewhen time >=  68500 then
        neutronEmsRateOut.nDot = 0;
      elsewhen time >= 75600 then
        neutronEmsRateOut.nDot = S;
      elsewhen time >=  86500 then
        neutronEmsRateOut.nDot = 0;
      elsewhen time >= 93600 then
        neutronEmsRateOut.nDot = S;
      elsewhen time >=  100900 then
        neutronEmsRateOut.nDot = 0;        
      end when;
      annotation(
        Diagram(graphics = {Rectangle(lineColor = {191, 64, 188}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Rectangle(lineColor = {191, 64, 188}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end ConstantNeutronSource;
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false)));
  end MSRRstartUp;
  
  package MSRRmaxAccident
    model MSRRmaxAccidentUHX
    MSRR.HeatExchanger heatExchanger(AcShell = 0.0072, AcTube = 0.0045, ArShell = 0.4132*12.5185,KF = 1E-6, KS = 1.1E-6, LF = 12.5185, LS = 12.5185*2, T_PN1_initial = 580, T_PN2_initial = 580, T_PN3_initial = 550, T_PN4_initial = 550, T_SN1_initial = 500, T_SN2_initial = 500, T_SN3_initial = 508, T_SN4_initial = 508, T_TN1_initial = 520, T_TN2_initial = 520, Tinf = 550,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0.1, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51/2, LF2 = 1.51/2, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 550, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.1, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 0.04, Ar = 0.7081*1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 1E-6, L = 1, T_0 = 570, Tinf = 550, e = 0.1, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5488.1, 2032.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 0.0027, Ar = 0.1853*1.0703,CpFluid = 2009.66E-6, KF = 1E-6, L = 1.0703, T_0 = 570, Tinf = 550, e = 0.1, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0029)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 0.0027, Ar = 0.1853*2.1407,CpFluid = 2009.66E-6, KF = 1E-6, L = 2.1407, T_0 = 570, Tinf = 550, e = 0.1, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0058)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 0.0027, Ar = 0.1853*3.2110, CpFluid = 2009.66E-6, KF = 1E-6, L = 3.2110, T_0 = 570, Tinf = 550, e = 0.1, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0088)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 0.0027, Ar = 0.1853*2.1407,CpFluid = 2.3900E-3, KF = 1.1E-6, L = 2.1407, T_0 = 508, Tinf = 550, e = 0.1, mDotNom = 53.39, mFracNode = 0.1, mPi = 1767.1*0.0088)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 0.0027, Ar = 0.1853*2.1407,CpFluid = 2.3900E-3, KF = 1.1E-6, L = 2.1407, T_0 = 500, Tinf = 550, e = 0.1, mDotNom = 53.39,mFracNode = 0.1, mPi = 1767.1*0.0088)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1, 0},numSteps = 2, stepTime = {0, 3600})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper2(amplitude = {37}, numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5172, 1754}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.01, secondaryPumpK = 50, tripSecondaryPump = 3600)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-6004, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, Tp_0 = 508, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
  MSRR.MSRRstartUp.ConstantNeutronSource constantNeutronSource(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5989, 2084}, extent = {{-27, -22}, {27, 22}}, rotation = 0)));
  SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.01, primaryPumpK = 50, tripPrimaryPump = 3600)  annotation(
        Placement(visible = true, transformation(origin = {-5578, 1722}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5512, 2040}, {-5511, 2040}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5473, 2040}, {-5451.5, 2040}, {-5451.5, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2055}, {-5511, 2055}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6003, 2001}, {-5996, 2001}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(constantNeutronSource.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5990, 2092}, {-5810, 2092}, {-5810, 2052}}, color = {191, 64, 188}));
  connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5578, 1698}, {-5770, 1698}, {-5770, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
  connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5578, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
  connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5578, 1698}, {-5510, 1698}, {-5510, 2024}}, color = {245, 121, 0}));
  connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5578, 1698}, {-5618, 1698}, {-5618, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
  connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5578, 1698}, {-5428, 1698}, {-5428, 2026}}, color = {245, 121, 0}));
  connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5578, 1698}, {-5452, 1698}, {-5452, 2078}, {-5340, 2078}, {-5340, 2046}}, color = {245, 121, 0}));
  connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5578, 1698}, {-5770, 1698}, {-5770, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRmaxAccidentUHX;
    
    model ConstantNeutronSource
      parameter SMD_MSR_Modelica.Units.NeutronEmissionRate S;
      output SMD_MSR_Modelica.PortsConnectors.NeutronEmsRateOut neutronEmsRateOut annotation(
        Placement(visible = true, transformation(origin = {-2, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-2, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      when time >= 3600 then
        neutronEmsRateOut.nDot = S;
      elsewhen time >=  46900 then
        neutronEmsRateOut.nDot = 0;
      elsewhen time >= 54000 then
        neutronEmsRateOut.nDot = S;
      elsewhen time >=  68500 then
        neutronEmsRateOut.nDot = 0;
      elsewhen time >= 75600 then
        neutronEmsRateOut.nDot = S;
      elsewhen time >=  86500 then
        neutronEmsRateOut.nDot = 0;
      elsewhen time >= 93600 then
        neutronEmsRateOut.nDot = S;
      elsewhen time >=  100900 then
        neutronEmsRateOut.nDot = 0;        
      end when;
      annotation(
        Diagram(graphics = {Rectangle(lineColor = {191, 64, 188}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Rectangle(lineColor = {191, 64, 188}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end ConstantNeutronSource;
    
    model MSRRmaxAccidentUHXScram
    MSRR.HeatExchanger heatExchanger(AcShell = 0.0072, AcTube = 0.0045, ArShell = 0.4132*12.5185,KF = 1E-6, KS = 1.1E-6, LF = 12.5185, LS = 12.5185*2, T_PN1_initial = 580, T_PN2_initial = 580, T_PN3_initial = 550, T_PN4_initial = 550, T_SN1_initial = 500, T_SN2_initial = 500, T_SN3_initial = 508, T_SN4_initial = 508, T_TN1_initial = 520, T_TN2_initial = 520, Tinf = 550,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0.1, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51/2, LF2 = 1.51/2, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 550, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.1, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 0.04, Ar = 0.7081*1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 1E-6, L = 1, T_0 = 570, Tinf = 550, e = 0.1, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5488.1, 2032.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 0.0027, Ar = 0.1853*1.0703,CpFluid = 2009.66E-6, KF = 1E-6, L = 1.0703, T_0 = 570, Tinf = 550, e = 0.1, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0029)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 0.0027, Ar = 0.1853*2.1407,CpFluid = 2009.66E-6, KF = 1E-6, L = 2.1407, T_0 = 570, Tinf = 550, e = 0.1, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0058)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 0.0027, Ar = 0.1853*3.2110, CpFluid = 2009.66E-6, KF = 1E-6, L = 3.2110, T_0 = 570, Tinf = 550, e = 0.1, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0088)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 0.0027, Ar = 0.1853*2.1407,CpFluid = 2.3900E-3, KF = 1.1E-6, L = 2.1407, T_0 = 508, Tinf = 550, e = 0.1, mDotNom = 53.39, mFracNode = 0.1, mPi = 1767.1*0.0088)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 0.0027, Ar = 0.1853*2.1407,CpFluid = 2.3900E-3, KF = 1.1E-6, L = 2.1407, T_0 = 500, Tinf = 550, e = 0.1, mDotNom = 53.39,mFracNode = 0.1, mPi = 1767.1*0.0088)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1, 0},numSteps = 2, stepTime = {0, 2000})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper2(amplitude = {37}, numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5172, 1754}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.01, secondaryPumpK = 50, tripSecondaryPump = 2000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {-1000}, numSteps = 1, stepTime = {2000})  annotation(
        Placement(visible = true, transformation(origin = {-6004, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, Tp_0 = 508, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    MSRR.MSRRstartUp.ConstantNeutronSource constantNeutronSource(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5989, 2084}, extent = {{-27, -22}, {27, 22}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.01, primaryPumpK = 50, tripPrimaryPump = 2000)  annotation(
        Placement(visible = true, transformation(origin = {-5578, 1722}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5512, 2040}, {-5511, 2040}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5473, 2040}, {-5451.5, 2040}, {-5451.5, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2055}, {-5511, 2055}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6003, 2001}, {-5996, 2001}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(constantNeutronSource.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5990, 2092}, {-5810, 2092}, {-5810, 2052}}, color = {191, 64, 188}));
    connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5578, 1698}, {-5770, 1698}, {-5770, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5578, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5578, 1698}, {-5510, 1698}, {-5510, 2024}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5578, 1698}, {-5618, 1698}, {-5618, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5578, 1698}, {-5428, 1698}, {-5428, 2026}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5578, 1698}, {-5452, 1698}, {-5452, 2078}, {-5340, 2078}, {-5340, 2046}}, color = {245, 121, 0}));
    connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5578, 1698}, {-5770, 1698}, {-5770, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRmaxAccidentUHXScram;
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false)));
  end MSRRmaxAccident;

  package MSRRpumpTrip
    
    model MSRRpump1
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5891, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {100000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.01, primaryPumpK = 1, tripPrimaryPump = 2000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5892, 2092}, {-5808, 2092}, {-5808, 2052}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRpump1;
    
    model MSRRpump5
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5891, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {100000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.01, primaryPumpK = 5, tripPrimaryPump = 2000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5892, 2092}, {-5808, 2092}, {-5808, 2052}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRpump5;
    
    model MSRRpump10
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5891, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {100000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.01, primaryPumpK = 10, tripPrimaryPump = 2000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5892, 2092}, {-5808, 2092}, {-5808, 2052}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRpump10;
    
    model MSRRpump15
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5891, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {100000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.01, primaryPumpK = 15, tripPrimaryPump = 2000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5892, 2092}, {-5808, 2092}, {-5808, 2052}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRpump15;
    
    model MSRRpump20
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5891, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {100000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.01, primaryPumpK = 20, tripPrimaryPump = 2000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5892, 2092}, {-5808, 2092}, {-5808, 2052}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRpump20;
    model MSRRpump25
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5891, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {100000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.01, primaryPumpK = 25, tripPrimaryPump = 2000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5892, 2092}, {-5808, 2092}, {-5808, 2052}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRpump25;
    
    model MSRRpump30
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5918, 2000}, extent = {{-66, -66}, {66, 66}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5891, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {100000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.01, primaryPumpK = 30, tripPrimaryPump = 2000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2027}, {-5912, 2027}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5911, 1987}, {-5911, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5892, 2092}, {-5808, 2092}, {-5808, 2052}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRpump30;
    
    model MSRRpump40
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5891, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {100000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.01, primaryPumpK = 40, tripPrimaryPump = 2000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5892, 2092}, {-5808, 2092}, {-5808, 2052}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRpump40;
    
    model MSRRpump50
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5891, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {100000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.01, primaryPumpK = 50, tripPrimaryPump = 2000)  annotation(
        Placement(visible = true, transformation(origin = {-5595, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5892, 2092}, {-5808, 2092}, {-5808, 2052}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5595, 1699}, {-5774, 1699}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5595, 1699}, {-5774, 1699}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5595, 1699}, {-5595, 1770}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5595, 1699}, {-5595, 2022}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5595, 1699}, {-5614, 1699}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5595, 1699}, {-5458, 1699}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5595, 1699}, {-5595, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRpump50;
    
    model MSRRpump75
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5891, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {100000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.01, primaryPumpK = 75, tripPrimaryPump = 2000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5892, 2092}, {-5808, 2092}, {-5808, 2052}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRpump75;
    
    model MSRRpump100
    MSRR.HeatExchanger heatExchanger(AcShell = 1, AcTube = 1, ArShell = 1,KF = 0, KS = 0, LF = 1, LS = 1, T_PN1_initial = 580.4100, T_PN2_initial = 580.4100, T_PN3_initial = 560, T_PN4_initial = 560, T_SN1_initial = 508, T_SN2_initial = 508, T_SN3_initial = 500, T_SN4_initial = 500, T_TN1_initial = 570, T_TN2_initial = 570, Tinf = 100,cP_Tube = 4.5000E-4,cP_pFluid = 2009.66E-6, cP_sFluid = 2.3900E-3, e = 0, hApNom = 0.02592, hAsNom = 0.02624, m_PN1 = (0.0383*2079.4)/4, m_PN2 = (0.0383*2079.4)/4, m_PN3 = (0.0383*2079.4)/4, m_PN4 = (0.0383*2079.4)/4, m_SN1 = (0.03473*1767.1)/4, m_SN2 = (0.03473*1767.1)/4, m_SN3 = (0.03473*1767.1)/4, m_SN4 = (0.03473*1767.1)/4, m_TN1 = (0.03226*7980)/2, m_TN2 = (0.03226*7980)/2, m_dot_pFluidNom = 23.9000, m_dot_sFluidNom = 53.4) annotation(
        Placement(visible = true, transformation(origin = {-5312.4, 2026.2}, extent = {{-47.6, 23.8}, {71.4, 71.4}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, addRho0 = true, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023},lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, n_0 = 1, nominalTauCore = 34.8025, nominalTauLoop = 8.7006, numGroups = 6)  annotation(
        Placement(visible = true, transformation(origin = {-5802, 2044}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
    
    MSRR.MSRRCoreRegion fuelChannel(Ac = 0.2636, ArF1 = 53.0317/2, ArF2 = 53.0317/2, KF = 1E-6, LF1 = 1.51, LF2 = 1.51, OutterRegion = true,TF1_0 = 570, TF2_0 = 580.4100, TG_0 = 572.8112, Tinf = 570, cP_fuel = 2009.66E-6, cP_graphite = 1.773E-3, e = 0.01, hAnom = 0.0249, kFN1 = 0.4650, kFN2 = 0.4650, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, m_FN1 = 2079.4*(0.4/2), m_FN2 = 2079.4*(0.4/2), m_GN = 3.1644e+03, mdot_fuelNom = 23.9, regionFlowFrac = 1) annotation(
        Placement(visible = true, transformation(origin = {-5676, 1864}, extent = {{-70, 56}, {42, 182}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.PowerBlock powerBlock(P = 1, TotalFuelMass = 2079.4*0.55)  annotation(
        Placement(visible = true, transformation(origin = {-5815, 1955}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.DecayHeat decayHeat(DHYG1 = 9.635981959409105E-01*2.464783802008740E-03, DHYG2 = 3.560674858154914E-02*2.464783802008740E-03, DHYG3 = 7.950554775404400E-04*2.464783802008740E-03, DHlamG1 = 0.09453, DHlamG2 = 0.004420, DHlamG3 = 8.6098E-5) annotation(
        Placement(visible = true, transformation(origin = {-5903, 2005}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Ac = 1, Ar = 1,Cp_fluid = 2009.66E-6, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 20000000, DHRS_timeConstant = 50, KF = 0, L = 1, T_0 = 580, Tinf = 100, e = 0, m_DHRS = 2079.4*0.04, m_dotNom = 23.9000) annotation(
        Placement(visible = true, transformation(origin = {-5490.1, 2030.1}, extent = {{-29.9048, 29.9048}, {22.4286, 74.7619}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 1, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5573, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe1(Ac = 1, Ar = 1,CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 1, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5411, 2027}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
    
    SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 570, FuelTempSetPointNode2 = 580.4100, GrapTempSetPoint = 572.8112, IF1 = 0.5, IF2 = 0.5, IG = 1, a_F = -6.26E-5, a_G = -5.16E-5) annotation(
        Placement(visible = true, transformation(origin = {-5569, 1841}, extent = {{-45, -45}, {45, 45}}, rotation = 90)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSourceZero constantNeutronSourceZero(S = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5891, 2081}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe2(Ac = 1, Ar = 1, CpFluid = 2009.66E-6, KF = 1, L = 1, T_0 = 550, Tinf = 100, e = 0, mDotNom = 23.9, mFracNode = 0.1, mPi = 2079.4*0.0072)  annotation(
        Placement(visible = true, transformation(origin = {-5396, 1772}, extent = {{46, -46}, {-46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe3(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 508, Tinf = 100, e = 0, mDotNom = 53.39, mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5103, 2003}, extent = {{41, -41}, {-41, 41}}, rotation = 0)));
    
    SMD_MSR_Modelica.HeatTransport.Pipe pipe4(Ac = 1, Ar = 1,CpFluid = 2.3900E-3, KF = 1, L = 1, T_0 = 500, Tinf = 100, e = 0, mDotNom = 53.39,mFracNode = 0.1, mPi = 53.39*1)  annotation(
        Placement(visible = true, transformation(origin = {-5262, 1866}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
    
    SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
        Placement(visible = true, transformation(origin = {-5176, 2082}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper1(amplitude = {1},numSteps = 1, stepTime = {0})  annotation(
        Placement(visible = true, transformation(origin = {-5127, 1733}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.SecondaryPump secondaryPump(freeConvectionFF = 0.02, secondaryPumpK = 100, tripSecondaryPump = 2000000)  annotation(
        Placement(visible = true, transformation(origin = {-5268, 1710}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
    SMD_MSR_Modelica.InputSignals.Stepper stepper(amplitude = {0}, numSteps = 1, stepTime = {100000})  annotation(
        Placement(visible = true, transformation(origin = {-6000, 1996}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.UHX uhx(P = 1, cP_pFluidUHX = 2.3900E-3, m_PN1UHX = (0.03473*1767.1), m_dot_pFluidUHXnom = 53.4)  annotation(
        Placement(visible = true, transformation(origin = {-5116, 1878}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    SMD_MSR_Modelica.HeatTransport.PrimaryPump primaryPump(freeConvectionFF = 0.01, primaryPumpK = 100, tripPrimaryPump = 2000)  annotation(
        Placement(visible = true, transformation(origin = {-5597, 1727}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
    equation
      connect(pke.n_population, powerBlock.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5806.4, 1956.4}}, color = {20, 36, 248}));
      connect(pke.n_population, decayHeat.nPop) annotation(
        Line(points = {{-5806.4, 2026.4}, {-5847.9, 2026.4}, {-5847.9, 2026}, {-5898, 2026}}, color = {20, 36, 248}));
      connect(decayHeat.decayHeat_Out, powerBlock.decayNP) annotation(
        Line(points = {{-5898, 1995}, {-5898, 1956.8}, {-5832.9, 1956.8}}, color = {0, 225, 255}));
      connect(powerBlock.fissionPower, fuelChannel.fissionPower) annotation(
        Line(points = {{-5807, 1930}, {-5806, 1930}, {-5806, 1903}, {-5726, 1903}}, color = {129, 61, 156}));
      connect(powerBlock.decayPowerM, fuelChannel.decayHeat) annotation(
        Line(points = {{-5831, 1930}, {-5832, 1930}, {-5832, 1872}, {-5726, 1872}}, color = {220, 138, 221}));
      connect(pipe.PiTemp_Out, dhrs.DHRS_tempIN) annotation(
        Line(points = {{-5560, 2039}, {-5512, 2039}, {-5513, 2038}}, color = {204, 0, 0}));
      connect(fuelChannel.fuelNode2, pipe.PiTemp_IN) annotation(
        Line(points = {{-5648, 1894}, {-5624, 1894}, {-5624, 2038}, {-5590, 2038}}, color = {204, 0, 0}));
      connect(dhrs.DHRS_TempOUT, pipe1.PiTemp_IN) annotation(
        Line(points = {{-5475, 2038}, {-5428, 2038}}, color = {204, 0, 0}));
      connect(pipe1.PiTemp_Out, heatExchanger.T_in_pFluid) annotation(
        Line(points = {{-5398, 2038}, {-5352, 2038}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode1, fuelChannel.fuelNode1) annotation(
        Line(points = {{-5587, 1823}, {-5648, 1823}, {-5648, 1822}}, color = {204, 0, 0}));
      connect(reactivityFeedback.fuelNode2, fuelChannel.fuelNode2) annotation(
        Line(points = {{-5587, 1836.5}, {-5614, 1836.5}, {-5614, 1894}, {-5648, 1894}}, color = {204, 0, 0}));
      connect(reactivityFeedback.grapNode, fuelChannel.grapNode) annotation(
        Line(points = {{-5587, 1850}, {-5648, 1850}, {-5648, 1858}}, color = {204, 0, 0}));
      connect(reactivityFeedback.feedback, pke.feedback) annotation(
        Line(points = {{-5560, 1836}, {-5542, 1836}, {-5542, 1954}, {-5768, 1954}, {-5768, 2094}, {-5802, 2094}, {-5802, 2052}}, color = {78, 154, 6}));
      connect(constantNeutronSourceZero.neutronEmsRateOut, pke.S) annotation(
        Line(points = {{-5892, 2092}, {-5808, 2092}, {-5808, 2052}, {-5810, 2052}}, color = {191, 64, 188}));
      connect(heatExchanger.T_out_pFluid, pipe2.PiTemp_IN) annotation(
        Line(points = {{-5248, 2038}, {-5218, 2038}, {-5218, 1784}, {-5378, 1784}}, color = {204, 0, 0}));
      connect(pipe2.PiTemp_Out, fuelChannel.temp_In) annotation(
        Line(points = {{-5409, 1784}, {-5726, 1784}, {-5726, 1814}}, color = {204, 0, 0}));
      connect(powerBlock.decayPowerM, pipe.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5610, 1930}, {-5610, 2054}, {-5590, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, dhrs.pDecay) annotation(
        Line(points = {{-5832, 1930}, {-5530, 1930}, {-5530, 2053}, {-5513, 2053}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe1.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5452, 1930}, {-5452, 2054}, {-5428, 2054}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, heatExchanger.P_decay) annotation(
        Line(points = {{-5832, 1930}, {-5366, 1930}, {-5366, 2064}, {-5300, 2064}, {-5300, 2046}}, color = {220, 138, 221}));
      connect(powerBlock.decayPowerM, pipe2.PiDecay_Heat) annotation(
        Line(points = {{-5832, 1930}, {-5378, 1930}, {-5378, 1800}}, color = {220, 138, 221}));
      connect(heatExchanger.T_in_sFluid, pipe3.PiTemp_Out) annotation(
        Line(points = {{-5248, 2014}, {-5114, 2014}}, color = {204, 0, 0}));
      connect(pipe4.PiTemp_IN, heatExchanger.T_out_sFluid) annotation(
        Line(points = {{-5280, 1878}, {-5352, 1878}, {-5352, 2014}}, color = {204, 0, 0}));
      connect(constantDecayPower.decayHeat_Out, pipe3.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5087, 2094}, {-5087, 2028}}, color = {220, 138, 221}));
      connect(constantDecayPower.decayHeat_Out, pipe4.PiDecay_Heat) annotation(
        Line(points = {{-5177, 2094}, {-5177, 1978}, {-5280, 1978}, {-5280, 1894}}, color = {220, 138, 221}));
      connect(secondaryPump.secondaryFlowFrac, pipe4.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5268, 1686}, {-5280, 1686}, {-5280, 1864}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, heatExchanger.secondaryFF) annotation(
        Line(points = {{-5268, 1688}, {-5260, 1688}, {-5260, 2006}}, color = {245, 121, 0}));
      connect(secondaryPump.secondaryFlowFrac, pipe3.fuelFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5086, 1688}, {-5086, 2002}}, color = {245, 121, 0}));
      connect(stepper.step, pke.ReactivityIn) annotation(
        Line(points = {{-6000, 2002}, {-5996, 2002}, {-5996, 2052}, {-5820, 2052}}));
      connect(pipe3.PiTemp_IN, uhx.UHXtemp_Out) annotation(
        Line(points = {{-5086, 2014}, {-5014, 2014}, {-5014, 1878}, {-5106, 1878}}, color = {204, 0, 0}));
      connect(uhx.UHXtemp_In, pipe4.PiTemp_Out) annotation(
        Line(points = {{-5148, 1878}, {-5249, 1878}}, color = {204, 0, 0}));
      connect(secondaryPump.secondaryFlowFrac, uhx.secondaryFlowFrac) annotation(
        Line(points = {{-5268, 1688}, {-5192, 1688}, {-5192, 1898}, {-5148, 1898}}, color = {245, 121, 0}));
      connect(stepper1.step, uhx.PowerDemand) annotation(
        Line(points = {{-5126, 1738}, {-5148, 1738}, {-5148, 1858}}));
      connect(primaryPump.primaryFlowFrac, fuelChannel.fuelFlowFraction) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 1844}, {-5724, 1844}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5774, 1698}, {-5774, 2052}, {-5794, 2052}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe2.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5378, 1698}, {-5378, 1770}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, dhrs.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5512, 1698}, {-5512, 2022}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5614, 1698}, {-5614, 2026}, {-5590, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, pipe1.fuelFlowFrac) annotation(
        Line(points = {{-5596, 1698}, {-5458, 1698}, {-5458, 2026}, {-5428, 2026}}, color = {245, 121, 0}));
      connect(primaryPump.primaryFlowFrac, heatExchanger.primaryFF) annotation(
        Line(points = {{-5596, 1698}, {-5456, 1698}, {-5456, 2082}, {-5340, 2082}, {-5340, 2046}}, color = {245, 121, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-6020, 2120}, {-5000, 1660}})));
    end MSRRpump100;
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false)));
  end MSRRpumpTrip;
end MSRR;