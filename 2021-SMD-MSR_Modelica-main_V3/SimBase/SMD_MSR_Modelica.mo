package SMD_MSR_Modelica
  package HeatTransport
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
      hApn = (hApNom/4)*(0.8215*primaryFF.FF^6 - 4.108*primaryFF.FF^5 + 7.848*primaryFF.FF^4 - 7.165*primaryFF.FF^3 + 3.004*primaryFF.FF^2 + 0.5903*primaryFF.FF + 0.008537);
      hAsn = (hAsNom/4)*(0.8215*secondaryFF.FF^6 - 4.108*secondaryFF.FF^5 + 7.848*secondaryFF.FF^4 - 7.165*secondaryFF.FF^3 + 3.004*secondaryFF.FF^2 + 0.5903*secondaryFF.FF + 0.008537);
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

    model Radiator
      parameter SMD_MSR_Modelica.Units.MassFlowRate m_dot_air;
      parameter SMD_MSR_Modelica.Units.MassFlowRate m_dot_saltNom;
      parameter SMD_MSR_Modelica.Units.Convection hA_radiatorNom;
      parameter SMD_MSR_Modelica.Units.Mass m_airN;
      parameter SMD_MSR_Modelica.Units.Mass m_saltN;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity Cp_salt;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity Cp_air;
      parameter SMD_MSR_Modelica.Units.Temperature Tp_0;
      parameter SMD_MSR_Modelica.Units.Temperature Tair_0;
      constant SMD_MSR_Modelica.Units.Temperature realToTemp = 1;
      SMD_MSR_Modelica.Units.Convection hA_radiator;
      SMD_MSR_Modelica.Units.Temperature AirTempOut;
      SMD_MSR_Modelica.Units.MassFlowRate m_dot_salt;
      SMD_MSR_Modelica.Units.Temperature AirTempInRad;
      input SMD_MSR_Modelica.PortsConnectors.Temp_In RadTempIN annotation(
        Placement(visible = true, transformation(origin = {6, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out RadTempOut annotation(
        Placement(visible = true, transformation(origin = {-10, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {28, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In secondaryFlowFrac annotation(
        Placement(visible = true, transformation(origin = {30, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.RealIn airTempIn annotation(
        Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-68, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.RealIn airFlowFrac annotation(
        Placement(visible = true, transformation(origin = {-22, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-28, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      RadTempOut.T = Tp_0;
      AirTempOut = Tair_0;
    equation
      AirTempInRad = airTempIn.R*realToTemp;
      m_dot_salt = m_dot_saltNom*secondaryFlowFrac.FF;
      
      hA_radiator = hA_radiatorNom*(0.8215*secondaryFlowFrac.FF^6 - 4.108*secondaryFlowFrac.FF^5 + 7.848*secondaryFlowFrac.FF^4 - 7.165*secondaryFlowFrac.FF^3 + 3.004*secondaryFlowFrac.FF^2 + 0.5903*secondaryFlowFrac.FF + 0.008537);
      
      der(RadTempOut.T)*m_saltN*Cp_salt = m_dot_salt*Cp_salt*(RadTempIN.T - RadTempOut.T) - hA_radiator*(RadTempOut.T - AirTempOut);
      der(AirTempOut)*m_airN*Cp_air = airFlowFrac.R*m_dot_air*Cp_air*(AirTempInRad - AirTempOut) + hA_radiator*(RadTempOut.T - AirTempOut);
      annotation(
        Diagram(graphics = {Rectangle(origin = {-10, -10}, lineThickness = 2, extent = {{-50, 50}, {50, -50}}), Text(origin = {7, 17}, extent = {{-15, 5}, {15, -5}}, textString = "TpIN", fontSize = 20), Text(origin = {23, -44}, extent = {{-13, 8}, {13, -8}}, textString = "Rad"), Text(origin = {-9, -51}, extent = {{-15, 5}, {15, -5}}, textString = "TpOUT", fontSize = 20), Text(origin = {31, 17}, extent = {{-15, 5}, {15, -5}}, textString = "FFp", fontSize = 20), Text(origin = {-49, 17}, extent = {{-15, 5}, {15, -5}}, textString = "TsIN", fontSize = 20), Text(origin = {-21, 17}, extent = {{-15, 5}, {15, -5}}, textString = "FFs", fontSize = 20), Rectangle(origin = {-32, -10}, lineColor = {224, 27, 36}, fillColor = {224, 27, 36}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-18, 10}, {18, -10}}), Rectangle(origin = {12, -10}, lineColor = {143, 240, 164}, fillColor = {143, 240, 164}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-18, 10}, {18, -10}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Text(origin = {-63, 7}, extent = {{-15, 5}, {15, -5}}, textString = "Temp IN"), Text(origin = {8, 37}, extent = {{-30, 7}, {30, -7}}, textString = "Radiator"), Rectangle(origin = {-19.3958, 14.7922}, extent = {{-63.2606, 33.1616}, {63.2606, -33.1616}}), Text(origin = {23, 9}, extent = {{-15, 5}, {15, -5}}, textString = "Temp OUT")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end Radiator;

    model DHRS
    
      parameter SMD_MSR_Modelica.Units.Mass m_DHRS;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity Cp_fluid;
      parameter SMD_MSR_Modelica.Units.MassFlowRate m_dotNom;
      parameter SMD_MSR_Modelica.Units.TimeConstant DHRS_timeConstant;
      parameter SMD_MSR_Modelica.Units.HeatFlowRate DHRS_MaxPowerRm;
      parameter SMD_MSR_Modelica.Units.HeatFlowRate DHRS_PowerBleed;
      parameter SMD_MSR_Modelica.Units.InitiationTime DHRS_time;
      parameter SMD_MSR_Modelica.Units.Conductivity KF;
      parameter SMD_MSR_Modelica.Units.Area Ac;
      parameter SMD_MSR_Modelica.Units.Area Ar;
      parameter SMD_MSR_Modelica.Units.Length L;
      parameter SMD_MSR_Modelica.Units.Emissivity e;
      parameter SMD_MSR_Modelica.Units.Temperature Tinf;
      parameter SMD_MSR_Modelica.Units.Temperature T_0;
      constant SMD_MSR_Modelica.Units.StefanBoltzmannConstant SigSBK = 5.670374419E-14;
      
      
      SMD_MSR_Modelica.Units.MassFlowRate m_dot;
      SMD_MSR_Modelica.Units.HeatFlowRate DHRS_PowerRm;
      
      input SMD_MSR_Modelica.PortsConnectors.volSpecificHeat_in pDecay annotation(
        Placement(visible = true, transformation(origin = {-60, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In fuelFlowFrac annotation(
        Placement(visible = true, transformation(origin = {-60, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.Temp_In DHRS_tempIN annotation(
        Placement(visible = true, transformation(origin = {-60, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out DHRS_TempOUT annotation(
        Placement(visible = true, transformation(origin = {40, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {40, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
    initial equation
      DHRS_TempOUT.T = T_0;
      
    equation
      m_dot = m_dotNom*fuelFlowFrac.FF;
      
      DHRS_PowerRm = (DHRS_MaxPowerRm - DHRS_PowerBleed)/(1 + exp(log(1/1E-3 - 1)*(1 - (time - DHRS_time)/DHRS_timeConstant))) + DHRS_PowerBleed;
      
      if fuelFlowFrac.FF > 0 then
      
        m_DHRS*Cp_fluid*der(DHRS_TempOUT.T) = m_dot*Cp_fluid*(DHRS_tempIN.T - DHRS_TempOUT.T) + pDecay.Q*m_DHRS- DHRS_PowerRm + e*SigSBK*Ar*((Tinf + 273.15)^4 - (DHRS_TempOUT.T + 273.15)^4);
        
        
      else
      
        m_DHRS*Cp_fluid*der(DHRS_TempOUT.T) = ((KF*Ac)/L)*(DHRS_tempIN.T - DHRS_TempOUT.T) + pDecay.Q*m_DHRS - DHRS_PowerRm + e*SigSBK*Ar*((Tinf + 273.15)^4 - (DHRS_TempOUT.T + 273.15)^4);
    
        
      end if;
      
      
    protected
      annotation(
        Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(origin = {-10, 20}, lineThickness = 2, extent = {{-70, 60}, {70, -60}}), Text(origin = {41, 9}, extent = {{-15, 11}, {15, -11}}, textString = "Tout", fontSize = 20), Text(origin = {-59, 7}, extent = {{-15, 11}, {15, -11}}, textString = "Tin", fontSize = 20), Text(origin = {21, 62}, extent = {{-31, 12}, {31, -12}}, textString = "DHRS"), Text(origin = {-59, -31}, extent = {{-15, 11}, {15, -11}}, textString = "FF", fontSize = 20), Text(origin = {-59, 48}, extent = {{-13, 8}, {13, -8}}, textString = "DH", fontSize = 20)}),
        Icon(coordinateSystem(extent = {{-80, 80}, {60, -40}}), graphics = {Rectangle(origin = {-10, 20}, lineThickness = 2, extent = {{-70, 60}, {70, -60}}), Text(origin = {21, 62}, extent = {{-31, 12}, {31, -12}}, textString = "DHRS"), Text(origin = {-59, 9}, extent = {{-15, 11}, {15, -11}}, textString = "Tin", fontSize = 20), Text(origin = {-59, 50}, extent = {{-13, 8}, {13, -8}}, textString = "DH", fontSize = 20), Text(origin = {41, 9}, extent = {{-15, 11}, {15, -11}}, textString = "Tout", fontSize = 20), Text(origin = {-59, -31}, extent = {{-15, 11}, {15, -11}}, textString = "FF", fontSize = 20)}));
    end DHRS;

    model Pipe
      parameter SMD_MSR_Modelica.Units.Mass mPi;
      parameter Real mFracNode = 0.1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity CpFluid;
      parameter SMD_MSR_Modelica.Units.MassFlowRate mDotNom;
      parameter SMD_MSR_Modelica.Units.Conductivity KF;
      parameter SMD_MSR_Modelica.Units.Area Ac;
      parameter SMD_MSR_Modelica.Units.Area Ar;
      parameter SMD_MSR_Modelica.Units.Length L;
      parameter SMD_MSR_Modelica.Units.Emissivity e;
      parameter SMD_MSR_Modelica.Units.Temperature Tinf;
      parameter SMD_MSR_Modelica.Units.Temperature T_0;
      constant SMD_MSR_Modelica.Units.StefanBoltzmannConstant SigSBK = 5.670374419E-14;
      
      SMD_MSR_Modelica.Units.Mass repMPn;
      SMD_MSR_Modelica.Units.ResidentTime TauPi;
      SMD_MSR_Modelica.Units.ResidentTime repTauPi;
      SMD_MSR_Modelica.Units.ResidentTime varTauToPn;
      SMD_MSR_Modelica.Units.ResidentTime varTauPnTo;
      SMD_MSR_Modelica.Units.Temperature PnTempIn;
      SMD_MSR_Modelica.Units.MassFlowRate m_dotPi;
      SMD_MSR_Modelica.Units.Temperature PnTemp;
      
      input SMD_MSR_Modelica.PortsConnectors.Temp_In PiTemp_IN annotation(
        Placement(visible = true, transformation(origin = {-40, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out PiTemp_Out annotation(
        Placement(visible = true, transformation(origin = {28, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {28, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.volSpecificHeat_in PiDecay_Heat annotation(
        Placement(visible = true, transformation(origin = {-40, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In fuelFlowFrac annotation(
        Placement(visible = true, transformation(origin = {-40, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
    initial equation
      PiTemp_Out.T = T_0;
    
    equation
      TauPi = mPi/mDotNom;
      repMPn = mPi*(1 - mFracNode);
      repTauPi = (mPi - repMPn)/mDotNom;
      m_dotPi = mDotNom*fuelFlowFrac.FF;
      
      if fuelFlowFrac.FF > 0 then
        varTauToPn = (repTauPi/2)/fuelFlowFrac.FF;
        varTauPnTo = (repTauPi/2)/fuelFlowFrac.FF;
        PnTempIn = delay(PiTemp_IN.T, varTauToPn, 50000);
        repMPn*CpFluid*der(PnTemp) = m_dotPi*CpFluid*(PnTempIn - PnTemp) + PiDecay_Heat.Q*mPi + e*SigSBK*Ar*((Tinf + 273.15)^4 - (PnTemp + 273.15)^4);
        PiTemp_Out.T = delay(PnTemp, varTauPnTo, 50000);
      else
        varTauToPn = 0;
        varTauPnTo = 0;
        PnTempIn = PiTemp_IN.T;
        mPi*CpFluid*der(PnTemp) = ((KF*Ac)/L)*(PnTempIn - PnTemp) + PiDecay_Heat.Q*mPi + e*SigSBK*Ar*((Tinf + 273.15)^4 - (PnTemp + 273.15)^4);
        PiTemp_Out.T = PnTemp;
      end if;
      annotation(
        Diagram(graphics = {Rectangle(origin = {-10, 24}, lineThickness = 1, extent = {{-54, 50}, {54, -50}}), Text(origin = {-40, 15}, extent = {{-14, 9}, {14, -9}}, textString = "Tin", fontSize = 20), Text(origin = {28, 13}, extent = {{-14, 9}, {14, -9}}, textString = "Tout", fontSize = 20), Text(origin = {19, 61}, extent = {{-21, 9}, {21, -9}}, textString = "Pipe"), Text(origin = {-40, 49}, extent = {{-14, 9}, {14, -9}}, textString = "DH", fontSize = 20), Text(origin = {-40, -15}, extent = {{-14, 9}, {14, -9}}, textString = "FF", fontSize = 20)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Text(origin = {-40, -15}, extent = {{-14, 9}, {14, -9}}, textString = "FF", fontSize = 20), Text(origin = {-40, 49}, extent = {{-14, 9}, {14, -9}}, textString = "DH", fontSize = 20), Rectangle(origin = {-10, 24}, lineThickness = 1, extent = {{-54, 50}, {54, -50}}), Text(origin = {19, 61}, extent = {{-21, 9}, {21, -9}}, textString = "PipeNode"), Text(origin = {28, 13}, extent = {{-14, 9}, {14, -9}}, textString = "Tout", fontSize = 20), Text(origin = {-40, 15}, extent = {{-14, 9}, {14, -9}}, textString = "Tin", fontSize = 20)}));
    end Pipe;

    model PrimaryPump
      parameter SMD_MSR_Modelica.Units.FlowFraction freeConvectionFF;
      parameter SMD_MSR_Modelica.Units.PumpConstant primaryPumpK;
      parameter SMD_MSR_Modelica.Units.InitiationTime tripPrimaryPump;
      output SMD_MSR_Modelica.PortsConnectors.FlowFraction_Out primaryFlowFrac annotation(
        Placement(visible = true, transformation(origin = {-1.77636e-15, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {1.77636e-15, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    initial equation
      primaryFlowFrac.FF = 1;
    equation
      primaryFlowFrac.FF = (1 - freeConvectionFF)*exp(-(1/primaryPumpK)*delay(time, tripPrimaryPump)) + freeConvectionFF;
      annotation(
        Diagram(graphics = {Rectangle(lineColor = {245, 121, 0}, lineThickness = 1, extent = {{-60, 60}, {60, -60}}), Text(origin = {2, 6}, extent = {{-42, 16}, {42, -16}}, textString = "Primary Pump")}),
        Icon(graphics = {Rectangle(origin = {0, -20}, lineColor = {245, 121, 0}, lineThickness = 1, extent = {{-60, 60}, {60, -60}}), Text(origin = {4, -20}, extent = {{-42, 16}, {42, -16}}, textString = "Primary Pump")}));
    end PrimaryPump;

    model SecondaryPump
      parameter SMD_MSR_Modelica.Units.FlowFraction freeConvectionFF;
      parameter SMD_MSR_Modelica.Units.PumpConstant secondaryPumpK;
      parameter SMD_MSR_Modelica.Units.InitiationTime tripSecondaryPump;
      output SMD_MSR_Modelica.PortsConnectors.FlowFraction_Out secondaryFlowFrac annotation(
        Placement(visible = true, transformation(origin = {-1.77636e-15, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {1.77636e-15, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    initial equation
//secondaryFlowFrac.FF = 1;
    equation
      secondaryFlowFrac.FF = (1 - freeConvectionFF)*exp(-(1/secondaryPumpK)*delay(time, tripSecondaryPump)) + freeConvectionFF;
      annotation(
        Diagram(graphics = {Rectangle(lineColor = {245, 121, 0}, lineThickness = 1, extent = {{-60, 60}, {60, -60}}), Text(origin = {2, 3}, extent = {{-52, 27}, {52, -27}}, textString = "Secondary pump")}),
        Icon(graphics = {Rectangle(origin = {0, -20}, lineColor = {245, 121, 0}, lineThickness = 1, extent = {{-60, 60}, {60, -60}}), Text(origin = {2, -19}, extent = {{-52, 27}, {52, -27}}, textString = "Secondary pump")}));
    end SecondaryPump;

    model UHX
      parameter SMD_MSR_Modelica.Units.Mass m_PN1UHX;
      parameter SMD_MSR_Modelica.Units.ReactorPower P;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_pFluidUHX;
      parameter SMD_MSR_Modelica.Units.MassFlowRate m_dot_pFluidUHXnom;
      parameter SMD_MSR_Modelica.Units.Temperature Tp_0;
      SMD_MSR_Modelica.Units.MassFlowRate m_dot_pFluidUHX;
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In secondaryFlowFrac annotation(
        Placement(visible = true, transformation(origin = {-60, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out UHXtemp_Out annotation(
        Placement(visible = true, transformation(origin = {20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.Temp_In UHXtemp_In annotation(
        Placement(visible = true, transformation(origin = {-60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.RealIn PowerDemand annotation(
        Placement(visible = true, transformation(origin = {-60, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      UHXtemp_Out.T = Tp_0;
    equation
      m_dot_pFluidUHX = m_dot_pFluidUHXnom*secondaryFlowFrac.FF;
      m_PN1UHX*cP_pFluidUHX*der(UHXtemp_Out.T) = m_dot_pFluidUHX*cP_pFluidUHX*(UHXtemp_In.T - UHXtemp_Out.T) - P*PowerDemand.R;
      annotation(
        Diagram(graphics = {Rectangle(origin = {-20, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {-60, -14}, extent = {{-14, 12}, {14, -12}}, textString = "Tin", fontSize = 20), Text(origin = {8, 44}, extent = {{-34, 18}, {34, -18}}, textString = "UHX"), Text(origin = {-60, 26}, extent = {{-14, 12}, {14, -12}}, textString = "FF", fontSize = 20), Text(origin = {-60, -52}, extent = {{-14, 12}, {14, -12}}, textString = "Qrm", fontSize = 20), Text(origin = {22, -14}, extent = {{-16, 12}, {16, -12}}, textString = "Tout", fontSize = 20)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Text(origin = {20, -14}, extent = {{-16, 12}, {16, -12}}, textString = "TpOUT", fontSize = 20), Text(origin = {-60, -52}, extent = {{-14, 12}, {14, -12}}, textString = "Qrm", fontSize = 20), Rectangle(origin = {-20, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {-60, 26}, extent = {{-14, 12}, {14, -12}}, textString = "FFp", fontSize = 20), Text(origin = {-60, -14}, extent = {{-14, 12}, {14, -12}}, textString = "TpIN", fontSize = 20), Text(origin = {8, 44}, extent = {{-34, 18}, {34, -18}}, textString = "UHX")}));
    end UHX;

    model MixingPot
      parameter Integer numInput = 4;
      parameter SMD_MSR_Modelica.Units.Mass m;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity Cp;
      parameter SMD_MSR_Modelica.Units.MassFlowRate mDotNom;
      parameter SMD_MSR_Modelica.Units.ResidentTime tauPotToDHRS;
      parameter SMD_MSR_Modelica.Units.Temperature T_0;
      parameter SMD_MSR_Modelica.Units.FlowFraction regionFlowFrac[numInput];
      parameter SMD_MSR_Modelica.Units.Conductivity KF;
      parameter SMD_MSR_Modelica.Units.Area Ac;
      parameter SMD_MSR_Modelica.Units.Area Ar;
      parameter SMD_MSR_Modelica.Units.Length L;
      parameter SMD_MSR_Modelica.Units.Emissivity e;
      parameter SMD_MSR_Modelica.Units.Temperature Tinf;
      constant SMD_MSR_Modelica.Units.StefanBoltzmannConstant SigSBK = 5.670374419E-14;
      SMD_MSR_Modelica.Units.MassFlowRate mDotIn[numInput];
      SMD_MSR_Modelica.Units.HeatFlowRate Qin[numInput];
      SMD_MSR_Modelica.Units.MassFlowRate mDotOut;
      SMD_MSR_Modelica.Units.ResidentTime varTauPotToDHRS;
      SMD_MSR_Modelica.Units.Temperature T;
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out temp_Out annotation(
        Placement(visible = true, transformation(origin = {60, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {60, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In flowFraction[numInput] annotation(
        Placement(visible = true, transformation(origin = {-58, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-58, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.Temp_In temp_In[numInput] annotation(
        Placement(visible = true, transformation(origin = {-60, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-60, -2}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.volSpecificHeat_in decayHeat annotation(
        Placement(visible = true, transformation(origin = {-58, 58}, extent = {{-22, -22}, {22, 22}}, rotation = 0), iconTransformation(origin = {-58, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      T = T_0;
    equation
      for i in 1:numInput loop
        mDotIn[i] = mDotNom*flowFraction[i].FF*regionFlowFrac[i];
        Qin[i] = mDotIn[i]*Cp*temp_In[i].T;
      end for;
      mDotOut = sum(mDotIn);
      if flowFraction[numInput].FF > 0 then
        varTauPotToDHRS = tauPotToDHRS/sum(flowFraction.FF);
        m*Cp*der(T) = sum(Qin) - mDotOut*Cp*T + decayHeat.Q*m;
        temp_Out.T = delay(T, varTauPotToDHRS, 5000);
      else
        varTauPotToDHRS = 0;
        m*Cp*der(T) = ((KF*Ac)/L)*((sum(temp_In[:].T)/numInput) - T) + e*SigSBK*Ar*((Tinf + 273.15)^4 - (T + 273.15)^4) + decayHeat.Q*m;
        temp_Out.T = T;
      end if;
      annotation(
        Diagram(graphics = {Rectangle(origin = {0, -6}, lineThickness = 2, extent = {{-80, 86}, {80, -86}}), Text(origin = {50, 69}, extent = {{-26, 15}, {26, -15}}, textString = "Mixing Pot"), Text(origin = {-58, -20}, extent = {{-16, 6}, {16, -6}}, textString = "T_In"), Text(origin = {-58, -80}, extent = {{-16, 6}, {16, -6}}, textString = "FF"), Text(origin = {-58, 38}, extent = {{-16, 6}, {16, -6}}, textString = "DH"), Text(origin = {62, -18}, extent = {{-16, 6}, {16, -6}}, textString = "T_Out")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Text(origin = {-56, -80}, extent = {{-16, 6}, {16, -6}}, textString = "FlowFrac"), Text(origin = {-58, -20}, extent = {{-16, 6}, {16, -6}}, textString = "TempIn"), Text(origin = {50, 21}, extent = {{-26, 15}, {26, -15}}, textString = "Mixing Pot"), Text(origin = {62, -18}, extent = {{-16, 6}, {16, -6}}, textString = "TempOut"), Text(origin = {58, -80}, extent = {{-16, 6}, {16, -6}}, textString = "DecayHeat"), Rectangle(origin = {1, -13}, extent = {{-79, 79}, {79, -79}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end MixingPot;

    model FlowDistributor
      parameter Integer numOutput;
      parameter SMD_MSR_Modelica.Units.FlowFraction freeConvectionFF;
      parameter SMD_MSR_Modelica.Units.PumpConstant coastDownK;
      parameter SMD_MSR_Modelica.Units.InitiationTime regionTripTime[numOutput];
      output SMD_MSR_Modelica.PortsConnectors.FlowFraction_Out flowFracOut[numOutput] annotation(
        Placement(visible = true, transformation(origin = {39, -39}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {39, -39}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In flowFracIn annotation(
        Placement(visible = true, transformation(origin = {-41, -39}, extent = {{-15, -15}, {15, 15}}, rotation = 0), iconTransformation(origin = {-41, -39}, extent = {{-15, -15}, {15, 15}}, rotation = 0)));
    equation
      for i in 1:numOutput loop
        flowFracOut[i].FF = (flowFracIn.FF*(1 - freeConvectionFF)*exp(-coastDownK*delay(time, regionTripTime[i]))) + freeConvectionFF;
      end for;
      annotation(
        Diagram(graphics = {Rectangle(origin = {0, -30}, lineColor = {230, 97, 0}, lineThickness = 1.5, extent = {{-60, 30}, {60, -30}}), Text(origin = {15, -9}, extent = {{-43, 23}, {43, -23}}, textString = "FlowDistributor")}),
        Icon(graphics = {Text(origin = {15, -9}, extent = {{-43, 23}, {43, -23}}, textString = "FlowDistributor"), Rectangle(origin = {0, -30}, lineColor = {230, 97, 0}, lineThickness = 1.5, extent = {{-60, 30}, {60, -30}})}));
    end FlowDistributor;

    model TransportDelay
      parameter SMD_MSR_Modelica.Units.ResidentTime tau;
      parameter SMD_MSR_Modelica.Units.FlowFraction freeFF;
      SMD_MSR_Modelica.Units.ResidentTime varTau;
      SMD_MSR_Modelica.Units.ResidentTime maxTau;
      SMD_MSR_Modelica.Units.Temperature tempIn;
      SMD_MSR_Modelica.Units.Temperature tempOut;
      SMD_MSR_Modelica.Units.FlowFraction flowFrac;
    equation
//Maximum transport delay @ lowest system FF -> tau/freeFF
      maxTau = tau/freeFF;
//Variable transport delay -> tau/FF
      varTau = tau/flowFrac;
      tempOut = delay(tempIn, varTau, maxTau);
    end TransportDelay;

    model Pump
      parameter Integer numRampUp=1;
      //parameter Integer numRampDown=1;
      parameter SMD_MSR_Modelica.Units.PumpConstant rampUpK[numRampUp];
      parameter SMD_MSR_Modelica.Units.FlowFraction rampUpTo[numRampUp];
      parameter SMD_MSR_Modelica.Units.InitiationTime rampUpTime[numRampUp];
      
      parameter SMD_MSR_Modelica.Units.PumpConstant tripK;
      parameter SMD_MSR_Modelica.Units.InitiationTime tripTime;
      //parameter SMD_MSR_Modelica.Units.PumpConstant rampDownK[numRampDown];
      //parameter SMD_MSR_Modelica.Units.FlowFraction rampDownTo[numRampDown];
      //parameter SMD_MSR_Modelica.Units.InitiationTime rampDownTime[numRampDown];
      parameter SMD_MSR_Modelica.Units.FlowFraction freeConvFF;
      
        
      constant Real epsilon = 1E-3;
      SMD_MSR_Modelica.Units.FlowFraction rampUp[numRampUp];
      SMD_MSR_Modelica.Units.FlowFraction rampDown;
      SMD_MSR_Modelica.Units.FlowFraction ramp;
      
      output SMD_MSR_Modelica.PortsConnectors.FlowFraction_Out primaryFlowFrac annotation(
        Placement(visible = true, transformation(origin = {-1.77636e-15, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {1.77636e-15, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    
    initial equation

    primaryFlowFrac.FF = 0;
    
    equation
    
      for i in 1:numRampUp loop
        if i == 1 then
          rampUp[i] = (rampUpTo[i]-freeConvFF)/(1 + exp(log(1/epsilon - 1)*(1 - (delay(time,rampUpTime[i]))/rampUpK[i])));
        else
          rampUp[i] = (rampUpTo[i] - rampUpTo[i - 1])/(1 + exp(log(1/epsilon - 1)*(1 - (delay(time,rampUpTime[i]))/rampUpK[i])));
        end if;
      end for;
      
      if time >= tripTime then
        rampDown = (-sum(rampUp))*exp(-(1/tripK)*delay(time, tripTime)) + sum(rampUp);
      else 
        rampDown = 0;
      end if;    
    
      /*for i in 1:numRampDown loop
          rampDown[i] = -(rampDownTo[i])*exp(-rampDownK[i]*delay(time, rampDownTime[i])) + rampDownTo[i]; //elseif i == numRampDown then
      end for;*/
      
    
      ramp = sum(rampUp) - rampDown + freeConvFF;  
        
        
      if ramp > 0.001 then
        primaryFlowFrac.FF = ramp;
      else
        primaryFlowFrac.FF = 0;
      end if;
      annotation(
        Diagram(graphics = {Rectangle(lineColor = {245, 121, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {-2, 6}, extent = {{-42, 16}, {42, -16}}, textString = "Pump")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Rectangle(origin = {0, -20}, lineColor = {245, 121, 0}, lineThickness = 1, extent = {{-60, 60}, {60, -60}}), Text(origin = {4, -20}, extent = {{-42, 16}, {42, -16}}, textString = "Primary Pump")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end Pump;
  end HeatTransport;

  package Nuclear
    model mPKE
      //Parameter decleration
      parameter Integer numGroups;
      parameter SMD_MSR_Modelica.Units.NeutronGenerationTime LAMBDA;
      parameter SMD_MSR_Modelica.Units.PrecursorDecayConstant lambda[numGroups];
      parameter SMD_MSR_Modelica.Units.DelayedNeutronFrac beta[numGroups];
      parameter SMD_MSR_Modelica.Units.ResidentTime tauCoreNom;
      parameter SMD_MSR_Modelica.Units.ResidentTime tauLoopNom;
      //Variable decleration
      SMD_MSR_Modelica.Units.ResidentTime varTauCore;
      SMD_MSR_Modelica.Units.ResidentTime varTauLoop;
      SMD_MSR_Modelica.Units.PrecursorConc CG[numGroups];
      SMD_MSR_Modelica.Units.PrecursorConc CG_0sta[numGroups];
      SMD_MSR_Modelica.Units.PrecursorReturnRate CGReturn[numGroups];
      SMD_MSR_Modelica.Units.PrecursorDecayRate CGDecay[numGroups];
      SMD_MSR_Modelica.Units.DelayedNeutronFrac rho_0;
      SMD_MSR_Modelica.Units.Reactivity reactivity;
      input SMD_MSR_Modelica.PortsConnectors.ReactivityIn feedback annotation(
        Placement(visible = true, transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.NomNeutronPopOut n_population annotation(
        Placement(visible = true, transformation(origin = {0, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In fuelFlowFrac annotation(
        Placement(visible = true, transformation(origin = {34, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {34, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.RealIn externalReact annotation(
        Placement(visible = true, transformation(origin = {-32, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-30, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      n_population.n = 1;
      for i in 1:numGroups loop
        CG_0sta[i] = beta[i]/(1 + (1/(lambda[i]*tauCoreNom))*(1 - exp(-lambda[i]*tauLoopNom)));
        CG[i] = CG_0sta[i];
      end for;
    equation
      for i in 1:numGroups loop
        CG_0sta[i] = beta[i]/(1 + (1/(lambda[i]*tauCoreNom))*(1 - exp(-lambda[i]*tauLoopNom)));
        CGDecay[i] = lambda[i]*CG[i];
        CGReturn[i] = delay(CG[i], varTauLoop, 5000)*exp(-lambda[i]*varTauLoop)/varTauCore;
        der(CG[i]) = (beta[i]*n_population.n)/LAMBDA - CGDecay[i] - (CG[i]/varTauCore) + CGReturn[i];
      end for;
      rho_0 = sum(beta) - sum(CG_0sta);
      varTauCore = tauCoreNom/fuelFlowFrac.FF;
      varTauLoop = tauLoopNom/fuelFlowFrac.FF;
      reactivity = feedback.rho + rho_0 + externalReact.R*1E-5;
      der(n_population.n) = ((reactivity - sum(beta))/LAMBDA)*n_population.n + sum(CGDecay);
      annotation(
        Diagram(graphics = {Rectangle(origin = {2, 0}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, lineThickness = 0.5, extent = {{-50, 50}, {50, -50}}), Text(origin = {3, 1}, extent = {{37, -27}, {-37, 27}}, textString = "mPKE")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Rectangle(origin = {2, 0}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, lineThickness = 2, extent = {{-50, 50}, {50, -50}}), Text(origin = {3, 1}, extent = {{37, -27}, {-37, 27}}, textString = "mPKE")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end mPKE;

    model DecayHeat
      parameter SMD_MSR_Modelica.Units.DecayHeatYield DHYG1;
      parameter SMD_MSR_Modelica.Units.DecayHeatYield DHYG2;
      parameter SMD_MSR_Modelica.Units.DecayHeatYield DHYG3;
      parameter SMD_MSR_Modelica.Units.DecayHeatPrecursorDecayConstant DHlamG1;
      parameter SMD_MSR_Modelica.Units.DecayHeatPrecursorDecayConstant DHlamG2;
      parameter SMD_MSR_Modelica.Units.DecayHeatPrecursorDecayConstant DHlamG3;
      SMD_MSR_Modelica.Units.HeatTransferFraction DHG1;
      SMD_MSR_Modelica.Units.HeatTransferFraction DHG2;
      SMD_MSR_Modelica.Units.HeatTransferFraction DHG3;
      input SMD_MSR_Modelica.PortsConnectors.NomNeutronPopIn nPop annotation(
        Placement(visible = true, transformation(origin = {1, 33}, extent = {{-19, -19}, {19, 19}}, rotation = 0), iconTransformation(origin = {9, 41}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.DecayHeat_Out decayHeat_Out annotation(
        Placement(visible = true, transformation(origin = {-4.44089e-16, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {10, -20}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
    initial equation
     
      DHG1 = nPop.n*DHYG1/DHlamG1;
      DHG2 = nPop.n*DHYG2/DHlamG2;
      DHG3 = nPop.n*DHYG3/DHlamG3;
    
    equation
        der(DHG1) = nPop.n*DHYG1 - DHlamG1*DHG1;
        der(DHG2) = nPop.n*DHYG2 - DHlamG2*DHG2;
        der(DHG3) = nPop.n*DHYG3 - DHlamG3*DHG3;
      decayHeat_Out.DH = (DHG1 + DHG2 + DHG3);
      annotation(
        Diagram(graphics = {Rectangle(origin = {2, 0}, lineColor = {0, 225, 255}, lineThickness = 0.75, extent = {{-50, 50}, {50, -50}}), Text(origin = {3, 0}, extent = {{-45, 46}, {45, -46}}, textString = "Decay Heat"), Text(origin = {1, 17}, extent = {{9, -9}, {-9, 9}}, textString = "n", fontSize = 20), Text(origin = {0, -45}, extent = {{12, -11}, {-12, 11}}, textString = "DH", fontSize = 20)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Text(origin = {9, 9}, extent = {{-43, 25}, {43, -25}}, textString = "DH"), Rectangle(origin = {10, 10}, lineColor = {0, 225, 255}, lineThickness = 2, extent = {{-50, 50}, {50, -50}}), Text(origin = {1, 36}, extent = {{1, 2}, {-1, -2}}, textString = "text")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end DecayHeat;

    model ReactivityFeedback
      parameter SMD_MSR_Modelica.Units.Temperature FuelTempSetPointNode1;
      parameter SMD_MSR_Modelica.Units.Temperature FuelTempSetPointNode2;
      parameter SMD_MSR_Modelica.Units.Temperature GrapTempSetPoint;
      parameter SMD_MSR_Modelica.Units.TemperatureReactivityCoef a_F;
      parameter SMD_MSR_Modelica.Units.TemperatureReactivityCoef a_G;
      parameter Real IF1;
      parameter Real IF2;
      parameter Real IG;
      SMD_MSR_Modelica.Units.Reactivity FuelTempFeedbackNode1;
      SMD_MSR_Modelica.Units.Reactivity FuelTempFeedbackNode2;
      SMD_MSR_Modelica.Units.Reactivity GrapTempFeedback;
      SMD_MSR_Modelica.Units.Reactivity TotalTempFeedback;
      input SMD_MSR_Modelica.PortsConnectors.Temp_In fuelNode1 annotation(
        Placement(visible = true, transformation(origin = {-40, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.Temp_In fuelNode2 annotation(
        Placement(visible = true, transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-10, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.Temp_In grapNode annotation(
        Placement(visible = true, transformation(origin = {-40, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {20, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.ReactivityOut feedback annotation(
        Placement(visible = true, transformation(origin = {32, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-10, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      FuelTempFeedbackNode1 = (fuelNode1.T - FuelTempSetPointNode1)*IF1*a_F;
      FuelTempFeedbackNode2 = (fuelNode2.T - FuelTempSetPointNode2)*IF2*a_F;
      GrapTempFeedback = (grapNode.T - GrapTempSetPoint)*a_G*IG;
      TotalTempFeedback = FuelTempFeedbackNode1 + FuelTempFeedbackNode2 + GrapTempFeedback;
      feedback.rho = TotalTempFeedback;
      annotation(
        Diagram(graphics = {Rectangle(origin = {0, 2}, extent = {{-50, 50}, {50, -50}}), Text(origin = {13, 40}, extent = {{-31, 24}, {31, -24}}, textString = "Feedback"), Text(origin = {-32, 22}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN1"), Text(origin = {-32, -12}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN2"), Text(origin = {-32, -40}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_GrapN1"), Text(origin = {24, 4}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_Feedback")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Text(origin = {-40, 50}, extent = {{-16, 6}, {16, -6}}, textString = "F1"), Text(origin = {-10, 8}, extent = {{-32, 24}, {32, -24}}, textString = "TRFB"), Text(origin = {-10, 50}, extent = {{-16, 6}, {16, -6}}, textString = "F2"), Rectangle(origin = {-10, 10}, lineColor = {11, 200, 36}, lineThickness = 2, extent = {{-50, 50}, {50, -50}}), Text(origin = {20, 50}, extent = {{-16, 6}, {16, -6}}, textString = "G")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end ReactivityFeedback;

    model FuelChannel
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
      hA = hAnom*(0.8215*fuelFlowFraction.FF^6 - 4.108*fuelFlowFraction.FF^5 + 7.848*fuelFlowFraction.FF^4 - 7.165*fuelFlowFraction.FF^3 + 3.004*fuelFlowFraction.FF^2 + 0.5903*fuelFlowFraction.FF + 0.008537);
      
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
    end FuelChannel;

    model Poisons
      parameter SMD_MSR_Modelica.Units.IsotopicDecayConstant Te135_lam = 1;
      parameter SMD_MSR_Modelica.Units.IsotopicDecayConstant I135_lam = 1;
      parameter SMD_MSR_Modelica.Units.IsotopicDecayConstant Xe135_lam = 1;
      parameter SMD_MSR_Modelica.Units.IsotopicFissYield Pm149_lam = 1;
      parameter SMD_MSR_Modelica.Units.IsotopicDecayConstant Sm149_lam = 1;
      parameter SMD_MSR_Modelica.Units.IsotopicFissYield Te135_gamma = 1;
      parameter SMD_MSR_Modelica.Units.IsotopicFissYield I135_gamma = 1;
      parameter SMD_MSR_Modelica.Units.IsotopicFissYield Xe135_gamma = 1;
      parameter SMD_MSR_Modelica.Units.IsotopicFissYield Pm149_gamma = 1;
      SMD_MSR_Modelica.Units.IsotopicConc Te135_Conc;
      SMD_MSR_Modelica.Units.IsotopicConc I135_Conc;
      SMD_MSR_Modelica.Units.IsotopicConc Xe135_Conc;
      SMD_MSR_Modleica.Units.IsotopicConc Pm149_Conc;
      SMD_MSR_Modelica.Units.IsotopicConc Sm149_Conc;
    initial equation

    equation
      der(Te135_Conc) = Te135_gamma*Sig_fission*phi0;
    end Poisons;

    model PowerBlock
      parameter SMD_MSR_Modelica.Units.ReactorPower P;
      parameter SMD_MSR_Modelica.Units.Mass TotalFuelMass;
      SMD_MSR_Modelica.Units.NominalPower nomFissionPower;
      SMD_MSR_Modelica.Units.NominalPower nomReactorPower;
      SMD_MSR_Modelica.Units.ReactorPower decayPower;
      SMD_MSR_Modelica.Units.ReactorPower reactorPower;
      input SMD_MSR_Modelica.PortsConnectors.NomNeutronPopIn nPop annotation(
        Placement(visible = true, transformation(origin = {-52, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.DecayHeat_In decayNP annotation(
        Placement(visible = true, transformation(origin = {32, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.PowerOut fissionPower annotation(
        Placement(visible = true, transformation(origin = {-16, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-20, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.volSpecificHeat_Out decayPowerM annotation(
        Placement(visible = true, transformation(origin = {38, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {40, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      nomFissionPower = nPop.n*(1 - 0.068);
      fissionPower.P = nomFissionPower*P;
      decayPower = decayNP.DH*P;
      decayPowerM.Q = decayPower/TotalFuelMass;
      nomReactorPower = nomFissionPower + decayNP.DH;
      reactorPower = fissionPower.P + decayPower;
      annotation(
        Diagram(graphics = {Rectangle(origin = {-10, 10}, lineColor = {224, 27, 36}, extent = {{-70, 70}, {70, -70}}), Text(origin = {-11, 10}, extent = {{61, 46}, {-61, -46}}, textString = "Power Block")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Rectangle(origin = {10, -30}, lineColor = {224, 27, 36}, lineThickness = 2, extent = {{-50, 50}, {50, -50}}), Text(origin = {11, -31}, extent = {{33, 21}, {-33, -21}}, textString = "PB"), Text(origin = {-19, -72}, extent = {{-9, 6}, {9, -6}}, textString = "FP"), Text(origin = {41, -72}, extent = {{-9, 6}, {9, -6}}, textString = "DP")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end PowerBlock;

    model SumReactivity
      parameter Integer numInput;
      input SMD_MSR_Modelica.PortsConnectors.ReactivityIn reactivityIn[numInput] annotation(
        Placement(visible = true, transformation(origin = {20, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {8, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.ReactivityOut reactivityOut annotation(
        Placement(visible = true, transformation(origin = {-20, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {10, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      reactivityOut.rho = sum(reactivityIn.rho);
      annotation(
        Diagram(graphics = {Rectangle(origin = {0, 20}, lineColor = {51, 209, 122}, extent = {{-40, 20}, {40, -20}}), Text(origin = {11, 33}, extent = {{-29, 7}, {29, -7}}, textString = "Sum React"), Text(origin = {21, 9}, extent = {{-9, 5}, {9, -5}}, textString = "In"), Text(origin = {-19, 9}, extent = {{-9, 7}, {9, -7}}, textString = "Out")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Text(origin = {10, -10}, extent = {{-40, 16}, {40, -16}}, textString = "SumR"), Rectangle(origin = {10, -10}, lineColor = {11, 200, 36}, lineThickness = 2, extent = {{-50, 50}, {50, -50}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end SumReactivity;

    model PKE
      //Parameter decleration
      parameter SMD_MSR_Modelica.Units.NominalNeutronPopulation n_0;
      parameter Integer numGroups= 6;
      parameter SMD_MSR_Modelica.Units.NeutronGenerationTime LAMBDA;
      parameter SMD_MSR_Modelica.Units.PrecursorDecayConstant lambda[numGroups];
      parameter SMD_MSR_Modelica.Units.DelayedNeutronFrac beta[numGroups];
      parameter SMD_MSR_Modelica.Units.ResidentTime nominalTauCore;
      parameter SMD_MSR_Modelica.Units.ResidentTime nominalTauLoop;
      parameter Boolean addRho0;
      //Variable decleration
      
      SMD_MSR_Modelica.Units.DelayedNeutronFrac sumBeta;
      SMD_MSR_Modelica.Units.DelayedNeutronFrac sumCG_0dyn;
      SMD_MSR_Modelica.Units.DelayedNeutronFrac rho_0dyn;
      SMD_MSR_Modelica.Units.DelayedNeutronFrac rho_0sta;
      //SMD_MSR_Modelica.Units.DelayedNeutronFrac beta_eff;
      SMD_MSR_Modelica.Units.ResidentTime varTauCore;
      SMD_MSR_Modelica.Units.ResidentTime varTauLoop;
      SMD_MSR_Modelica.Units.PrecursorConc CG[numGroups];
      SMD_MSR_Modelica.Units.PrecursorConc CG_0dyn[numGroups];
      SMD_MSR_Modelica.Units.PrecursorConc CG_0sta[numGroups];
      SMD_MSR_Modelica.Units.PrecursorReturnRate CGReturn[numGroups];
      SMD_MSR_Modelica.Units.PrecursorDecayRate CGDecay[numGroups];
      SMD_MSR_Modelica.Units.Reactivity reactivity;
      SMD_MSR_Modelica.Units.Reactivity externalReactivityIn;
      SMD_MSR_Modelica.Units.NomalizedNeutronEmissionRate nomS;
      input SMD_MSR_Modelica.PortsConnectors.ReactivityIn feedback annotation(
        Placement(visible = true, transformation(origin = {18, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In fuelFlowFrac annotation(
        Placement(visible = true, transformation(origin = {42, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {20, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.NeutronEmsRateIn S annotation(
        Placement(visible = true, transformation(origin = {-6, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-20, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.NomNeutronPopOut n_population annotation(
        Placement(visible = true, transformation(origin = {0, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-10, -40}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.RealIn ReactivityIn annotation(
        Placement(visible = true, transformation(origin = {-36, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
    initial equation
      n_population.n = n_0;
      
      if (n_0 >= 0 and fuelFlowFrac.FF > 0) then
          for i in 1:numGroups loop
            CG[i] = (beta[i]/LAMBDA)*n_population.n*(1.0/(lambda[i] - (exp(-lambda[i]*nominalTauLoop) - 1.0)/nominalTauCore));
          end for;
          varTauCore = nominalTauCore;
          varTauLoop = nominalTauLoop;
               
      end if;
    
    equation
    
      der(n_population.n) = ((reactivity - sum(beta))/LAMBDA)*n_population.n + sum(CGDecay) + nomS;
      //beta_eff = sum(CG);
      
      if addRho0== true then
        reactivity = feedback.rho + externalReactivityIn + rho_0sta;
      else 
        reactivity = feedback.rho + externalReactivityIn;
      end if;  
      
      sumBeta = sum(beta);
      sumCG_0dyn = sum(CG_0dyn);
      
      externalReactivityIn = ReactivityIn.R*1E-5;
      nomS = S.nDot/1.58E20;
      for i in 1:numGroups loop
        CG_0sta[i] = beta[i]/(1 + (1/(lambda[i]*nominalTauCore))*(1 - exp(-lambda[i]*nominalTauLoop)));
      end for;
      rho_0sta = sumBeta - sum(CG_0sta);
      
//noEvent is used to prevent compilor from symbolically eveluating both branches of the if statement
//if noEvent removed the model tends to break when flow fraction cross zero
      if noEvent(fuelFlowFrac.FF > 0) then
        for i in 1:numGroups loop
          CGDecay[i] = lambda[i]*CG[i];
          CGReturn[i] = (delay(CG[i], varTauLoop, 5000000)*exp(-lambda[i]*varTauLoop))/varTauCore;
          der(CG[i]) = (beta[i]*n_population.n)/LAMBDA - CGDecay[i] - (CG[i]/varTauCore) + CGReturn[i];
          CG_0dyn[i] = beta[i]/(1 + (1/(lambda[i]*varTauCore))*(1 - exp(-lambda[i]*varTauLoop)));
        end for;
        varTauCore = nominalTauCore/fuelFlowFrac.FF;
        varTauLoop = nominalTauLoop/fuelFlowFrac.FF;
        
        rho_0dyn = sumBeta - sumCG_0dyn;
      
      else
        for i in 1:numGroups loop
          CGDecay[i] = lambda[i]*CG[i];
          CGReturn[i] = 0;
          der(CG[i]) = (beta[i]*n_population.n)/LAMBDA - CGDecay[i];
          CG_0dyn[i] = 0;
        end for;
        rho_0dyn = 0;
        varTauCore = 0;
        varTauLoop = 0;
        
      end if;
      annotation(
        Diagram(graphics = {Rectangle(origin = {2, 0}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, lineThickness = 2, extent = {{-50, 50}, {50, -50}}), Text(origin = {3, 1}, extent = {{37, -27}, {-37, 27}}, textString = "PKE")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Rectangle(origin = {-10, -10}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, lineThickness = 1.5, extent = {{-50, 50}, {50, -50}}), Text(origin = {-11, -11}, extent = {{37, -27}, {-37, 27}}, textString = "PKE")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end PKE;
  end Nuclear;

  package Units
    type Area = Real(unit = "m2", min = 0);
    //Area
    type Length = Real(unit = "m", min = 0);
    //Length
    type Volume = Real(unit = "m3", min = 0);
    type Frequency = Real(unit = "rad/s");
    type Emissivity = Real(unit = "1", min = 0);
    type StefanBoltzmannConstant = Real(unit = "MW/(m2.K4)");
    type AvagadroNumber = Real(unit = "1/mol");
    type MolarMass = Real(unit = "kg/mol", min = 0);
    type MoleFraction = Real(unit = "1", min = 0);
    type Enrichment = Real(unit = "1", min = 0);
    type EnergyPerFission = Real(unit = "J", min = 0);
    type CrossSection = Real(unit = "m2", min = 0);
    type Density = Real(unit = "kg/m3", min = 0);
    type Speed = Real(unit = "m/s", min = 0);
    type MultipicationFactor = Real(unit = "1", min = 0);
    //Frequency [rad/s]
    // Nuclear - PKE
    type NeutronDensity = Real(unit = "1/m3", min = 0);
    //Neutron Density [#n/m^3]
    type NominalNeutronPopulation = Real(unit = "1", min = 0);
    type NeutronPopulation = Real(unit = "1", min = 0);
    //Nominl Neutron Population
    type NeutronEmissionRate = Real(unit = "1/s", min = 0);
    type NomalizedNeutronEmissionRate = Real(unit = "1/s", min = 0);
    type AtomicDensity = Real(unit = "1/m3", min = 0);
    type NeutronGenerationTime = Real(unit = "s", min = 0);
    //Neutron Generation Time (Lambda) [s]
    type Reactivity = Real(unit = "1");
    //Absolute Reactivity [unitless]
    type PrecursorDecayConstant = Real(unit = "1/s", min = 0);
    //Delayed Precursor Decay Constant [s^-1]
    type PrecursorConc = Real(unit = "1", min = 0);
    //Delayed Neutron Precursor
    type PrecursorDecayRate = Real(unit = "1/s", min = 0);
    type PrecursorReturnRate = Real(unit = "1/s", min = 0);
    type DelayedNeutronFrac = Real(unit = "1", min = 0);
    //Delayed Neutron Fraction [unitless]
    type VolumeImportance = Real(unit = "1", min = 0);
    //Fraction of Fission Power Deposited in Core Nodes [m^3]
    // Nuclear - Temperature feedback
    type TemperatureReactivityCoef = Real(unit = "1/C");
    //Absolute Temperature Reactivity Coef [C^-1]
    // Nuclear - Poisons
    type IsotopicConc = Real(unit = "1/m^3", min = 0);
    //Isotopic Concentration [#atoms/m^3]
    type IsotopicDecayConstant = Real(unit = "1/s", min = 0);
    //Isotopic Decay Constant [s^-1]
    type IsotopicFissYield = Real(unit = "1", min = 0);
    //Isotopic Yeild per fission
    // Nuclear - Model Specs
    type NeutronFlux = Real(unit = "n/(cm2.s)");
    // Flow
    type MassFlowRate = Real(unit = "kg/s", min = 0);
    //Mass Flow Rate [kg/m^3]
    type ResidentTime = Real(unit = "s", min = 0, max = 1E6);
    //Fuel Salt Resident Time (tau) [s]
    type FlowFraction = Real(unit = "1", min = 0);
    //Flow Fraction [unitless]
    type SpecificHeat = Real(unit = "MJ/(kg.s)");
    type HeatFlowRate = Real(unit = "MJ/s");
    //Heat Flow Rate [MJ/kg.s]
    // Heat Transfer
    type NominalPower = Real(unit = "1", min = 0);
    //Nominal reactor power as a fraction
    type RadPower = Real(unit = "MW");
    type ReactorPower = Real(unit = "MW", min = 0);
    //Nominal Reactor Thermal Power
    type Conductivity = Real(unit = "MW/(m.C)", min = 0);
    //Heat conduction [MW/(m.C) == MJ/(s.m.C)]
    type Convection = Real(unit = "MJ/(s.C)", min = 0);
    //Rate of heat transfer by convection [J/(s.C) == W/C]
    type HeatTransferFraction = Real(unit = "1", min = 0);
    //Fraction of Heat Trasfer from Node to an adjacent Node [unitless]
    // Decay Heat
    type DecayHeatPrecursorDecayConstant = Real(unit = "1/s", min = 0);
    //Decay heat precoursor decay constant [s^-1]
    type DecayHeatYield = Real(unit = "1/s", min = 0);
    type DecayHeatFraction = Real(unit = "1", min = 0);
    // Physical Properties
    type Mass = Real(unit = "kg", min = 0);
    type SpecificHeatCapacity = Real(unit = "MJ/(kg.C)", min = 0);
    type Temperature = Real(unit = "C");
    // Component Related
    type TimeConstant = Real(unit = "1/s", min = 0);
    //DHRS time constant [per seconds]
    type InitiationTime = Real(unit = "s", min = 0);
    //Start component [seconds]
    type Demand = Real(unit = "1", min = 0);
    //Component demand
    type PumpConstant = Real(unit = "1/s", min = 0);
  end Units;

  package PortsConnectors
    connector Temp_In
      SMD_MSR_Modelica.Units.Temperature T;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {204, 0, 0}, fillColor = {204, 0, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {204, 0, 0}, fillColor = {204, 0, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end Temp_In;

    connector Temp_Out
      SMD_MSR_Modelica.Units.Temperature T;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {204, 0, 0}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {204, 0, 0}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end Temp_Out;

    connector DecayHeat_In
      SMD_MSR_Modelica.Units.NominalPower DH;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end DecayHeat_In;

    connector DecayHeat_Out
      SMD_MSR_Modelica.Units.NominalPower DH;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end DecayHeat_Out;

    connector NomNeutronPopIn
      SMD_MSR_Modelica.Units.NominalNeutronPopulation n;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end NomNeutronPopIn;

    connector ReactivityOut
      SMD_MSR_Modelica.Units.Reactivity rho;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end ReactivityOut;

    connector ReactivityIn
      SMD_MSR_Modelica.Units.Reactivity rho;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end ReactivityIn;

    connector FlowFraction_In
      SMD_MSR_Modelica.Units.FlowFraction FF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {255, 120, 0}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {255, 120, 0}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end FlowFraction_In;

    connector FlowFraction_Out
      SMD_MSR_Modelica.Units.FlowFraction FF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {245, 121, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {245, 121, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end FlowFraction_Out;

    connector NeutronEmsRateIn
      SMD_MSR_Modelica.Units.NeutronEmissionRate nDot;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {191, 64, 188}, fillColor = {191, 64, 188}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {191, 64, 188}, fillColor = {191, 64, 188}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end NeutronEmsRateIn;

    connector NomNeutronPopOut
      SMD_MSR_Modelica.Units.NominalNeutronPopulation n;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Ellipse(origin = {1, -1}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end NomNeutronPopOut;

    connector NeutronEmsRateOut
      SMD_MSR_Modelica.Units.NeutronEmissionRate nDot;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {191, 64, 188}, fillColor = {191, 64, 188}, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {191, 64, 188}, fillColor = {191, 64, 188}, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end NeutronEmsRateOut;

    connector RealIn
      Real R;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}}), Ellipse(origin = {-8, 36}, extent = {{0, -2}, {0, 2}}), Text(extent = {{2, -10}, {2, -10}}, textString = "text")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end RealIn;

    connector RealOut
      Real R;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end RealOut;

    connector PowerIn
      SMD_MSR_Modelica.Units.ReactorPower P;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {129, 61, 156}, fillColor = {129, 61, 156}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Ellipse(origin = {-1, -1}, lineColor = {129, 61, 156}, fillColor = {129, 61, 156}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end PowerIn;

    connector PowerOut
      SMD_MSR_Modelica.Units.ReactorPower P;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Ellipse(origin = {-1, -1}, lineColor = {129, 61, 156}, fillColor = {129, 61, 156}, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end PowerOut;

    connector volSpecificHeat_in
      SMD_MSR_Modelica.Units.SpecificHeat Q;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {220, 138, 221}, fillColor = {220, 138, 221}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {220, 138, 221}, fillColor = {220, 138, 221}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end volSpecificHeat_in;

    connector volSpecificHeat_Out
      SMD_MSR_Modelica.Units.SpecificHeat Q;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {220, 138, 221}, fillColor = {0, 225, 255}, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {220, 138, 221}, fillColor = {0, 225, 255}, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
    end volSpecificHeat_Out;
  end PortsConnectors;

  package QAtoolBox
    package TestIO
      model NeutronPopulation
        parameter SMD_MSR_Modelica.Units.NominalNeutronPopulation relNpop;
        output SMD_MSR_Modelica.PortsConnectors.NominalNeutronPopulation n_population annotation(
          Placement(visible = true, transformation(origin = {-1, 29}, extent = {{-29, -29}, {29, 29}}, rotation = 0), iconTransformation(origin = {4, 38}, extent = {{-28, -28}, {28, 28}}, rotation = 0)));
      initial equation
        n_population.n = 0.999;
      equation
        n_population.n = 0.999 + (if time < 8000 then 0 else 0.1*sin(0.01*time));
        annotation(
          Diagram(graphics = {Rectangle(origin = {0, -2}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-60, 60}, {60, -60}}), Text(origin = {3, -28}, extent = {{-51, 16}, {51, -16}}, textString = "n")}),
          Icon(graphics = {Text(origin = {7, -35}, extent = {{-53, 21}, {53, -21}}, textString = "n"), Rectangle(origin = {4, -6}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, lineThickness = 2, extent = {{-76, 72}, {76, -72}})}));
      end NeutronPopulation;

      model SteadyStateTempFeedback
        parameter SMD_MSR_Modelica.Units.Reactivity Feedback;
        output SMD_MSR_Modelica.PortsConnectors.ReactivityOut reactivity annotation(
          Placement(visible = true, transformation(origin = {0, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        when time > 0 then
//reactivity.rho = -0.00264445;
          reactivity.rho = -1533E-5;
//step1 k=0.104
        elsewhen time >= 3600 then
          reactivity.rho = -1364E-5;
//step2 k=0.2211
        elsewhen time >= 7200 then
          reactivity.rho = -1241E-5;
//step3 k=0.3083
        elsewhen time >= 10800 then
          reactivity.rho = -1076.7E-5;
//step4 k=0.4247
        elsewhen time >= 14400 then
          reactivity.rho = -935.7E-5;
//step5 k=0.5239
        elsewhen time >= 18000 then
          reactivity.rho = -786.2E-5;
//step6 k=0.6304
        elsewhen time >= 21600 then
          reactivity.rho = -638.15E-5;
//step7 k=0.7353
        elsewhen time >= 25200 then
          reactivity.rho = -482.1596E-5;
//step8 k=0.845
        elsewhen time >= 28800 then
          reactivity.rho = -312.3573E-5;
//step9 k=0.9661
        elsewhen time >= 32400 then
          reactivity.rho = -264.445E-5;
//step10 k=0.99
        elsewhen time >= 36000 then
          reactivity.rho = -199.4E-5;
//step11 k=0.9059 after pump ramp to 25% --> k=0.952
        elsewhen time >= 37500 then
          reactivity.rho = -158.8E-5;
//step12 k=0.9808
        elsewhen time >= 39500 then
          reactivity.rho = -140.415E-5;
//step13 k=0.9881
        elsewhen time >= 43000 then
          reactivity.rho = -131.678E-5;
//step14 k=0.9957
        elsewhen time >= 48500 then
          reactivity.rho = -89.1478E-5;
//step15 k=0.9613 after pump ramp to 50% --> k=0.9920
        elsewhen time >= 51000 then
          reactivity.rho = -77.9088E-5;
//step16 k=0.996
        elsewhen time >= 54500 then
          reactivity.rho = -61.9895E-5;
//step17 k=0.9575 after pump ramp to 100% --> 0.9688
        elsewhen time >= 55500 then
          reactivity.rho = -17.9305E-5;
//step18 k=0.996
//elsewhen time >= 70000 then
//reactivity.rho = 0;
//elsewhen time >= 80000 then
//reactivity.rho = (time-80000)*1E-6; //step18 k=0.996
//elsewhen time >= 82470 then
//reactivity.rho = 247E-5;
//elsewhen time >= 80000  then
//reactivity.rho = 247E-5;
        end when;
        annotation(
          Diagram(graphics = {Rectangle(lineColor = {78, 154, 6}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {2, -36}, extent = {{-40, 20}, {40, -20}}, textString = "rho_{tempFB}")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
          Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Text(origin = {2, -36}, extent = {{-40, 20}, {40, -20}}, textString = "rho_{tempFB}"), Rectangle(lineColor = {78, 154, 6}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}));
      end SteadyStateTempFeedback;

      model ConstantTempIn
        parameter SMD_MSR_Modelica.Units.Temperature TestTempIn;
        output SMD_MSR_Modelica.PortsConnectors.Temp_Out temp_Out annotation(
          Placement(visible = true, transformation(origin = {-1.77636e-15, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-1.77636e-15, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      equation
        temp_Out.T = TestTempIn;
        annotation(
          Diagram(graphics = {Rectangle(lineColor = {204, 0, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {1, -31}, extent = {{-25, 13}, {25, -13}}, textString = "T")}),
          Icon(graphics = {Rectangle(lineColor = {204, 0, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {1, -31}, extent = {{-25, 13}, {25, -13}}, textString = "T")}));
      end ConstantTempIn;

      model ConstantTempOut
        parameter SMD_MSR_Modelica.Units.Temperature ConstantTempOut;
        output SMD_MSR_Modelica.PortsConnectors.Temp_Out temp_Out annotation(
          Placement(visible = true, transformation(origin = {-1.77636e-15, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-1.77636e-15, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      equation
        temp_Out.T = ConstantTempOut;
        annotation(
          Diagram(graphics = {Rectangle(lineColor = {204, 0, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}),
          Icon(graphics = {Rectangle(lineColor = {204, 0, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}));
      end ConstantTempOut;

      model ConstantFlowFrac
        parameter SMD_MSR_Modelica.Units.FlowFraction ConstantFlowFrac;
        output SMD_MSR_Modelica.PortsConnectors.FlowFraction_Out flowFraction_Out annotation(
          Placement(visible = true, transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 38}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
      initial equation
        flowFraction_Out.FF = ConstantFlowFrac;
      equation
        flowFraction_Out.FF = ConstantFlowFrac;
        annotation(
          Diagram(graphics = {Rectangle(lineColor = {245, 121, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {3, -32}, extent = {{-45, 10}, {45, -10}}, textString = "FlowFrac")}),
          Icon(graphics = {Rectangle(lineColor = {245, 121, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {3, -32}, extent = {{-45, 10}, {45, -10}}, textString = "FlowFrac")}));
      end ConstantFlowFrac;

      model ConstantDecayPower
        parameter SMD_MSR_Modelica.Units.SpecificHeat DecayPower = 0;
    output SMD_MSR_Modelica.PortsConnectors.volSpecificHeat_Out decayHeat_Out annotation(
          Placement(visible = true, transformation(origin = {-4, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-4, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        decayHeat_Out.Q = DecayPower;
        annotation(
          Diagram(graphics = {Rectangle(lineColor = {220, 138, 221}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {-2, -33}, extent = {{-42, 13}, {42, -13}}, textString = "P_{Decay}")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
          Icon(graphics = {Rectangle(lineColor = {220, 138, 221}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {2, -35}, extent = {{-42, 13}, {42, -13}}, textString = "P_{Decay}")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
      end ConstantDecayPower;

      model ConstantHeatProduction
        parameter SMD_MSR_Modelica.Units.MassFlowRate mDot;
        parameter SMD_MSR_Modelica.Units.Mass nodeMass;
        parameter SMD_MSR_Modelica.Units.HeatFlowRate qDotIn;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP;
        input SMD_MSR_Modelica.PortsConnectors.Temp_In temp_In annotation(
          Placement(visible = true, transformation(origin = {-40, -1.77636e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-40, -1.9984e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        output SMD_MSR_Modelica.PortsConnectors.Temp_Out temp_Out annotation(
          Placement(visible = true, transformation(origin = {40, -1.77636e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {40, -1.77636e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      equation
        nodeMass*cP*der(temp_Out.T) = mDot*cP*(temp_In.T - temp_Out.T) + qDotIn;
        annotation(
          Diagram(graphics = {Rectangle(lineColor = {204, 0, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {0, -42}, extent = {{-48, -30}, {48, 30}}, textString = "Q Production")}),
          Icon(graphics = {Text(origin = {0, -42}, extent = {{-48, -30}, {48, 30}}, textString = "Q Production"), Rectangle(lineColor = {204, 0, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}));
      end ConstantHeatProduction;

      model ConstantHeatRemoval
        parameter SMD_MSR_Modelica.Units.MassFlowRate mDot;
        parameter SMD_MSR_Modelica.Units.Mass nodeMass;
        parameter SMD_MSR_Modelica.Units.HeatFlowRate qDotOut;
        parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP;
        parameter SMD_MSR_Modelica.Units.Temperature tempInSetPoint;
        parameter SMD_MSR_Modelica.Units.Temperature tempoutSetPoint;
        input SMD_MSR_Modelica.PortsConnectors.Temp_In temp_In annotation(
          Placement(visible = true, transformation(origin = {-40, -1.77636e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-40, -1.9984e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        output SMD_MSR_Modelica.PortsConnectors.Temp_Out temp_Out annotation(
          Placement(visible = true, transformation(origin = {40, -1.77636e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {40, -1.77636e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      initial equation
        temp_In.T = tempInSetPoint;
        temp_Out.T = tempoutSetPoint;
      equation
        nodeMass*cP*der(temp_Out.T) = mDot*cP*(temp_In.T - temp_Out.T) - qDotOut;
        annotation(
          Diagram(graphics = {Rectangle(lineColor = {204, 0, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {0, -42}, extent = {{-48, -30}, {48, 30}}, textString = "Q Remove")}),
          Icon(graphics = {Text(origin = {0, -42}, extent = {{-48, -30}, {48, 30}}, textString = "Q Remove"), Rectangle(lineColor = {204, 0, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}));
      end ConstantHeatRemoval;

      model ConstantNeutronSource
        parameter SMD_MSR_Modelica.Units.NeutronEmissionRate S;
        output SMD_MSR_Modelica.PortsConnectors.NeutronEmsRateOut neutronEmsRateOut annotation(
          Placement(visible = true, transformation(origin = {-2, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-2, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        when time >= 10 then
          neutronEmsRateOut.nDot = S;
//elsewhen time >= 32400 then
        elsewhen time >= 32500 then
          neutronEmsRateOut.nDot = 0;
        elsewhen time >= 39600 then
//34200
          neutronEmsRateOut.nDot = S;
        elsewhen time >= 54100 then
//43200
          neutronEmsRateOut.nDot = 0;
        elsewhen time >= 61200 then
          neutronEmsRateOut.nDot = S;
        elsewhen time >= 68500 then
//68400
          neutronEmsRateOut.nDot = 0;
        elsewhen time >= 75600 then
          neutronEmsRateOut.nDot = S;
        elsewhen time >= 82900 then
//82800
          neutronEmsRateOut.nDot = 0;
        end when;
        annotation(
          Diagram(graphics = {Rectangle(lineColor = {191, 64, 188}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
          Icon(graphics = {Rectangle(lineColor = {191, 64, 188}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
      end ConstantNeutronSource;

      model ConstantNeutronSourceZero
        parameter SMD_MSR_Modelica.Units.NeutronEmissionRate S;
        output SMD_MSR_Modelica.PortsConnectors.NeutronEmsRateOut neutronEmsRateOut annotation(
          Placement(visible = true, transformation(origin = {-2, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-2, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        neutronEmsRateOut.nDot = S;
        annotation(
          Diagram(graphics = {Rectangle(lineColor = {191, 64, 188}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
          Icon(graphics = {Rectangle(lineColor = {191, 64, 188}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
      end ConstantNeutronSourceZero;
    end TestIO;

    package TestCases
      model TestPHX
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac primaryFF(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {-67, 71}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac secondaryFF(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {67, -67}, extent = {{-27, -27}, {27, 27}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0) annotation(
          Placement(visible = true, transformation(origin = {92, 70}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
        SMD_MSR_Modelica.HeatTransport.PrimaryHeatExchanger primaryHeatExchanger(m_PN1 = 87.0959, m_PN2 = 87.0959, m_PN3 = 87.0959, m_PN4 = 87.0959, m_TN1 = 101.1662, m_TN2 = 101.1662, m_SN1 = 24.9200, m_SN2 = 24.9200, m_SN3 = 24.9200, m_SN4 = 24.9200, cP_pFluid = 1.9665E-3, cP_Tube = 5.778E-04, cP_sFluid = 2.39E-3, m_dot_pFluidNom = 7.5708E-02*2.14647E+03, m_dot_sFluidNom = 5.36265E-02*1.922E3, hApnNom = 0.1620, hAsnNom = 0.0765, tauPHXtoPipe = 3.85, tauPHXtoUHX = 4.71) annotation(
          Placement(visible = true, transformation(origin = {-9, 7}, extent = {{-55, -55}, {55, 55}}, rotation = 0)));
      equation
        connect(primaryFF.flowFraction_Out, primaryHeatExchanger.primaryFF) annotation(
          Line(points = {{-66, 80}, {-48, 80}, {-48, 34}}, color = {245, 121, 0}));
        connect(secondaryFF.flowFraction_Out, primaryHeatExchanger.secondaryFF) annotation(
          Line(points = {{68, -56}, {62, -56}, {62, -20}}, color = {245, 121, 0}));
        connect(constantDecayPower.decayHeat_Out, primaryHeatExchanger.P_decay) annotation(
          Line(points = {{92, 80}, {6, 80}, {6, 34}}, color = {0, 225, 255}));
        connect(primaryHeatExchanger.T_in_pFluid, primaryHeatExchanger.T_out_pFluid) annotation(
          Line(points = {{-66, 24}, {80, 24}}));
        connect(primaryHeatExchanger.T_out_sFluid, primaryHeatExchanger.T_in_sFluid) annotation(
          Line(points = {{-64, -10}, {80, -10}}));
      end TestPHX;

      model TestMPKE
        SMD_MSR_Modelica.Nuclear.mPKE mPKE(Lam = 2.400E-04, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, lam = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, omega = 0, rho_0 = 0.002465140767843, sinInsertionTime = 0, sin_mag = 0, stepInsertionTime = 0, step_mag = 0, tauCoreNom = 8.46, tauLoopNom = 16.73) annotation(
          Placement(visible = true, transformation(origin = {3, -25}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.SteadyStateTempFeedback zeroFeedback(Feedback = 0) annotation(
          Placement(visible = true, transformation(origin = {-81, 21}, extent = {{-29, -29}, {29, 29}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac constantFlowFrac(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {74, 26}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
      equation
        connect(constantFlowFrac.flowFraction_Out, mPKE.fuelFlowFrac) annotation(
          Line(points = {{74, 34}, {19, 34}, {19, -6}}, color = {245, 121, 0}));
        connect(zeroFeedback.reactivity, mPKE.feedback) annotation(
          Line(points = {{-81, 29}, {4, 29}, {4, -6}}, color = {78, 154, 6}));
      end TestMPKE;

      model TestFeedback
        SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback annotation(
          Placement(visible = true, transformation(origin = {21, 19}, extent = {{-63, -63}, {63, 63}}, rotation = 0)));
      equation

      end TestFeedback;

      model TestFuelChannel
        SMD_MSR_Modelica.Nuclear.FuelChannel fuelChannel(P = 8, TotalFuelMass = 4.093499709644400e+03, hAnom = 0.0180, m_FN1 = 687.3959, m_FN2 = 687.3959, m_GN = 3.6342e+03, cP_fuel = 1.9665E-3, cP_graphite = 1.773E-3, mdot_fuelNom = 7.5708E-02*2.14647E+03, kFN1 = 0.465, kFN2 = 0.465, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, tauCoreToDHRS = 1.385) annotation(
          Placement(visible = true, transformation(origin = {8, 2}, extent = {{-48, -48}, {48, 48}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.NeutronPopulation neutronPopulation(relNpop = 0) annotation(
          Placement(visible = true, transformation(origin = {-121, 63}, extent = {{-27, -27}, {27, 27}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac constantFlowFrac(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {-120, -60}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0) annotation(
          Placement(visible = true, transformation(origin = {-120, -7.99361e-15}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
      equation
        connect(neutronPopulation.n_population, fuelChannel.nPop) annotation(
          Line(points = {{-120, 74}, {-30, 74}, {-30, 30}}, color = {20, 36, 248}));
        connect(constantDecayPower.decayHeat_Out, fuelChannel.P_decay) annotation(
          Line(points = {{-120, 12}, {-70, 12}, {-70, 0}, {-30, 0}}, color = {0, 225, 255}));
        connect(constantFlowFrac.flowFraction_Out, fuelChannel.fuelFlowFraction) annotation(
          Line(points = {{-120, -48}, {-12, -48}, {-12, -36}}, color = {245, 121, 0}));
        connect(fuelChannel.fuelNode1, fuelChannel.temp_In) annotation(
          Line(points = {{32, -32}, {-30, -32}}));
      protected
      end TestFuelChannel;

      model TestRadiator
        SMD_MSR_Modelica.HeatTransport.Radiator radiator(AirTempIn = 37.78, Cp_air = 1.0085E-3, Cp_salt = 2.39E-3, hA_radiatorNom = 0.0170, m_airN = 3.549878645002705, m_dot_air = 94.389*1.1237, m_dot_saltNom = 103.0701, m_saltN = 3.924401289158446e+02, tauUHXtoPHX = 8.24, tripRadiator = 2000000) annotation(
          Placement(visible = true, transformation(origin = {8, -10}, extent = {{-88, -88}, {88, 88}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac constantFlowFrac(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {-119, 11}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
      equation
        connect(constantFlowFrac.flowFraction_Out, radiator.secondaryFlowFrac) annotation(
          Line(points = {{-119, 24}, {-54, 24}, {-54, 26}}, color = {245, 121, 0}));
        connect(radiator.RadTempOut, radiator.RadTempIN) annotation(
          Line(points = {{32, 8}, {-54, 8}}));
      end TestRadiator;

      model TestDHRS
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs(Cp_fluid = 1.9665E-3, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 10, DHRS_timeConstant = 10, m_DHRS = 1.625049507600000e+02, m_Pi_C_PHX = 1, m_dotNom = 7.5708E-02*2.14647E+03, tauDHRStoPHX = 1.385) annotation(
          Placement(visible = true, transformation(origin = {36, -2}, extent = {{-70, -70}, {70, 70}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac constantFlowFrac(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {-80, -40}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0) annotation(
          Placement(visible = true, transformation(origin = {-87, 61}, extent = {{-25, -25}, {25, 25}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantTempIn constantTempIn(TestTempIn = 100) annotation(
          Placement(visible = true, transformation(origin = {-88, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
//connect(dhrs.DHRS_TempOUT, dhrs.DHRS_tempIN) annotation(
//    Line(points = {{67, 13}, {-6, 13}}));
        connect(constantFlowFrac.flowFraction_Out, dhrs.fuelFlowFrac) annotation(
          Line(points = {{-80, -28}, {-80, -16}, {-6, -16}}, color = {245, 121, 0}));
        connect(constantDecayPower.decayHeat_Out, dhrs.DHRS_DecayHeat) annotation(
          Line(points = {{-86, 72}, {-6, 72}, {-6, 40}}, color = {0, 225, 255}));
        connect(constantTempIn.temp_Out, dhrs.DHRS_tempIN) annotation(
          Line(points = {{-88, 12}, {-47, 12}, {-47, 14}, {-6, 14}}));
      end TestDHRS;

      model TestPipe
        SMD_MSR_Modelica.HeatTransport.Pipe pipe(Cp_fluidPi = 1000, m_Pi_PHX_C = 100, m_dotPiNom = 10, m_pi = 1, tauPiToCore = 10) annotation(
          Placement(visible = true, transformation(origin = {4, -2}, extent = {{-74, -74}, {74, 74}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0) annotation(
          Placement(visible = true, transformation(origin = {-132, 74}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac constantFlowFrac(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {-126, -66}, extent = {{-54, -54}, {54, 54}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantTempIn constantTempIn(TestTempIn = 100) annotation(
          Placement(visible = true, transformation(origin = {-131, 11}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
      equation
        connect(constantDecayPower.decayHeat_Out, pipe.PiDecay_Heat) annotation(
          Line(points = {{-132, 83}, {-32, 83}, {-32, 42}}, color = {0, 225, 255}));
        connect(constantFlowFrac.flowFraction_Out, pipe.fuelFlowFrac) annotation(
          Line(points = {{-126, -46}, {-32, -46}, {-32, -10}}, color = {245, 121, 0}));
        connect(constantTempIn.temp_Out, pipe.PiTemp_IN) annotation(
          Line(points = {{-130, 24}, {-32, 24}, {-32, 18}}));
      end TestPipe;

      model TestUHX
        SMD_MSR_Modelica.HeatTransport.UHX uhx(SetDemand = 1, UHXP = 1, cP_pFluidUHX = 1000, m_PN1UHX = 100, m_dot_pFluidUHXnom = 10, tauUHXtoPHX = 10, tripUHX = 0) annotation(
          Placement(visible = true, transformation(origin = {-1, 1}, extent = {{-55, -55}, {55, 55}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac constantFlowFrac(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {-133, 11}, extent = {{-39, -39}, {39, 39}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantTempIn constantTempIn(TestTempIn = 100) annotation(
          Placement(visible = true, transformation(origin = {-135, -73}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
      equation
        connect(constantFlowFrac.flowFraction_Out, uhx.secondaryFlowFrac) annotation(
          Line(points = {{-132, 26}, {-40, 26}, {-40, 24}}, color = {245, 121, 0}));
        connect(constantTempIn.temp_Out, uhx.UHXtemp_In) annotation(
          Line(points = {{-134, -64}, {-76, -64}, {-76, 0}, {-40, 0}}));
      end TestUHX;

      model TestPHX2
        SMD_MSR_Modelica.HeatTransport.PrimaryHeatExchanger PHX(m_PN1 = 87.0959, m_PN2 = 87.0959, m_PN3 = 87.0959, m_PN4 = 87.0959, m_TN1 = 101.1662, m_TN2 = 101.1662, m_SN1 = 24.9200, m_SN2 = 24.9200, m_SN3 = 24.9200, m_SN4 = 24.9200, cP_pFluid = 1.9665E-3, cP_Tube = 5.778E-04, cP_sFluid = 2.39E-3, m_dot_pFluidNom = 7.5708E-02*2.14647E+03, m_dot_sFluidNom = 5.36265E-02*1.922E3, hApnNom = 0.1620, hAsnNom = 0.0765, tauPHXtoPipe = 3.85, tauPHXtoUHX = 4.71) annotation(
          Placement(visible = true, transformation(origin = {-16, -6}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac primaryFlowFrac(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {-76, 54}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac secondaryFlowFrac(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {52, -74}, extent = {{-28, -28}, {28, 28}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0) annotation(
          Placement(visible = true, transformation(origin = {137, 29}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantHeatProduction constantHeatProduction(cP = 1.9665E-3, mDot = 7.5708E-02*2.14647E+03, nodeMass = 100, qDotIn = 10) annotation(
          Placement(visible = true, transformation(origin = {-140, 10}, extent = {{-34, -34}, {34, 34}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantHeatRemoval constantHeatRemoval(cP = 2.39E-3, mDot = 5.36265E-02*1.922E3, nodeMass = 100, qDotOut = 10) annotation(
          Placement(visible = true, transformation(origin = {131, -53}, extent = {{-33, -33}, {33, 33}}, rotation = 0)));
      equation
        connect(primaryFlowFrac.flowFraction_Out, PHX.primaryFF) annotation(
          Line(points = {{-76, 64}, {-52, 64}, {-52, 20}}, color = {245, 121, 0}));
        connect(secondaryFlowFrac.flowFraction_Out, PHX.secondaryFF) annotation(
          Line(points = {{52, -63}, {52, -32}}, color = {245, 121, 0}));
        connect(constantDecayPower.decayHeat_Out, PHX.P_decay) annotation(
          Line(points = {{137, 43}, {137, 43.5}, {-1, 43.5}, {-1, 20}}, color = {0, 225, 255}));
        connect(PHX.T_in_pFluid, constantHeatProduction.temp_Out) annotation(
          Line(points = {{-69, 10}, {-126, 10}}));
        connect(constantHeatProduction.temp_In, PHX.T_out_pFluid) annotation(
          Line(points = {{-154, 10}, {-168, 10}, {-168, 74}, {98, 74}, {98, 10}, {68, 10}}));
        connect(PHX.T_in_sFluid, constantHeatRemoval.temp_Out) annotation(
          Line(points = {{68, -23}, {170, -23}, {170, -52}, {144, -52}}));
        connect(constantHeatRemoval.temp_In, PHX.T_out_sFluid) annotation(
          Line(points = {{118, -52}, {-94, -52}, {-94, -22}, {-68, -22}}));
        annotation(
          Diagram(coordinateSystem(extent = {{-180, 80}, {180, -100}})));
      end TestPHX2;

      model TestFuelChannel2
        SMD_MSR_Modelica.Nuclear.FuelChannel fuelChannel(P = 8, TotalFuelMass = 4.093499709644400e+03, hAnom = 0.0180, m_FN1 = 687.3959, m_FN2 = 687.3959, m_GN = 3.6342e+03, cP_fuel = 1.9665E-3, cP_graphite = 1.773E-3, mdot_fuelNom = 7.5708E-02*2.14647E+03, kFN1 = 0.465, kFN2 = 0.465, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, tauCoreToDHRS = 1.385) annotation(
          Placement(visible = true, transformation(origin = {-72, -12}, extent = {{-52, -52}, {52, 52}}, rotation = 180)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantHeatRemoval constantHeatRemoval(cP = 1.9665E-3, mDot = 7.5708E-02*2.14647E+03, nodeMass = 100, qDotOut = 8) annotation(
          Placement(visible = true, transformation(origin = {-69, -87}, extent = {{-37, -37}, {37, 37}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac fuelFlowFrac(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {32, 32}, extent = {{-36, -36}, {36, 36}}, rotation = 0)));
        SMD_MSR_Modelica.Nuclear.mPKE mPKE(Lam = 2.400E-04, lam = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, tauCoreNom = 8.46, tauLoopNom = 16.73) annotation(
          Placement(visible = true, transformation(origin = {33, -55}, extent = {{-41, -41}, {41, 41}}, rotation = 0)));
        SMD_MSR_Modelica.Nuclear.ReactivityFeedback reactivityFeedback(FuelTempSetPointNode1 = 644.5000, FuelTempSetPointNode2 = 6.57E+02, GrapTempSetPoint = 6.756111E+02, a_F = -8.71E-05, a_G = -6.66E-05, omega = 1, sin_mag = 0, step_mag = 0, stepInsertionTime = 2000000, sinInsertionTime = 2000000) annotation(
          Placement(visible = true, transformation(origin = {-158, 4}, extent = {{-48, -48}, {48, 48}}, rotation = 180)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0) annotation(
          Placement(visible = true, transformation(origin = {33, -17}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
      equation
        connect(fuelChannel.fuelNode2, reactivityFeedback.fuelNode2) annotation(
          Line(points = {{-98, 4}, {-139, 4}}));
        connect(fuelChannel.fuelNode1, reactivityFeedback.fuelNode1) annotation(
          Line(points = {{-98, 26}, {-116, 26}, {-116, -12}, {-138, -12}}));
        connect(fuelChannel.grapNode, reactivityFeedback.grapNode) annotation(
          Line(points = {{-98, -20}, {-118, -20}, {-118, 18}, {-138, 18}}));
        connect(fuelChannel.temp_Out, constantHeatRemoval.temp_In) annotation(
          Line(points = {{-98, -40}, {-132, -40}, {-132, -86}, {-84, -86}}));
        connect(constantHeatRemoval.temp_Out, fuelChannel.temp_In) annotation(
          Line(points = {{-54, -86}, {-12, -86}, {-12, 26}, {-30, 26}, {-30, 24}}));
        connect(fuelChannel.nPop, mPKE.n_population) annotation(
          Line(points = {{-30, -42}, {0, -42}, {0, -84}, {34, -84}, {34, -70}}, color = {20, 36, 248}));
        connect(constantDecayPower.decayHeat_Out, fuelChannel.P_decay) annotation(
          Line(points = {{33, -9}, {-30, -9}, {-30, -8}}, color = {0, 225, 255}));
        connect(fuelFlowFrac.flowFraction_Out, fuelChannel.fuelFlowFraction) annotation(
          Line(points = {{32, 46}, {-52, 46}, {-52, 30}}, color = {245, 121, 0}));
        connect(mPKE.feedback, reactivityFeedback.feedback) annotation(
          Line(points = {{34, -38}, {0, -38}, {0, 58}, {-176, 58}, {-176, -4}, {-174, -4}}, color = {78, 154, 6}));
        connect(fuelFlowFrac.flowFraction_Out, mPKE.fuelFlowFrac) annotation(
          Line(points = {{32, 46}, {68, 46}, {68, -38}, {46, -38}}, color = {245, 121, 0}));
        annotation(
          Diagram(coordinateSystem(extent = {{-220, 80}, {80, -140}})));
      end TestFuelChannel2;

      model testSumReact
        SMD_MSR_Modelica.Nuclear.SumReactivity sumReactivity(numInput = 3) annotation(
          Placement(visible = true, transformation(origin = {-31, -17}, extent = {{-39, -39}, {39, 39}}, rotation = 0)));
        TestIO.SteadyStateTempFeedback R1(Feedback = 1) annotation(
          Placement(visible = true, transformation(origin = {30, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        TestIO.SteadyStateTempFeedback R2(Feedback = 2) annotation(
          Placement(visible = true, transformation(origin = {48, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        TestIO.SteadyStateTempFeedback R3(Feedback = 3) annotation(
          Placement(visible = true, transformation(origin = {44, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(R3.reactivity, sumReactivity.reactivityIn[3]) annotation(
          Line(points = {{44, -58}, {-24, -58}, {-24, -10}}, color = {78, 154, 6}));
        connect(R2.reactivity, sumReactivity.reactivityIn[2]) annotation(
          Line(points = {{48, -2}, {-24, -2}, {-24, -10}}, color = {78, 154, 6}));
        connect(R1.reactivity, sumReactivity.reactivityIn[1]) annotation(
          Line(points = {{30, 52}, {-24, 52}, {-24, -10}}, color = {78, 154, 6}));
      end testSumReact;

      model testFlowDistributor
        SMD_MSR_Modelica.HeatTransport.FlowDistributor flowDistributor(coastDownK = 0.02, freeConvectionFF = 0.005, numOutput = 4, regionFlowFrac = {0.05, 0.2, 0.1, 0.11}, regionTripTime = {200, 20000, 20000, 20000}) annotation(
          Placement(visible = true, transformation(origin = {-1, 27}, extent = {{-39, -39}, {39, 39}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac constantFlowFrac(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {-86, -50}, extent = {{-32, -32}, {32, 32}}, rotation = 0)));
      equation
        connect(constantFlowFrac.flowFraction_Out, flowDistributor.flowFracIn) annotation(
          Line(points = {{-86, -38}, {-16, -38}, {-16, 12}}, color = {245, 121, 0}));
      end testFlowDistributor;

      model TestPKE
        SMD_MSR_Modelica.QAtoolBox.TestIO.SteadyStateTempFeedback zeroFeedback(Feedback = 0) annotation(
          Placement(visible = true, transformation(origin = {14, 50}, extent = {{-48, -48}, {48, 48}}, rotation = 0)));
        SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, energyFiss = 200E6*1.6E-19, enrichment = 0.05, fissCrossSec = 5.851e-24, fuelMolFrac = 0.009, fuelSaltDensity = 2.14647E+03, fuelSlatMolMass = 0.0416, fuelVol = 0.6405, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, nThermSpeed = 2000, n_0 = 0, nomPower = 8, nominalTauCore = 8.46, nominalTauLoop = 16.73, numGroups = 6) annotation(
          Placement(visible = true, transformation(origin = {-58, -52}, extent = {{-64, -64}, {64, 64}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSource constantNeutronSource(S = 1E6) annotation(
          Placement(visible = true, transformation(origin = {-190, -38}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.SteadyStateTempFeedback extReact annotation(
          Placement(visible = true, transformation(origin = {-161, 53}, extent = {{-55, -55}, {55, 55}}, rotation = 0)));
        TestIO.ConstantFlowFrac constantFlowFrac(ConstantFlowFrac = 0) annotation(
          Placement(visible = true, transformation(origin = {90, 60}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
      equation
        connect(zeroFeedback.reactivity, pke.feedback) annotation(
          Line(points = {{14, 67}, {14, 39}, {-46, 39}, {-46, -26}}, color = {78, 154, 6}));
        connect(constantNeutronSource.neutronEmsRateOut, pke.S) annotation(
          Line(points = {{-191, -19}, {-191, -26}, {-64, -26}}, color = {191, 64, 188}));
        connect(extReact.reactivity, pke.externalReactivityIn) annotation(
          Line(points = {{-161, 73}, {-63, 73}, {-63, -26}, {-78, -26}}, color = {78, 154, 6}));
        connect(constantFlowFrac.flowFraction_Out, pke.fuelFlowFrac) annotation(
          Line(points = {{90, 69}, {90, -26}, {-34, -26}}, color = {245, 121, 0}));
      protected
        annotation(
          Diagram);
      end TestPKE;

      model TestPump
        SMD_MSR_Modelica.HeatTransport.Pump pump(
        freeConvFF = 0, numRampUp = 2, rampUpK = {50, 50}, 
        rampUpTime = {1, 3000}, 
        rampUpTo = {0.25, 0.5}, 
        tripK = 50, 
        tripTime = 7000) annotation(
          Placement(visible = true, transformation(origin = {-50, -10}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
      equation

      end TestPump;

      model TestPKEandPump
        SMD_MSR_Modelica.QAtoolBox.TestIO.SteadyStateTempFeedback zeroFeedback(Feedback = 0) annotation(
          Placement(visible = true, transformation(origin = {14, 50}, extent = {{-48, -48}, {48, 48}}, rotation = 0)));
        SMD_MSR_Modelica.Nuclear.PKE pke(LAMBDA = 2.400E-04, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, energyFiss = 200E6*1.6E-19, enrichment = 0.05, fissCrossSec = 5.851e-24, fuelMolFrac = 0.009, fuelSaltDensity = 2.14647E+03, fuelSlatMolMass = 0.0416, fuelVol = 0.6405, lambda = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, nThermSpeed = 2000, n_0 = 0, nomPower = 8, nominalTauCore = 8.46, nominalTauLoop = 16.73, numGroups = 6) annotation(
          Placement(visible = true, transformation(origin = {-58, -52}, extent = {{-64, -64}, {64, 64}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantNeutronSource constantNeutronSource(S = 1E6) annotation(
          Placement(visible = true, transformation(origin = {-190, -38}, extent = {{-56, -56}, {56, 56}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.SteadyStateTempFeedback extReact annotation(
          Placement(visible = true, transformation(origin = {-161, 53}, extent = {{-55, -55}, {55, 55}}, rotation = 0)));
        SMD_MSR_Modelica.HeatTransport.Pump pump(freeConvFF = 0, initialFF = 0, numRampUp = 2, rampTo = {0.5, 1}, rampUpK = {100, 100}, rampUpTime = {40000, 50000}) annotation(
          Placement(visible = true, transformation(origin = {134, -14}, extent = {{-28, -28}, {28, 28}}, rotation = 0)));
      equation
        connect(zeroFeedback.reactivity, pke.feedback) annotation(
          Line(points = {{14, 67}, {14, 39}, {-46, 39}, {-46, -26}}, color = {78, 154, 6}));
        connect(constantNeutronSource.neutronEmsRateOut, pke.S) annotation(
          Line(points = {{-191, -19}, {-191, -26}, {-64, -26}}, color = {191, 64, 188}));
        connect(extReact.reactivity, pke.externalReactivityIn) annotation(
          Line(points = {{-161, 73}, {-63, 73}, {-63, -26}, {-78, -26}}, color = {78, 154, 6}));
        connect(pump.primaryFlowFrac, pke.fuelFlowFrac) annotation(
          Line(points = {{134, -30}, {-34, -30}, {-34, -26}}, color = {245, 121, 0}));
      protected
        annotation(
          Diagram);
      end TestPKEandPump;
    end TestCases;
  end QAtoolBox;

  package InputSignals
    model Sigmoid
      parameter Real signalMax;
      parameter Real signalMin;
      parameter Real epsilon = 1E-3;
      parameter SMD_MSR_Modelica.Units.InitiationTime activationTime;
      parameter SMD_MSR_Modelica.Units.TimeConstant timeConstant;
      output SMD_MSR_Modelica.PortsConnectors.RealOut signal annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-24, -24}, {24, 24}}, rotation = 0), iconTransformation(origin = {1, -3}, extent = {{-27, -27}, {27, 27}}, rotation = 0)));
    equation
      signal.R = (signalMax - signalMin)/(1 + exp(log(1/epsilon - 1)*(1 - (time - activationTime)/timeConstant))) + signalMin;
      annotation(
        Diagram(graphics = {Rectangle(lineThickness = 1.5, extent = {{-60, 60}, {60, -60}})}));
    end Sigmoid;

    model Stepper
      parameter Integer numSteps;
      parameter SMD_MSR_Modelica.Units.InitiationTime stepTime[numSteps];
      parameter Real amplitude[numSteps];
      Real amp[numSteps];
      output SMD_MSR_Modelica.PortsConnectors.RealOut step annotation(
        Placement(visible = true, transformation(origin = {0, 22}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {4, 22}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    equation
      for i in 1:numSteps loop
        if time > stepTime[i] and i == 1 then
          amp[i] = amplitude[i];
        elseif time > stepTime[i] and i > 1 then
          amp[i] = amplitude[i] - amplitude[i - 1];
        else
          amp[i] = 0;
        end if;
      end for;
      step.R = sum(amp);
      annotation(
        Diagram(graphics = {Rectangle(origin = {3, -3}, lineThickness = 1, extent = {{-59, 57}, {59, -57}})}),
        Icon(graphics = {Rectangle(origin = {3, -3}, lineThickness = 1, extent = {{-59, 57}, {59, -57}})}));
    end Stepper;

    model Ramper
      parameter Integer numRamps;
      parameter SMD_MSR_Modelica.Units.InitiationTime rampTime[numRamps];
      parameter Real rate[numRamps];
      parameter Real offset;
      parameter Real maxValue;
      Real mag[numRamps];
      output SMD_MSR_Modelica.PortsConnectors.RealOut ramp annotation(
        Placement(visible = true, transformation(origin = {0, 22}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {4, 22}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    equation
      for i in 1:numRamps loop
        if time > rampTime[i] and i == 1 then
          mag[i] = rate[i]*(time - rampTime[i]);
        elseif time > rampTime[i] and i > 1 then
          mag[i] = (rate[i] - rate[i - 1])*(time - rampTime[i]);
        else
          mag[i] = 0;
        end if;
      end for;
      if maxValue > sum(mag) then
        ramp.R = sum(mag) + offset;
      else
        ramp.R = maxValue;
      end if;
      annotation(
        Diagram(graphics = {Rectangle(origin = {3, -3}, lineThickness = 1, extent = {{-59, 57}, {59, -57}})}),
        Icon(graphics = {Rectangle(origin = {3, -3}, lineThickness = 1, extent = {{-59, 57}, {59, -57}})}));
    end Ramper;

    model SumSignals
      parameter Integer numInput;
      input SMD_MSR_Modelica.PortsConnectors.RealIn realIn[numInput] annotation(
        Placement(visible = true, transformation(origin = {20, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {20, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.RealOut realOut annotation(
        Placement(visible = true, transformation(origin = {-18, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-18, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      realOut.R = sum(realIn.R);
      annotation(
        Diagram(graphics = {Rectangle(origin = {0, 20}, extent = {{-40, 20}, {40, -20}}), Text(origin = {11, 33}, extent = {{-29, 7}, {29, -7}}, textString = "Sum"), Text(origin = {21, 9}, extent = {{-9, 5}, {9, -5}}, textString = "In"), Text(origin = {-19, 9}, extent = {{-9, 7}, {9, -7}}, textString = "Out")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Text(origin = {-19, 9}, extent = {{-9, 7}, {9, -7}}, textString = "Out"), Text(origin = {21, 9}, extent = {{-9, 5}, {9, -5}}, textString = "In"), Text(origin = {11, 33}, extent = {{-29, 7}, {29, -7}}, textString = "Sum"), Rectangle(origin = {0, 20}, extent = {{-40, 20}, {40, -20}})}));
    end SumSignals;
    annotation(
      Diagram);
  end InputSignals;
  annotation(
    Documentation(info = "<html><head></head><body>SMD-MSR_ModelicaV1</body></html>"),
    uses(Modelica(version = "4.0.0")));
end SMD_MSR_Modelica;