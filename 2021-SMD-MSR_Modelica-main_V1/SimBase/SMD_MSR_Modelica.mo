package SMD_MSR_Modelica
  package HeatTransport
    model PrimaryHeatExchanger
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
      parameter SMD_MSR_Modelica.Units.Convection hApnNom;
      parameter SMD_MSR_Modelica.Units.Convection hAsnNom;
      parameter SMD_MSR_Modelica.Units.ResidentTime tauPHXtoPipe;
      parameter SMD_MSR_Modelica.Units.ResidentTime tauPHXtoUHX;
      parameter SMD_MSR_Modelica.Units.Temperature T_PN1_initial = 650.7500;
      parameter SMD_MSR_Modelica.Units.Temperature T_PN2_initial = 644.5000;
      parameter SMD_MSR_Modelica.Units.Temperature T_PN3_initial = 638.2500;
      parameter SMD_MSR_Modelica.Units.Temperature T_PN4_initial = 632;
      parameter SMD_MSR_Modelica.Units.Temperature T_TN1_initial = 625.1923;
      parameter SMD_MSR_Modelica.Units.Temperature T_TN2_initial = 611.3643;
      parameter SMD_MSR_Modelica.Units.Temperature T_SN1_initial = 554.4300;
      parameter SMD_MSR_Modelica.Units.Temperature T_SN2_initial = 562.7500;
      parameter SMD_MSR_Modelica.Units.Temperature T_SN3_initial = 571.0700;
      parameter SMD_MSR_Modelica.Units.Temperature T_SN4_initial = 579.3900;
      SMD_MSR_Modelica.Units.ResidentTime varTauPHXtoPipe;
      SMD_MSR_Modelica.Units.ResidentTime varTauPHXtoUHX;
      //SMD_MSR_Modelica.Units.ResidentTime maxTauPHXtoPipe;
      //SMD_MSR_Modelica.Units.ResidentTime maxTauPHXtoUHX;
      SMD_MSR_Modelica.Units.MassFlowRate m_dot_pFluid;
      SMD_MSR_Modelica.Units.MassFlowRate m_dot_sFluid;
      SMD_MSR_Modelica.Units.Convection hApn;
      SMD_MSR_Modelica.Units.Convection hAsn;
      SMD_MSR_Modelica.Units.Temperature T_PN1;
      SMD_MSR_Modelica.Units.Temperature T_PN2;
      SMD_MSR_Modelica.Units.Temperature T_PN3;
      SMD_MSR_Modelica.Units.Temperature T_PN4;
      SMD_MSR_Modelica.Units.Temperature T_TN1;
      SMD_MSR_Modelica.Units.Temperature T_TN2;
      SMD_MSR_Modelica.Units.Temperature T_SN1;
      SMD_MSR_Modelica.Units.Temperature T_SN2;
      SMD_MSR_Modelica.Units.Temperature T_SN3;
      SMD_MSR_Modelica.Units.Temperature T_SN4;
      input SMD_MSR_Modelica.PortsConnectors.Temp_In T_in_pFluid annotation(
        Placement(visible = true, transformation(origin = {-100, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-102, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.Temp_In T_in_sFluid annotation(
        Placement(visible = true, transformation(origin = {160, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {162, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out T_out_sFluid annotation(
        Placement(visible = true, transformation(origin = {-100, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out T_out_pFluid annotation(
        Placement(visible = true, transformation(origin = {160, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {162, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.DecayHeat_In P_decay annotation(
        Placement(visible = true, transformation(origin = {30, 50}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {28, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In primaryFF annotation(
        Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In secondaryFF annotation(
        Placement(visible = true, transformation(origin = {130, -50}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {130, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      T_PN1 = T_PN1_initial;
      T_PN2 = T_PN2_initial;
      T_PN3 = T_PN3_initial;
      T_PN4 = T_PN4_initial;
      T_TN1 = T_TN1_initial;
      T_TN2 = T_TN2_initial;
      T_SN1 = T_SN1_initial;
      T_SN2 = T_SN2_initial;
      T_SN3 = T_SN3_initial;
      T_SN4 = T_SN4_initial;
/*T_PN1 = 644.5625;
      T_PN2 = 644.5000;
      T_PN3 = 644.4375;
      T_PN4 = 644.3750;
      T_TN1 = 625.1923;
      T_TN2 = 611.3643;
      T_SN1 = 562.6668;
      T_SN2 = 562.7500;
      T_SN3 = 562.8332;
      T_SN4 = 562.9164;*/
    equation
      m_dot_pFluid = m_dot_pFluidNom * primaryFF.FF;
      m_dot_sFluid = m_dot_sFluidNom * secondaryFF.FF;
      varTauPHXtoPipe = tauPHXtoPipe / primaryFF.FF;
      varTauPHXtoUHX = tauPHXtoUHX / secondaryFF.FF;
//maxTauPHXtoPipe = tauPHXtoPipe/primaryFF.FF_min;
//maxTauPHXtoUHX = tauPHXtoUHX/secondaryFF.FF_min;
      hApn = hApnNom * (0.8215 * primaryFF.FF ^ 6 - 4.108 * primaryFF.FF ^ 5 + 7.848 * primaryFF.FF ^ 4 - 7.165 * primaryFF.FF ^ 3 + 3.004 * primaryFF.FF ^ 2 + 0.5903 * primaryFF.FF + 0.008537);
      hAsn = hAsnNom * (0.8215 * secondaryFF.FF ^ 6 - 4.108 * secondaryFF.FF ^ 5 + 7.848 * secondaryFF.FF ^ 4 - 7.165 * secondaryFF.FF ^ 3 + 3.004 * secondaryFF.FF ^ 2 + 0.5903 * secondaryFF.FF + 0.008537);
//hApn = hApnNom;
//hAsn = hAsnNom;
      m_PN1 * cP_pFluid * der(T_PN1) = m_dot_pFluid * cP_pFluid * (T_in_pFluid.T - T_PN1) - hApn * (T_PN1 - T_TN1) + P_decay.Q * m_PN1;
      m_PN2 * cP_pFluid * der(T_PN2) = m_dot_pFluid * cP_pFluid * (T_PN1 - T_PN2) - hApn * (T_PN1 - T_TN1) + P_decay.Q * m_PN2;
      m_PN3 * cP_pFluid * der(T_PN3) = m_dot_pFluid * cP_pFluid * (T_PN2 - T_PN3) - hApn * (T_PN3 - T_TN2) + P_decay.Q * m_PN3;
      m_PN4 * cP_pFluid * der(T_PN4) = m_dot_pFluid * cP_pFluid * (T_PN3 - T_PN4) - hApn * (T_PN3 - T_TN2) + P_decay.Q * m_PN4;
      m_TN1 * cP_Tube * der(T_TN1) = 2 * hApn * (T_PN1 - T_TN1) - 2 * hAsn * (T_TN1 - T_SN3);
      m_TN2 * cP_Tube * der(T_TN2) = 2 * hApn * (T_PN3 - T_TN2) - 2 * hAsn * (T_TN2 - T_SN1);
      m_SN1 * cP_sFluid * der(T_SN1) = m_dot_sFluid * cP_sFluid * (T_in_sFluid.T - T_SN1) + hAsn * (T_TN2 - T_SN1);
      m_SN2 * cP_sFluid * der(T_SN2) = m_dot_sFluid * cP_sFluid * (T_SN1 - T_SN2) + hAsn * (T_TN2 - T_SN1);
      m_SN3 * cP_sFluid * der(T_SN3) = m_dot_sFluid * cP_sFluid * (T_SN2 - T_SN3) + hAsn * (T_TN1 - T_SN3);
      m_SN4 * cP_sFluid * der(T_SN4) = m_dot_sFluid * cP_sFluid * (T_SN3 - T_SN4) + hAsn * (T_TN1 - T_SN3);
//T_out_pFluid.T = delay(T_PN4,varPHXtoPipe);
//T_out_sFluid.T = delay(T_SN4,varPHXtoUHX);
      T_out_pFluid.T = delay(T_PN4, varTauPHXtoPipe, 5000);
      T_out_sFluid.T = delay(T_SN4, varTauPHXtoUHX, 5000);
      annotation(
        Diagram(graphics = {Rectangle(origin = {30.2864, 0.0340212}, extent = {{-149.242, 59.993}, {149.242, -59.993}}), Rectangle(origin = {-60, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {120, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {0, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-59.95, -29.61}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.63, -9.61}, {19.63, 9.61}}), Rectangle(origin = {-0.12, -29.82}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.88, 9.82}, {19.88, -9.82}}), Rectangle(origin = {59.8, -29.85}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.8, 9.85}, {19.8, -9.85}}), Rectangle(origin = {119.85, -29.89}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.85, 9.89}, {19.85, -9.89}}), Rectangle(origin = {89.92, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {60, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-30.08, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Text(origin = {-100, 8}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_IN"), Text(origin = {-100, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_OUT"), Text(origin = {160, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_IN"), Text(origin = {162, 12}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_OUT"), Text(origin = {147, 50}, extent = {{59, -14}, {-59, 14}}, textString = "PHX")}, coordinateSystem(extent = {{-120, 80}, {220, -80}})),
        Icon(graphics = {Rectangle(origin = {30.2864, 0.0340212}, extent = {{-149.242, 59.993}, {149.242, -59.993}}), Rectangle(origin = {-30.08, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {60, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {89.92, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {-0.12, -29.82}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.88, 9.82}, {19.88, -9.82}}), Rectangle(origin = {0, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-60, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {120, 30}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {59.8, -29.85}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.8, 9.85}, {19.8, -9.85}}), Rectangle(origin = {119.85, -29.89}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.85, 9.89}, {19.85, -9.89}}), Rectangle(origin = {-59.95, -29.61}, lineColor = {239, 41, 41}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-19.63, -9.61}, {19.63, 9.61}}), Text(origin = {-100, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_OUT"), Text(origin = {162, 12}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_OUT"), Text(origin = {-100, 8}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_IN"), Text(origin = {160, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_IN"), Text(origin = {147, 50}, extent = {{59, -14}, {-59, 14}}, textString = "PHX")}));
    end PrimaryHeatExchanger;

    model SecondaryHeatExchanger
      // Parameter decleration
      parameter SMD_MSR_Modelica.Units.Mass m_PN1 = 1;
      parameter SMD_MSR_Modelica.Units.Mass m_PN2 = 1;
      parameter SMD_MSR_Modelica.Units.Mass m_PN3 = 1;
      parameter SMD_MSR_Modelica.Units.Mass m_PN4 = 1;
      parameter SMD_MSR_Modelica.Units.Mass m_TN1 = 1;
      parameter SMD_MSR_Modelica.Units.Mass m_TN2 = 1;
      parameter SMD_MSR_Modelica.Units.Mass m_SN1 = 1;
      parameter SMD_MSR_Modelica.Units.Mass m_SN2 = 1;
      parameter SMD_MSR_Modelica.Units.Mass m_SN3 = 1;
      parameter SMD_MSR_Modelica.Units.Mass m_SN4 = 1;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_pFluid = 20;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_Tube = 10;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_sFluid = 20;
      parameter SMD_MSR_Modelica.Units.MassFlowRate m_dot_pFluid = 10;
      parameter SMD_MSR_Modelica.Units.MassFlowRate m_dot_sFluid = 10;
      parameter SMD_MSR_Modelica.Units.Convection hA = 2;
      parameter SMD_MSR_Modelica.Units.Temperature T_in_pFluid = 400;
      parameter SMD_MSR_Modelica.Units.Temperature T_in_sFluid = 300;
      // Variable decleration
      SMD_MSR_Modelica.Units.Temperature T_PN1;
      SMD_MSR_Modelica.Units.Temperature T_PN2;
      SMD_MSR_Modelica.Units.Temperature T_PN3;
      SMD_MSR_Modelica.Units.Temperature T_PN4;
      SMD_MSR_Modelica.Units.Temperature T_TN1;
      SMD_MSR_Modelica.Units.Temperature T_TN2;
      SMD_MSR_Modelica.Units.Temperature T_SN1;
      SMD_MSR_Modelica.Units.Temperature T_SN2;
      SMD_MSR_Modelica.Units.Temperature T_SN3;
      SMD_MSR_Modelica.Units.Temperature T_SN4;
    initial equation
      T_PN1 = 10;
      T_PN2 = 10;
      T_PN3 = 10;
      T_PN4 = 10;
      T_TN1 = 10;
      T_TN2 = 10;
      T_SN1 = 10;
      T_SN2 = 10;
      T_SN3 = 10;
      T_SN4 = 10;
    equation
      m_PN1 * cP_pFluid * der(T_PN1) = m_dot_pFluid * cP_pFluid * (T_in_pFluid - T_PN1) - hA * (T_PN1 - T_TN1);
      m_PN2 * cP_pFluid * der(T_PN2) = m_dot_pFluid * cP_pFluid * (T_PN1 - T_PN2) - hA * (T_PN1 - T_TN1);
      m_PN3 * cP_pFluid * der(T_PN3) = m_dot_pFluid * cP_pFluid * (T_PN2 - T_PN3) - hA * (T_PN3 - T_TN2);
      m_PN4 * cP_pFluid * der(T_PN4) = m_dot_pFluid * cP_pFluid * (T_PN3 - T_PN4) - hA * (T_PN3 - T_TN2);
      m_TN1 * cP_Tube * der(T_TN1) = hA * (T_PN1 - T_TN1) + hA * (T_PN1 - T_TN1) - hA * (T_TN1 - T_SN3) - hA * (T_TN1 - T_SN3);
      m_TN2 * cP_Tube * der(T_TN2) = hA * (T_PN3 - T_TN2) + hA * (T_PN3 - T_TN2) - hA * (T_TN2 - T_SN1) - hA * (T_TN2 - T_SN1);
      m_SN1 * cP_sFluid * der(T_SN1) = m_dot_sFluid * cP_sFluid * (T_in_sFluid - T_SN1) + hA * (T_TN2 - T_SN1);
      m_SN2 * cP_sFluid * der(T_SN2) = m_dot_sFluid * cP_sFluid * (T_SN1 - T_SN2) - hA * (T_TN2 - T_SN1);
      m_SN3 * cP_sFluid * der(T_SN3) = m_dot_sFluid * cP_sFluid * (T_SN2 - T_SN3) - hA * (T_TN1 - T_SN3);
      m_SN4 * cP_sFluid * der(T_SN4) = m_dot_sFluid * cP_sFluid * (T_SN3 - T_SN4) - hA * (T_TN1 - T_SN3);
    end SecondaryHeatExchanger;

    model Radiator
      parameter SMD_MSR_Modelica.Units.MassFlowRate m_dot_air;
      parameter SMD_MSR_Modelica.Units.MassFlowRate m_dot_saltNom;
      parameter SMD_MSR_Modelica.Units.Convection hA_radiatorNom;
      parameter SMD_MSR_Modelica.Units.Mass m_airN;
      parameter SMD_MSR_Modelica.Units.Mass m_saltN;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity Cp_salt;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity Cp_air;
      parameter SMD_MSR_Modelica.Units.Temperature AirTempIn;
      parameter SMD_MSR_Modelica.Units.ResidentTime tauUHXtoPHX;
      parameter SMD_MSR_Modelica.Units.InitiationTime tripRadiator;
      SMD_MSR_Modelica.Units.ResidentTime varTauUHXtoPHX;
      //SMD_MSR_Modelica.Units.ResidentTime maxTauUHXtoPHX;
      SMD_MSR_Modelica.Units.Convection hA_radiator;
      SMD_MSR_Modelica.Units.Temperature Rad_TempOUT;
      SMD_MSR_Modelica.Units.Temperature AirTempOut;
      SMD_MSR_Modelica.Units.Demand PowerDemand;
      SMD_MSR_Modelica.Units.MassFlowRate m_dot_salt;
      input SMD_MSR_Modelica.PortsConnectors.Temp_In RadTempIN annotation(
        Placement(visible = true, transformation(origin = {-68, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out RadTempOut annotation(
        Placement(visible = true, transformation(origin = {28, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {28, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In secondaryFlowFrac annotation(
        Placement(visible = true, transformation(origin = {-68, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      RadTempOut.T = 546.1100;
      AirTempOut = 148.9;
    equation
      m_dot_salt = m_dot_saltNom * secondaryFlowFrac.FF;
      hA_radiator = hA_radiatorNom * (0.8215 * secondaryFlowFrac.FF ^ 6 - 4.108 * secondaryFlowFrac.FF ^ 5 + 7.848 * secondaryFlowFrac.FF ^ 4 - 7.165 * secondaryFlowFrac.FF ^ 3 + 3.004 * secondaryFlowFrac.FF ^ 2 + 0.5903 * secondaryFlowFrac.FF + 0.008537);
//hA_radiator = hA_radiatorNom;
      varTauUHXtoPHX = tauUHXtoPHX / secondaryFlowFrac.FF;
//maxTauUHXtoPHX = tauUHXtoPHX/secondaryFlowFrac.FF_min;
      PowerDemand = if time < tripRadiator then 1 else 0;
      der(Rad_TempOUT) * m_saltN * Cp_salt = m_dot_salt * Cp_salt * (RadTempIN.T - Rad_TempOUT) - hA_radiator * (Rad_TempOUT - AirTempOut);
      der(AirTempOut) * m_airN * Cp_air = PowerDemand * m_dot_air * Cp_air * (AirTempIn - AirTempOut) + hA_radiator * (RadTempOut.T - AirTempOut);
      RadTempOut.T = delay(Rad_TempOUT, varTauUHXtoPHX, 5000);
      annotation(
        Diagram(graphics = {Rectangle(origin = {-19.3958, 14.7922}, extent = {{-63.2606, 33.1616}, {63.2606, -33.1616}}), Text(origin = {-63, 7}, extent = {{-15, 5}, {15, -5}}, textString = "Temp IN"), Text(origin = {23, 9}, extent = {{-15, 5}, {15, -5}}, textString = "Temp OUT"), Text(origin = {8, 37}, extent = {{-30, 7}, {30, -7}}, textString = "Radiator")}),
        Icon(graphics = {Text(origin = {-63, 7}, extent = {{-15, 5}, {15, -5}}, textString = "Temp IN"), Text(origin = {8, 37}, extent = {{-30, 7}, {30, -7}}, textString = "Radiator"), Rectangle(origin = {-19.3958, 14.7922}, extent = {{-63.2606, 33.1616}, {63.2606, -33.1616}}), Text(origin = {23, 9}, extent = {{-15, 5}, {15, -5}}, textString = "Temp OUT")}));
    end Radiator;

    model DHRS
      parameter SMD_MSR_Modelica.Units.Mass m_DHRS;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity Cp_fluid;
      parameter SMD_MSR_Modelica.Units.MassFlowRate m_dotNom;
      parameter SMD_MSR_Modelica.Units.Mass m_Pi_C_PHX;
      parameter SMD_MSR_Modelica.Units.DHRStimeConstant DHRS_timeConstant;
      parameter SMD_MSR_Modelica.Units.HeatFlowRate DHRS_MaxPowerRm;
      parameter SMD_MSR_Modelica.Units.HeatFlowRate DHRS_PowerBleed;
      parameter SMD_MSR_Modelica.Units.InitiationTime DHRS_time;
      parameter SMD_MSR_Modelica.Units.ResidentTime tauDHRStoPHX;
      SMD_MSR_Modelica.Units.ResidentTime varTauDHRStoPHX;
      //SMD_MSR_Modelica.Units.ResidentTime maxTauDHRStoPHX;
      SMD_MSR_Modelica.Units.MassFlowRate m_dot;
      SMD_MSR_Modelica.Units.Temperature DHRS_Temp;
      SMD_MSR_Modelica.Units.HeatFlowRate DHRS_PowerRm;
      input SMD_MSR_Modelica.PortsConnectors.Temp_In DHRS_tempIN annotation(
        Placement(visible = true, transformation(origin = {-60, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.DecayHeat_In DHRS_DecayHeat annotation(
        Placement(visible = true, transformation(origin = {-60, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out DHRS_TempOUT annotation(
        Placement(visible = true, transformation(origin = {44, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {44, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In fuelFlowFrac annotation(
        Placement(visible = true, transformation(origin = {-60, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      DHRS_TempOUT.T = 650.7500;
    equation
      m_dot = m_dotNom * fuelFlowFrac.FF;
      varTauDHRStoPHX = tauDHRStoPHX / fuelFlowFrac.FF;
//maxTauDHRStoPHX = tauDHRStoPHX/fuelFlowFrac.FF_min;
      DHRS_PowerRm = (DHRS_MaxPowerRm - DHRS_PowerBleed) / (1 + exp(log(1 / 1E-3 - 1) * (1 - (time - DHRS_time) / DHRS_timeConstant))) + DHRS_PowerBleed;
//DHRS_PowerRm = (DHRS_MaxPowerRm - DHRS_PowerBleed) / (1 + exp(1 - (time - DHRS_time) / DHRS_timeConstant * log(1 / 1E-3 - 1))) + DHRS_PowerBleed;
      m_DHRS * Cp_fluid * der(DHRS_Temp) = m_dot * Cp_fluid * (DHRS_tempIN.T - DHRS_Temp) + DHRS_DecayHeat.Q * m_Pi_C_PHX - DHRS_PowerRm;
      DHRS_TempOUT.T = delay(DHRS_Temp, varTauDHRStoPHX, 5000);
      annotation(
        Diagram(graphics = {Rectangle(origin = {-10, 20}, extent = {{-70, 60}, {70, -60}}), Text(origin = {-57, 5}, extent = {{-15, 11}, {15, -11}}, textString = "Temp IN"), Text(origin = {41, 5}, extent = {{-15, 11}, {15, -11}}, textString = "Temp OUT"), Text(origin = {-57, 47}, extent = {{-15, 11}, {15, -11}}, textString = "Decay Heat"), Text(origin = {31, 66}, extent = {{-31, 12}, {31, -12}}, textString = "DHRS")}),
        Icon(graphics = {Text(origin = {31, 66}, extent = {{-31, 12}, {31, -12}}, textString = "DHRS"), Rectangle(origin = {-10, 20}, extent = {{-70, 60}, {70, -60}}), Text(origin = {-57, 47}, extent = {{-15, 11}, {15, -11}}, textString = "Decay Heat"), Text(origin = {-57, 5}, extent = {{-15, 11}, {15, -11}}, textString = "Temp IN"), Text(origin = {41, 5}, extent = {{-15, 11}, {15, -11}}, textString = "Temp OUT")}));
    end DHRS;

    model Pipe
      parameter SMD_MSR_Modelica.Units.Mass m_pi;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity Cp_fluidPi;
      parameter SMD_MSR_Modelica.Units.MassFlowRate m_dotPiNom;
      parameter SMD_MSR_Modelica.Units.Mass m_Pi_PHX_C;
      parameter SMD_MSR_Modelica.Units.ResidentTime tauPiToCore;
      SMD_MSR_Modelica.Units.ResidentTime varTauPiToCore;
      //SMD_MSR_Modelica.Units.ResidentTime maxTauPiToCore;
      SMD_MSR_Modelica.Units.MassFlowRate m_dotPi;
      SMD_MSR_Modelica.Units.Temperature Pipe_Temp;
      input SMD_MSR_Modelica.PortsConnectors.Temp_In PiTemp_IN annotation(
        Placement(visible = true, transformation(origin = {-48, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-48, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.DecayHeat_In PiDecay_Heat annotation(
        Placement(visible = true, transformation(origin = {-48, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-48, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out PiTemp_Out annotation(
        Placement(visible = true, transformation(origin = {28, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {28, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In fuelFlowFrac annotation(
        Placement(visible = true, transformation(origin = {-48, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-48, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      PiTemp_Out.T = 632;
    equation
      m_dotPi = m_dotPiNom * fuelFlowFrac.FF;
      varTauPiToCore = tauPiToCore / fuelFlowFrac.FF;
//maxTauPiToCore = tauPiToCore/fuelFlowFrac.FF_min;
      m_pi * Cp_fluidPi * der(Pipe_Temp) = m_dotPi * Cp_fluidPi * (PiTemp_IN.T - Pipe_Temp) + PiDecay_Heat.Q * m_Pi_PHX_C;
      PiTemp_Out.T = delay(Pipe_Temp, varTauPiToCore, 5000);
      annotation(
        Diagram(graphics = {Rectangle(origin = {-10, 24}, extent = {{-54, 50}, {54, -50}}), Text(origin = {-46, 11}, extent = {{-14, 9}, {14, -9}}, textString = "Temp IN"), Text(origin = {26, 11}, extent = {{-14, 9}, {14, -9}}, textString = "Temp OUT"), Text(origin = {12, 66}, extent = {{-24, 30}, {24, -30}}, textString = "PipeNode"), Text(origin = {-46, 47}, extent = {{-14, 9}, {14, -9}}, textString = "Decay Heat")}),
        Icon(graphics = {Text(origin = {26, 11}, extent = {{-14, 9}, {14, -9}}, textString = "Temp OUT"), Text(origin = {12, 66}, extent = {{-24, 30}, {24, -30}}, textString = "PipeNode"), Text(origin = {-46, 11}, extent = {{-14, 9}, {14, -9}}, textString = "Temp IN"), Text(origin = {-46, 47}, extent = {{-14, 9}, {14, -9}}, textString = "Decay Heat"), Rectangle(origin = {-10, 24}, extent = {{-54, 50}, {54, -50}})}));
    end Pipe;

    model PrimaryPump
      parameter SMD_MSR_Modelica.Units.FlowFraction freeConvectionFF;
      parameter SMD_MSR_Modelica.Units.PumpConstant primaryPumpK;
      parameter SMD_MSR_Modelica.Units.InitiationTime tripPrimaryPump;
      output SMD_MSR_Modelica.PortsConnectors.FlowFraction_Out primaryFlowFrac annotation(
        Placement(visible = true, transformation(origin = {-1.77636e-15, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {1.77636e-15, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    initial equation
//primaryFlowFrac.FF = 1;
    equation
      primaryFlowFrac.FF = (1 - freeConvectionFF) * exp(-primaryPumpK * delay(time, tripPrimaryPump)) + freeConvectionFF;
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
      secondaryFlowFrac.FF = (1 - freeConvectionFF) * exp(-secondaryPumpK * delay(time, tripSecondaryPump)) + freeConvectionFF;
      annotation(
        Diagram(graphics = {Rectangle(lineColor = {245, 121, 0}, lineThickness = 1, extent = {{-60, 60}, {60, -60}}), Text(origin = {2, 3}, extent = {{-52, 27}, {52, -27}}, textString = "Secondary pump")}),
        Icon(graphics = {Rectangle(origin = {0, -20}, lineColor = {245, 121, 0}, lineThickness = 1, extent = {{-60, 60}, {60, -60}}), Text(origin = {2, -19}, extent = {{-52, 27}, {52, -27}}, textString = "Secondary pump")}));
    end SecondaryPump;

    model UHX
      parameter SMD_MSR_Modelica.Units.Mass m_PN1UHX;
      parameter SMD_MSR_Modelica.Units.ReactorPower UHXP;
      parameter SMD_MSR_Modelica.Units.SpecificHeatCapacity cP_pFluidUHX;
      parameter SMD_MSR_Modelica.Units.MassFlowRate m_dot_pFluidUHXnom;
      parameter SMD_MSR_Modelica.Units.ResidentTime tauUHXtoPHX;
      parameter SMD_MSR_Modelica.Units.InitiationTime tripUHX;
      parameter SMD_MSR_Modelica.Units.Demand SetDemand;
      SMD_MSR_Modelica.Units.MassFlowRate m_dot_pFluidUHX;
      SMD_MSR_Modelica.Units.ResidentTime varTauUHXtoPHX;
      SMD_MSR_Modelica.Units.Temperature UHXtemp;
      SMD_MSR_Modelica.Units.Demand PowerDemand;
      input SMD_MSR_Modelica.PortsConnectors.Temp_In UHXtemp_In annotation(
        Placement(visible = true, transformation(origin = {-70, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out UHXtemp_Out annotation(
        Placement(visible = true, transformation(origin = {26, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {26, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In secondaryFlowFrac annotation(
        Placement(visible = true, transformation(origin = {-68, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      UHXtemp_Out.T = 546.1100;
    equation
      PowerDemand = if time < tripUHX then 1 else SetDemand;
      varTauUHXtoPHX = tauUHXtoPHX / secondaryFlowFrac.FF;
      m_dot_pFluidUHX = m_dot_pFluidUHXnom * secondaryFlowFrac.FF;
      m_PN1UHX * cP_pFluidUHX * der(UHXtemp) = m_dot_pFluidUHX * cP_pFluidUHX * (UHXtemp_In.T - UHXtemp) - UHXP * PowerDemand;
      UHXtemp_Out.T = delay(UHXtemp, varTauUHXtoPHX, 5000);
      annotation(
        Diagram(graphics = {Rectangle(origin = {-20, 0}, extent = {{-60, 60}, {60, -60}}), Text(origin = {-64, -14}, extent = {{-14, 12}, {14, -12}}, textString = "Temp In"), Text(origin = {24, -14}, extent = {{-14, 12}, {14, -12}}, textString = "Temp Out"), Text(origin = {8, 44}, extent = {{-34, 18}, {34, -18}}, textString = "UHX")}),
        Icon(graphics = {Text(origin = {8, 44}, extent = {{-34, 18}, {34, -18}}, textString = "UHX"), Rectangle(origin = {-20, 0}, extent = {{-60, 60}, {60, -60}}), Text(origin = {24, -14}, extent = {{-14, 12}, {14, -12}}, textString = "Temp Out"), Text(origin = {-64, -14}, extent = {{-14, 12}, {14, -12}}, textString = "Temp In")}));
    end UHX;
  end HeatTransport;

  package Nuclear
    model mPKE
      //Parameter decleration
      parameter SMD_MSR_Modelica.Units.NeutronGenerationTime Lam;
      parameter SMD_MSR_Modelica.Units.PrecursorDecayConstant lam[6];
      parameter SMD_MSR_Modelica.Units.DelayedNeutronFrac beta[6];
      parameter SMD_MSR_Modelica.Units.ResidentTime tauCoreNom;
      parameter SMD_MSR_Modelica.Units.ResidentTime tauLoopNom;
      //Variable decleration
      SMD_MSR_Modelica.Units.ResidentTime varTauCore;
      //SMD_MSR_Modelica.Units.ResidentTime maxTauCore;
      SMD_MSR_Modelica.Units.ResidentTime varTauLoop;
      //SMD_MSR_Modelica.Units.ResidentTime maxTauLoop;
      SMD_MSR_Modelica.Units.PrecursorConc CG1;
      SMD_MSR_Modelica.Units.PrecursorConc CG2;
      SMD_MSR_Modelica.Units.PrecursorConc CG3;
      SMD_MSR_Modelica.Units.PrecursorConc CG4;
      SMD_MSR_Modelica.Units.PrecursorConc CG5;
      SMD_MSR_Modelica.Units.PrecursorConc CG6;
      SMD_MSR_Modelica.Units.PrecursorReturn CG1Return;
      SMD_MSR_Modelica.Units.PrecursorReturn CG2Return;
      SMD_MSR_Modelica.Units.PrecursorReturn CG3Return;
      SMD_MSR_Modelica.Units.PrecursorReturn CG4Return;
      SMD_MSR_Modelica.Units.PrecursorReturn CG5Return;
      SMD_MSR_Modelica.Units.PrecursorReturn CG6Return;
      SMD_MSR_Modelica.Units.Reactivity reactivity;
      input SMD_MSR_Modelica.PortsConnectors.Reactivity feedback annotation(
        Placement(visible = true, transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.NominalNeutronPopulation n_population annotation(
        Placement(visible = true, transformation(origin = {0, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In fuelFlowFrac annotation(
        Placement(visible = true, transformation(origin = {34, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {34, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      n_population.n = 1;
      CG1 = (beta[1] / Lam) * (1.0 / (lam[1] - (exp(-lam[1] * tauLoopNom) - 1.0) / tauCoreNom));
      CG2 = (beta[2] / Lam) * (1.0 / (lam[2] - (exp(-lam[2] * tauLoopNom) - 1.0) / tauCoreNom));
      CG3 = (beta[3] / Lam) * (1.0 / (lam[3] - (exp(-lam[3] * tauLoopNom) - 1.0) / tauCoreNom));
      CG4 = (beta[4] / Lam) * (1.0 / (lam[4] - (exp(-lam[4] * tauLoopNom) - 1.0) / tauCoreNom));
      CG5 = (beta[5] / Lam) * (1.0 / (lam[5] - (exp(-lam[5] * tauLoopNom) - 1.0) / tauCoreNom));
      CG6 = (beta[6] / Lam) * (1.0 / (lam[6] - (exp(-lam[6] * tauLoopNom) - 1.0) / tauCoreNom));
    equation
      reactivity = feedback.rho + 0.002465140767843;
      der(n_population.n) = (reactivity-sum(beta))/Lam*n_population.n+lam[1]*CG1+lam[2]*CG2+lam[3]*CG3+lam[4]*CG4+lam[5]*CG5+lam[6]*CG6;
      varTauCore = tauCoreNom / fuelFlowFrac.FF;
      varTauLoop = tauLoopNom / fuelFlowFrac.FF;
//maxTauCore = tauCoreNom/fuelFlowFrac.FF_min;
//maxTauLoop = tauLoopNom/fuelFlowFrac.FF_min;
      CG1Return = delay(CG1,varTauLoop,5000)*exp(-lam[1]*varTauLoop)/varTauCore;
      CG2Return = delay(CG2,varTauLoop,5000)*exp(-lam[2]*varTauLoop)/varTauCore;
      CG3Return = delay(CG3,varTauLoop,5000)*exp(-lam[3]*varTauLoop)/varTauCore;
      CG4Return = delay(CG4,varTauLoop,5000)*exp(-lam[4]*varTauLoop)/varTauCore;
      CG5Return = delay(CG5,varTauLoop,5000)*exp(-lam[5]*varTauLoop)/varTauCore;
      CG6Return = delay(CG6,varTauLoop,5000)*exp(-lam[6]*varTauLoop)/varTauCore;
      der(CG1) = (beta[1]*n_population.n)/Lam - (lam[1]*CG1) - (CG1/varTauCore) + CG1Return;
      der(CG2) = (beta[2]*n_population.n)/Lam - (lam[2]*CG2) - (CG2/varTauCore) + CG2Return;
      der(CG3) = (beta[3]*n_population.n)/Lam - (lam[3]*CG3) - (CG3/varTauCore) + CG3Return;
      der(CG4) = (beta[4]*n_population.n)/Lam - (lam[4]*CG4) - (CG4/varTauCore) + CG4Return;
      der(CG5) = (beta[5]*n_population.n)/Lam - (lam[5]*CG5) - (CG5/varTauCore) + CG5Return;
      der(CG6) = (beta[6]*n_population.n)/Lam - (lam[6]*CG6) - (CG6/varTauCore) + CG6Return;
      annotation(
        Diagram(graphics = {Rectangle(origin = {2, 0}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, lineThickness = 0.5, extent = {{-50, 50}, {50, -50}}), Text(origin = {3, 1}, extent = {{37, -27}, {-37, 27}}, textString = "mPKE")}),
        Icon(graphics = {Rectangle(origin = {2, 0}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, lineThickness = 0.5, extent = {{-50, 50}, {50, -50}}), Text(origin = {3, 1}, extent = {{37, -27}, {-37, 27}}, textString = "mPKE")}));
    end mPKE;

    model DecayHeat
      parameter SMD_MSR_Modelica.Units.ReactorPower P;
      parameter SMD_MSR_Modelica.Units.DecayHeatYield DHYG1;
      parameter SMD_MSR_Modelica.Units.DecayHeatYield DHYG2;
      parameter SMD_MSR_Modelica.Units.DecayHeatYield DHYG3;
      parameter SMD_MSR_Modelica.Units.Mass TotalFuelMass;
      parameter SMD_MSR_Modelica.Units.DecayHeatPrecursorDecayConstant DHlamG1;
      parameter SMD_MSR_Modelica.Units.DecayHeatPrecursorDecayConstant DHlamG2;
      parameter SMD_MSR_Modelica.Units.DecayHeatPrecursorDecayConstant DHlamG3;
      SMD_MSR_Modelica.Units.HeatTransferFraction DHG1;
      SMD_MSR_Modelica.Units.HeatTransferFraction DHG2;
      SMD_MSR_Modelica.Units.HeatTransferFraction DHG3;
      input SMD_MSR_Modelica.PortsConnectors.NominalNeutronPopulation nPop annotation(
        Placement(visible = true, transformation(origin = {1, 33}, extent = {{-19, -19}, {19, 19}}, rotation = 0), iconTransformation(origin = {-32, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.DecayHeat_Out decayHeat_Out annotation(
        Placement(visible = true, transformation(origin = {1, -31}, extent = {{-23, -23}, {23, 23}}, rotation = 0), iconTransformation(origin = {36, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      DHG1 = 0.025125;
      DHG2 = 0.0198578;
      DHG3 = 0.0227606;
    equation
      der(DHG1) = nPop.n * DHYG1 - DHlamG1 * DHG1;
      der(DHG2) = nPop.n * DHYG2 - DHlamG2 * DHG2;
      der(DHG3) = nPop.n * DHYG3 - DHlamG3 * DHG3;
      decayHeat_Out.Q = (DHG1 + DHG2 + DHG3) * P / TotalFuelMass;
//decayHeat_Out.Q = 0.066*P/TotalFuelMass;
//decayHeat_Out.Q = 0;
      annotation(
        Diagram(graphics = {Rectangle(origin = {2, 0}, lineColor = {0, 225, 255}, lineThickness = 0.75, extent = {{-50, 50}, {50, -50}}), Text(origin = {3, 0}, extent = {{-45, 46}, {45, -46}}, textString = "Decay Heat")}),
        Icon(graphics = {Text(origin = {19, 41}, extent = {{-29, 33}, {29, -33}}, textString = "Decay Heat"), Rectangle(origin = {2, 0}, lineColor = {0, 225, 255}, lineThickness = 0.75, extent = {{-50, 50}, {50, -50}}), Text(origin = {1, 36}, extent = {{1, 2}, {-1, -2}}, textString = "text"), Text(origin = {-31, 9}, extent = {{-13, 5}, {13, -5}}, textString = "n(t)/n_0"), Text(origin = {19, -43}, extent = {{-31, 5}, {31, -5}}, textString = "Specific Decay Heat")}));
    end DecayHeat;

    model ReactivityFeedback
      parameter SMD_MSR_Modelica.Units.Temperature FuelTempSetPointNode1;
      parameter SMD_MSR_Modelica.Units.Temperature FuelTempSetPointNode2;
      parameter SMD_MSR_Modelica.Units.Temperature GrapTempSetPoint;
      parameter SMD_MSR_Modelica.Units.TemperatureReactivityCoef a_F;
      parameter SMD_MSR_Modelica.Units.TemperatureReactivityCoef a_G;
      parameter SMD_MSR_Modelica.Units.Reactivity step_mag;
      parameter SMD_MSR_Modelica.Units.Frequency omega;
      parameter SMD_MSR_Modelica.Units.Reactivity sin_mag;
      parameter SMD_MSR_Modelica.Units.InitiationTime stepInsertionTime;
      parameter SMD_MSR_Modelica.Units.InitiationTime sinInsertionTime;
      SMD_MSR_Modelica.Units.Reactivity FuelTempFeedbackNode1;
      SMD_MSR_Modelica.Units.Reactivity FuelTempFeedbackNode2;
      SMD_MSR_Modelica.Units.Reactivity GrapTempFeedback;
      SMD_MSR_Modelica.Units.Reactivity TotalTempFeedback;
      SMD_MSR_Modelica.Units.Reactivity ReactExternalStep;
      SMD_MSR_Modelica.Units.Reactivity ReactExternalSin;
      input SMD_MSR_Modelica.PortsConnectors.Temp_In fuelNode1 annotation(
        Placement(visible = true, transformation(origin = {-40, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.Temp_In fuelNode2 annotation(
        Placement(visible = true, transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.Temp_In grapNode annotation(
        Placement(visible = true, transformation(origin = {-40, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Reactivity feedback annotation(
        Placement(visible = true, transformation(origin = {32, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {34, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      ReactExternalStep = 0 + (if time < stepInsertionTime then 0 else step_mag);
      ReactExternalSin = 0 + (if time < sinInsertionTime then 0 else sin_mag * sin(omega * time));
      FuelTempFeedbackNode1 = (fuelNode1.T - FuelTempSetPointNode1) * (a_F / 2);
      FuelTempFeedbackNode2 = (fuelNode2.T - FuelTempSetPointNode2) * (a_F / 2);
      GrapTempFeedback = (grapNode.T - GrapTempSetPoint) * a_G;
      TotalTempFeedback = FuelTempFeedbackNode1 + FuelTempFeedbackNode2 + GrapTempFeedback + ReactExternalSin + ReactExternalStep;
      feedback.rho = TotalTempFeedback;
      annotation(
        Diagram(graphics = {Rectangle(origin = {0, 2}, extent = {{-50, 50}, {50, -50}}), Text(origin = {13, 40}, extent = {{-31, 24}, {31, -24}}, textString = "Feedback"), Text(origin = {-32, 22}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN1"), Text(origin = {-32, -12}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN2"), Text(origin = {-32, -40}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_GrapN1"), Text(origin = {24, 4}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_Feedback")}),
        Icon(graphics = {Rectangle(origin = {0, 2}, extent = {{-50, 50}, {50, -50}}), Text(origin = {-32, 22}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN1"), Text(origin = {13, 40}, extent = {{-31, 24}, {31, -24}}, textString = "Feedback"), Text(origin = {-32, -12}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN2"), Rectangle(origin = {0, 2}, extent = {{-50, 50}, {50, -50}}), Text(origin = {-32, -40}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_GrapN1"), Text(origin = {30, 6}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_Feedback")}));
    end ReactivityFeedback;

    model FuelChannel
    
      parameter SMD_MSR_Modelica.Units.ReactorPower P;
      parameter SMD_MSR_Modelica.Units.Mass TotalFuelMass;
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
      parameter SMD_MSR_Modelica.Units.ResidentTime tauCoreToDHRS;
    
      SMD_MSR_Modelica.Units.ResidentTime varTauCoreToDHRS;
      //SMD_MSR_Modelica.Units.ResidentTime maxTauCoreToDHRS;
      SMD_MSR_Modelica.Units.NominalPower FissionPower;
      SMD_MSR_Modelica.Units.NominalPower DecayPower;
      SMD_MSR_Modelica.Units.NominalPower NomPower;
      SMD_MSR_Modelica.Units.MassFlowRate mdot_fuel;
      SMD_MSR_Modelica.Units.Convection hA;
    
      input SMD_MSR_Modelica.PortsConnectors.Temp_In temp_In annotation(
        Placement(visible = true, transformation(origin = {-76, -70}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-80, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.DecayHeat_In P_decay annotation(
        Placement(visible = true, transformation(origin = {-77, -7}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {-80, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.NominalNeutronPopulation nPop annotation(
        Placement(visible = true, transformation(origin = {-80, 50}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-80, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out temp_Out annotation(
        Placement(visible = true, transformation(origin = {44, 60}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {50, 53}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out fuelNode1 annotation(
        Placement(visible = true, transformation(origin = {45, -78}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {50, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out fuelNode2 annotation(
        Placement(visible = true, transformation(origin = {45, -40}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {50, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output SMD_MSR_Modelica.PortsConnectors.Temp_Out grapNode annotation(
        Placement(visible = true, transformation(origin = {44, 16}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {50, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input SMD_MSR_Modelica.PortsConnectors.FlowFraction_In fuelFlowFraction annotation(
        Placement(visible = true, transformation(origin = {-40, -80}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-40, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
    initial equation
      fuelNode2.T = 657;
      fuelNode1.T = 644.5000;
      grapNode.T = 675.6111;
    
    equation
      FissionPower = nPop.n * (1 - 0.068);
      DecayPower = P_decay.Q / P * TotalFuelMass;
      NomPower = FissionPower + DecayPower;
      mdot_fuel = mdot_fuelNom * fuelFlowFraction.FF;
    
      hA = hAnom * (0.8215 * fuelFlowFraction.FF ^ 6 - 4.108 * fuelFlowFraction.FF ^ 5 + 7.848 * fuelFlowFraction.FF ^ 4 - 7.165 * fuelFlowFraction.FF ^ 3 + 3.004 * fuelFlowFraction.FF ^ 2 + 0.5903 * fuelFlowFraction.FF + 0.008537);
//hA = hAnom;
      varTauCoreToDHRS = tauCoreToDHRS / fuelFlowFraction.FF;
//maxTauCoreToDHRS = tauCoreToDHRS/fuelFlowFraction.FF_min;
    
      m_FN1 * cP_fuel * der(fuelNode1.T) = mdot_fuel * cP_fuel * (temp_In.T - fuelNode1.T) + kFN1 * FissionPower * P - hA * kHT_FN1 * (fuelNode1.T - grapNode.T) + P_decay.Q * m_FN1;
      m_FN2 * cP_fuel * der(fuelNode2.T) = mdot_fuel * cP_fuel * (fuelNode1.T - fuelNode2.T) + kFN2 * FissionPower * P - hA * kHT_FN2 * (fuelNode1.T - grapNode.T) + P_decay.Q * m_FN2;
      m_GN * cP_graphite * der(grapNode.T) = hA * (fuelNode1.T - grapNode.T) + kG * FissionPower * P;
      temp_Out.T = delay(fuelNode2.T, varTauCoreToDHRS, 5000);
    
      annotation(
        Diagram(graphics = {Rectangle(origin = {-15.4297, -10.014}, extent = {{-75.8099, 90.0986}, {75.8099, -90.0986}}), Rectangle(origin = {-40, 20}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}}), Rectangle(origin = {-40, -40}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}}), Rectangle(origin = {15, -10}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-15, 50}, {15, -50}}), Text(origin = {10, 74}, extent = {{-47, 12}, {47, -12}}, textString = "Core Heat Transfer"), Text(origin = {-76, 35}, extent = {{-13, 13}, {13, -13}}, textString = "n(t)/n_0"), Text(origin = {-76, -20}, extent = {{-13, 13}, {13, -13}}, textString = "Decay Heat"), Text(origin = {-75, -85}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_In"), Text(origin = {46, 45}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_Out"), Text(origin = {47, 2}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_G1"), Text(origin = {45, -54}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_FN2"), Text(origin = {45, -92}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_FN1")}),
        Icon(graphics = {Rectangle(origin = {-40, 20}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}}), Rectangle(origin = {15, -10}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-15, 50}, {15, -50}}), Rectangle(origin = {-40, -40}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}}), Rectangle(origin = {-15.4297, -10.014}, extent = {{-75.8099, 90.0986}, {75.8099, -90.0986}}), Text(origin = {9, 72}, extent = {{-47, 12}, {47, -12}}, textString = "Core Heat Transfer"), Text(origin = {-76, -20}, extent = {{-13, 13}, {13, -13}}, textString = "Decay Heat"), Text(origin = {46, 4}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_G1"), Text(origin = {46, -84}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_FN1"), Text(origin = {46, -42}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_FN2"), Text(origin = {-76, -85}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_In"), Text(origin = {46, 41}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_Out"), Text(origin = {-76, 45}, extent = {{-13, 13}, {13, -13}}, textString = "n(t)/n_0")}));
    end FuelChannel;

    model Poisons
      //Design specs
      //MaterialProperties - Locate to MaterialProperties Library
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
      der(Te135_Conc) = Te135_gamma * Sig_fission * phi0;
    end Poisons;
  end Nuclear;

  package Units
    type Frequency = Real(unit = "rad/s");
    //Frequency [rad/s]
    // Nuclear - PKE
    type NeutronDensity = Real(unit = "1/m^3", min = 0);
    //Neutron Density [#n/m^3]
    type NominalNeutronPopulation = Real(unit = "1", min = 0);
    //Nominl Neutron Population
    type NeutronGenerationTime = Real(unit = "s", min = 0);
    //Neutron Generation Time (Lambda) [s]
    type Reactivity = Real(unit = "1");
    //Absolute Reactivity [unitless]
    type PrecursorDecayConstant = Real(unit = "1/s", min = 0);
    //Delayed Precursor Decay Constant [s^-1]
    type PrecursorConc = Real(unit = "1", min = 0);
    //Delayed Neutron Precursor
    type PrecursorReturn = Real(unit = "1/s", min = 0);
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
    type NominalNeutronFlux = Real(unit = "n/(cm^2.s)");
    // Flow
    type MassFlowRate = Real(unit = "kg/s", min = 0);
    //Mass Flow Rate [kg/m^3]
    type ResidentTime = Real(unit = "s", min = 0);
    //Fuel Salt Resident Time (tau) [s]
    type FlowFraction = Real(unit = "1", min = 0);
    //Flow Fraction [unitless]
    type SpecificHeat = Real(unit = "MJ/(kg.s)");
    type HeatFlowRate = Real(unit = "MJ/s");
    //Heat Flow Rate [MJ/kg.s]
    // Heat Transfer
    type NominalPower = Real(unit = "1", min = 0);
    //Nominal reactor power as a fraction
    type ReactorPower = Real(unit = "MW", min = 0);
    //Nominal Reactor Thermal Power
    type Convection = Real(unit = "MJ/(s.C)", min = 0);
    //Rate of heat transfer by convection [J/(s.C) == W/C]
    type HeatTransferFraction = Real(unit = "1", min = 0);
    //Fraction of Heat Trasfer from Node to an adjacent Node [unitless]
    // Decay Heat
    type DecayHeatPrecursorDecayConstant = Real(unit = "1/s", min = 0);
    //Decay heat precoursor decay constant [s^-1]
    type DecayHeatYield = Real(unit = "1/s", min = 0);
    type DecayHeatFraction = Real(tinits = "1", min = 0);
    // Physical Properties
    type Mass = Real(unit = "kg", min = 0);
    type SpecificHeatCapacity = Real(unit = "MJ/(kg.C)", min = 0);
    type Temperature = Real(unit = "C", min = 0);
    // Component Related
    type DHRStimeConstant = Real(unit = "1/s", min = 0);
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
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}}), Ellipse(origin = {-8, 36}, extent = {{0, -2}, {0, 2}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}));
    end Temp_In;

    connector Temp_Out
      SMD_MSR_Modelica.Units.Temperature T;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end Temp_Out;

    connector DecayHeat_In
      SMD_MSR_Modelica.Units.SpecificHeat Q;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end DecayHeat_In;

    connector DecayHeat_Out
      SMD_MSR_Modelica.Units.SpecificHeat Q;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end DecayHeat_Out;

    connector NominalNeutronPopulation
      SMD_MSR_Modelica.Units.NominalNeutronPopulation n;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}));
    end NominalNeutronPopulation;

    connector Reactivity
      SMD_MSR_Modelica.Units.Reactivity rho;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end Reactivity;

    connector ExternalReact
      SMD_MSR_Modelica.Units.Reactivity rho_Ex;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end ExternalReact;

    connector FlowFraction_In
      SMD_MSR_Modelica.Units.FlowFraction FF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}));
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end FlowFraction_In;

    connector FlowFraction_Out
      SMD_MSR_Modelica.Units.FlowFraction FF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {245, 121, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {245, 121, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end FlowFraction_Out;
  end PortsConnectors;

  package QAtoolBox
    package TestIO
      model NeutronPopulation
        parameter SMD_MSR_Modelica.Units.NominalNeutronPopulation relNpop;
        output SMD_MSR_Modelica.PortsConnectors.NominalNeutronPopulation n_population annotation(
          Placement(visible = true, transformation(origin = {-1, 29}, extent = {{-29, -29}, {29, 29}}, rotation = 0), iconTransformation(origin = {4, 38}, extent = {{-28, -28}, {28, 28}}, rotation = 0)));
      equation
        n_population.n = relNpop;
        annotation(
          Diagram(graphics = {Rectangle(origin = {0, -2}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, extent = {{-60, 60}, {60, -60}}), Text(origin = {3, -28}, extent = {{-51, 16}, {51, -16}}, textString = "n")}),
          Icon(graphics = {Text(origin = {7, -35}, extent = {{-53, 21}, {53, -21}}, textString = "n"), Rectangle(origin = {4, -6}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, lineThickness = 2, extent = {{-76, 72}, {76, -72}})}));
      end NeutronPopulation;

      model SteadyStateTempFeedback
        parameter SMD_MSR_Modelica.Units.Reactivity Feedback;
        SMD_MSR_Modelica.PortsConnectors.Reactivity reactivity annotation(
          Placement(visible = true, transformation(origin = {0, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 28}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
      equation
        reactivity.rho = Feedback;
        annotation(
          Diagram(graphics = {Rectangle(lineColor = {78, 154, 6}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {2, -36}, extent = {{-40, 20}, {40, -20}}, textString = "rho_{tempFB}")}),
          Icon(graphics = {Rectangle(lineColor = {78, 154, 6}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {2, -36}, extent = {{-40, 20}, {40, -20}}, textString = "rho_{tempFB}")}));
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
      equation
        flowFraction_Out.FF = ConstantFlowFrac;
        annotation(
          Diagram(graphics = {Rectangle(lineColor = {245, 121, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {3, -32}, extent = {{-45, 10}, {45, -10}}, textString = "FlowFrac")}),
          Icon(graphics = {Rectangle(lineColor = {245, 121, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {3, -32}, extent = {{-45, 10}, {45, -10}}, textString = "FlowFrac")}));
      end ConstantFlowFrac;

      model ConstantDecayPower
        parameter SMD_MSR_Modelica.Units.SpecificHeat DecayPower = 0;
        output SMD_MSR_Modelica.PortsConnectors.DecayHeat_Out decayHeat_Out annotation(
          Placement(visible = true, transformation(origin = {2, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        decayHeat_Out.Q = DecayPower;
        annotation(
          Diagram(graphics = {Rectangle(lineColor = {0, 225, 255}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {-2, -33}, extent = {{-42, 13}, {42, -13}}, textString = "P_{Decay}")}),
          Icon(graphics = {Rectangle(lineColor = {0, 225, 255}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {2, -35}, extent = {{-42, 13}, {42, -13}}, textString = "P_{Decay}")}));
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
        
        input SMD_MSR_Modelica.PortsConnectors.Temp_In temp_In annotation(
          Placement(visible = true, transformation(origin = {-40, -1.77636e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-40, -1.9984e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        output SMD_MSR_Modelica.PortsConnectors.Temp_Out temp_Out annotation(
          Placement(visible = true, transformation(origin = {40, -1.77636e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {40, -1.77636e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      
      equation
      
        nodeMass*cP*der(temp_Out.T) = mDot*cP*(temp_In.T - temp_Out.T) - qDotOut;
      
      annotation(
          Diagram(graphics = {Rectangle(lineColor = {204, 0, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}}), Text(origin = {0, -42}, extent = {{-48, -30}, {48, 30}}, textString = "Q Remove")}),
          Icon(graphics = {Text(origin = {0, -42}, extent = {{-48, -30}, {48, 30}}, textString = "Q Remove"), Rectangle(lineColor = {204, 0, 0}, lineThickness = 2, extent = {{-60, 60}, {60, -60}})}));
      
      end ConstantHeatRemoval;
    end TestIO;

    package TestCases
      model TestPHX
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac primaryFF(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {-67, 71}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac secondaryFF(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {67, -67}, extent = {{-27, -27}, {27, 27}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0) annotation(
          Placement(visible = true, transformation(origin = {92, 70}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
        SMD_MSR_Modelica.HeatTransport.PrimaryHeatExchanger primaryHeatExchanger(m_PN1 = 87.0959, m_PN2 = 87.0959, m_PN3 = 87.0959, m_PN4 = 87.0959, m_TN1 = 101.1662, m_TN2 = 101.1662, m_SN1 = 24.9200, m_SN2 = 24.9200, m_SN3 = 24.9200, m_SN4 = 24.9200, cP_pFluid = 1.9665E-3, cP_Tube = 5.778E-04, cP_sFluid = 2.39E-3, m_dot_pFluidNom = 7.5708E-02 * 2.14647E+03, m_dot_sFluidNom = 5.36265E-02 * 1.922E3, hApnNom = 0.1620, hAsnNom = 0.0765, tauPHXtoPipe = 3.85, tauPHXtoUHX = 4.71) annotation(
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
        SMD_MSR_Modelica.Nuclear.mPKE mPKE(Lam = 2.400E-04, lam = {1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00}, beta = {0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023}, tauCoreNom = 8.46, tauLoopNom = 16.73) annotation(
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
        SMD_MSR_Modelica.Nuclear.FuelChannel fuelChannel(P = 8, TotalFuelMass = 4.093499709644400e+03, hAnom = 0.0180, m_FN1 = 687.3959, m_FN2 = 687.3959, m_GN = 3.6342e+03, cP_fuel = 1.9665E-3, cP_graphite = 1.773E-3, mdot_fuelNom = 7.5708E-02 * 2.14647E+03, kFN1 = 0.465, kFN2 = 0.465, kG = 0.07, kHT_FN1 = 0.5, kHT_FN2 = 0.5, tauCoreToDHRS = 1.385) annotation(
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
        SMD_MSR_Modelica.HeatTransport.Radiator radiator(AirTempIn = 37.78, Cp_air = 1.0085E-3, Cp_salt = 2.39E-3, hA_radiatorNom = 0.0170, m_airN = 3.549878645002705, m_dot_air = 94.389 * 1.1237, m_dot_saltNom = 103.0701, m_saltN = 3.924401289158446e+02, tauUHXtoPHX = 8.24, tripRadiator = 2000000) annotation(
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
        SMD_MSR_Modelica.HeatTransport.DHRS dhrs( Cp_fluid = 1.9665E-3, DHRS_MaxPowerRm = 0, DHRS_PowerBleed = 0, DHRS_time = 10, DHRS_timeConstant = 10,m_DHRS = 1.625049507600000e+02, m_Pi_C_PHX = 1, m_dotNom = 7.5708E-02 * 2.14647E+03, tauDHRStoPHX = 1.385) annotation(
          Placement(visible = true, transformation(origin = {36, -2}, extent = {{-70, -70}, {70, 70}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac constantFlowFrac(ConstantFlowFrac = 1) annotation(
          Placement(visible = true, transformation(origin = {-80, -40}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0) annotation(
          Placement(visible = true, transformation(origin = {-87, 61}, extent = {{-25, -25}, {25, 25}}, rotation = 0)));
  SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantTempIn constantTempIn(TestTempIn = 100)  annotation(
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
          tauPHXtoUHX = 4.71)  annotation(
          Placement(visible = true, transformation(origin = {-16, -6}, extent = {{-52, -52}, {52, 52}}, rotation = 0)));
    
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac primaryFlowFrac(ConstantFlowFrac = 1)  annotation(
          Placement(visible = true, transformation(origin = {-76, 54}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
    
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac secondaryFlowFrac(ConstantFlowFrac = 1)  annotation(
          Placement(visible = true, transformation(origin = {52, -74}, extent = {{-28, -28}, {28, 28}}, rotation = 0)));
    
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
          Placement(visible = true, transformation(origin = {137, 29}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantHeatProduction constantHeatProduction(
          cP = 1.9665E-3, 
          mDot = 7.5708E-02*2.14647E+03, 
          nodeMass = 100, 
          qDotIn = 10)  annotation(
          Placement(visible = true, transformation(origin = {-140, 10}, extent = {{-34, -34}, {34, 34}}, rotation = 0)));
    
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantHeatRemoval constantHeatRemoval(
          cP = 2.39E-3,
          mDot = 5.36265E-02*1.922E3,
          nodeMass = 100,
          qDotOut = 10) annotation(
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
          Diagram(coordinateSystem(extent = {{-180, 80}, {180, -100}})));end TestPHX2;

      model TestFuelChannel2
        SMD_MSR_Modelica.Nuclear.FuelChannel fuelChannel(
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
          Placement(visible = true, transformation(origin = {-72, -12}, extent = {{-52, -52}, {52, 52}}, rotation = 180)));
      
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantHeatRemoval constantHeatRemoval(
          cP = 1.9665E-3, 
          mDot = 7.5708E-02*2.14647E+03, 
          nodeMass = 100, 
          qDotOut = 8) annotation(
          Placement(visible = true, transformation(origin = {-69, -87}, extent = {{-37, -37}, {37, 37}}, rotation = 0)));
      
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantFlowFrac fuelFlowFrac(ConstantFlowFrac = 1)  annotation(
          Placement(visible = true, transformation(origin = {32, 32}, extent = {{-36, -36}, {36, 36}}, rotation = 0)));
      
        SMD_MSR_Modelica.Nuclear.mPKE mPKE(
          Lam = 2.400E-04,
          lam = {1.240E-02,3.05E-02,1.11E-01,3.01E-01,1.140E+00,3.014E+00},
          beta = {0.000223,0.001457,0.001307,0.002628,0.000766,0.00023},
          tauCoreNom = 8.46,
          tauLoopNom = 16.73) annotation(
          Placement(visible = true, transformation(origin = {33, -55}, extent = {{-41, -41}, {41, 41}}, rotation = 0)));
      
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
          Placement(visible = true, transformation(origin = {-158, 4}, extent = {{-48, -48}, {48, 48}}, rotation = 180)));
        SMD_MSR_Modelica.QAtoolBox.TestIO.ConstantDecayPower constantDecayPower(DecayPower = 0)  annotation(
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
    end TestCases;
  end QAtoolBox;
  annotation(
    Documentation(info = "<html><head></head><body>SMD-MSR_ModelicaV1</body></html>"));
end SMD_MSR_Modelica;