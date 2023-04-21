%%csv results order
%time
%Nom power
%Nom fission power
%Nom decay power
%Core inlet temp
%core outlet temp
%core graphite 
%total feedback
%fuel feedback
%grap feedback

clearvars -except tout Temp_mux rho_fb_tot rho_fb_f rho_fb_g

Results = readmatrix('Transient1UHX.csv');

timeSlink = tout/3600;
timeModelica = Results(:,1)/3600;

zeroStamp = 2000/3600;
start_plot = -500/3600;
stop_plot = 8000/3600;

plot_width = stop_plot - zeroStamp; 

%Figure 1 - Nominal total power and constituants
figure(1)
subplot(3,1,1)
grid on
box on 
hold on
plot(timeSlink-zeroStamp,(Temp_mux(:,1) + Temp_mux(:,2)),'color','#ff0000','LineWidth',1)
plot(timeModelica-zeroStamp,Results(:,2),'color','#0000ff','LineWidth',1)
title('Normalized Total Power')
ylabel('Relative Power')
legend('Simulink', 'Modelica')
xlim([start_plot plot_width]) 

subplot(3,1,2)
grid on
box on 
hold on
plot(timeSlink-zeroStamp,Temp_mux(:,1),'color','#ff0000','LineWidth',1)
plot(timeModelica-zeroStamp,Results(:,3),'color','#0000ff','LineWidth',1)
title('Normalized Fission Power')
ylabel('Relative Power')
legend('Simulink', 'Modelica')
xlim([start_plot plot_width]) 

subplot(3,1,3)
grid on
box on 
hold on
plot(timeSlink-zeroStamp,Temp_mux(:,2),'color','#ff0000','LineWidth',1)
plot(timeModelica-zeroStamp,Results(:,4),'color','#0000ff','LineWidth',1)
title('Normalized Decay Power')
ylabel('Relative Power')
legend('Simulink', 'Modelica')
xlim([start_plot plot_width]) 
xlabel('Time [h]')

x0=10;
y0=10;
width=1100;
height=1050;
set(gcf,'position',[x0,y0,width,height])

% Save plot as fig and png
saveas(gcf,'powerTransient1UHX.png')
savefig('powerTransient1UHX.fig')

%Figure 2 - Core inlet, outlet and graphite temp
figure(2)
subplot(3,1,1)
grid on
box on 
hold on
plot(timeSlink-zeroStamp,Temp_mux(:,6),'color','#ff0000','LineWidth',1)
plot(timeModelica-zeroStamp,Results(:,6),'color','#0000ff','LineWidth',1)
title('Core Outlet Temperatures')
legend('Simulink','Modelica')
ylabel('Temperature [^{\circ}C]')
xlim([start_plot plot_width]) 

subplot(3,1,2)
grid on
box on 
hold on
plot(timeSlink-zeroStamp,Temp_mux(:,3),'color','#ff0000','LineWidth',1)
plot(timeModelica-zeroStamp,Results(:,5),'color','#0000ff','LineWidth',1)
title('Core Inlet Temperatures')
legend('Simulink','Modelica')
ylabel('Temperature [^{\circ}C]')
xlim([start_plot plot_width]) 

subplot(3,1,3)
grid on
box on 
hold on
plot(timeSlink-zeroStamp,Temp_mux(:,4),'color','#ff0000','LineWidth',1)
plot(timeModelica-zeroStamp,Results(:,7),'color','#0000ff','LineWidth',1)
title('Graphite Temperatures')
legend('Simulink','Modelica')
ylabel('Temperature [^{\circ}C]')
xlim([start_plot plot_width]) 
xlabel('Time [h]')

x0=10;
y0=10;
width=1100;
height=1050;
set(gcf,'position',[x0,y0,width,height])

% Save plot as fig and png
saveas(gcf,'tempsTransient1UHX.png')
savefig('tempsTransient1UHX.fig')

%Figure 3 - Fuel and grap feedback
figure(3)
subplot(3,1,1)
grid on
box on 
hold on
plot(timeSlink-zeroStamp,rho_fb_tot*1E5,'color','#ff0000','LineWidth',1)
plot(timeModelica-zeroStamp,(Results(:,8)+Results(:,9)+Results(:,10))*1E5,'color','#0000ff','LineWidth',1)
title('Total Temperature Feedback')
legend('Simulink','Modelica')
ylabel('Reactivity [pcm]')
xlim([start_plot plot_width]) 

subplot(3,1,2)
grid on
box on 
hold on
plot(timeSlink-zeroStamp,rho_fb_f*1E5,'color','#ff0000','LineWidth',1)
plot(timeModelica-zeroStamp,(Results(:,8)+Results(:,9))*1E5,'color','#0000ff','LineWidth',1)
title('Fuel Temperature Feedback')
legend('Simulink','Modelica')
ylabel('Reactivity [pcm]')
xlim([start_plot plot_width]) 


subplot(3,1,3)
grid on
box on 
hold on
plot(timeSlink-zeroStamp,rho_fb_g*1E5,'color','#ff0000','LineWidth',1)
plot(timeModelica-zeroStamp,(Results(:,10))*1E5,'color','#0000ff','LineWidth',1)
title('Graphite Temperature Feedback')
legend('Simulink','Modelica')
ylabel('Reactivity [pcm]')
xlim([start_plot plot_width]) 
xlabel('Time [h]')

x0=10;
y0=10;
width=1100;
height=1050;
set(gcf,'position',[x0,y0,width,height])

% Save plot as fig and png
saveas(gcf,'reactTransient1UHX.png')
savefig('reactTransient1UHX.fig')

%%

% timeOffset = zeroStamp + start_plot;
% time_range = (0: 0.1 : 2400).';
% % 
% % tout = tout - timeOffset;
% % Results(:,1) = Results(:,1) - timeOffset;
% 
% timeOffsetIndexSlink = find(tout == 0);
% timeCutoffIndexSlink = find(tout == 2400);
% 
% timeOffsetIndexOMC = find(Results(:,1) == 0);
% timeCutoffIndexOMC = find(Results(:,1) == 2400);
% 
% timeSlink = tout(timeOffsetIndexSlink:timeCutoffIndexSlink);
% 
% fissPowerSlink = Temp_mux(timeOffsetIndexSlink:timeCutoffIndexSlink,1);
% decayPowerSlink = Temp_mux(timeOffsetIndexSlink:timeCutoffIndexSlink,2);
% nomPowerSlink = fissPowerSlink + decayPowerSlink;
% 
% inletTempSlink = Temp_mux(timeOffsetIndexSlink:timeCutoffIndexSlink,3);
% outletTempSlink = Temp_mux(timeOffsetIndexSlink:timeCutoffIndexSlink,6);
% grapTempSlink = Temp_mux(timeOffsetIndexSlink:timeCutoffIndexSlink,4);
% 
% totalFBslink = rho_fb_tot(timeOffsetIndexSlink:timeCutoffIndexSlink)*1E5;
% fuelFBslink = rho_fb_f(timeOffsetIndexSlink:timeCutoffIndexSlink)*1E5;
% grapFBslink = rho_fb_g(timeOffsetIndexSlink:timeCutoffIndexSlink)*1E5;
% 
% timeOMC = Results(timeOffsetIndexOMC:timeCutoffIndexOMC,1);
% nomPowerOMC = Results(timeOffsetIndexOMC:timeCutoffIndexOMC,2);
% fissPowerOMC = Results(timeOffsetIndexOMC:timeCutoffIndexOMC,3);
% decayPowerOMC = Results(timeOffsetIndexOMC:timeCutoffIndexOMC,4);
% 
% inletTempOMC = Results(timeOffsetIndexOMC:timeCutoffIndexOMC,5);
% outletTempOMC = Results(timeOffsetIndexOMC:timeCutoffIndexOMC,6);
% grapTempOMC = Results(timeOffsetIndexOMC:timeCutoffIndexOMC,7);
% 
% fuelFBOMC = (Results(timeOffsetIndexOMC:timeCutoffIndexOMC,8)+Results(timeOffsetIndexOMC:timeCutoffIndexOMC,9))*1E5;
% grapFBOMC = Results(timeOffsetIndexOMC:timeCutoffIndexOMC,10)*1E5;
% totalFBOMC = fuelFBOMC + grapFBOMC;
% 
% nomPowerSlinkIntpl = interp1(timeSlink,nomPowerSlink,time_range,'linear');
% fissPowerSlinkIntpl = interp1(timeSlink,fissPowerSlink,time_range,'linear');
% decayPowerSlinkIntpl = interp1(timeSlink,decayPowerSlink,time_range,'linear');
% 
% inletTempSlinkIntpl = interp1(timeSlink,inletTempSlink,time_range,'linear');
% outletTempSlinkIntpl = interp1(timeSlink,outletTempSlink,time_range,'linear');
% grapTempSlinkIntpl = interp1(timeSlink,grapTempSlink,time_range,'linear');
% 
% totalFBslinkIntpl = interp1(timeSlink,totalFBslink,time_range,'linear');
% fuelFBslinkIntpl = interp1(timeSlink,fuelFBslink,time_range,'linear');
% grapFBslinkIntpl = interp1(timeSlink,grapFBslink,time_range,'linear');
% 
% 
% % %Make OMC results unique
% % [timeOMCun, timeOMCindexUN] = unique(timeOMC); 
% % timeOMCintpl = interp1(timeOMC(timeOMCindexUN), timeOMC(timeOMCindexUN),time_range,'linear');
% % 
% % nomPowerOMCintpl = interp1(timeOMC(timeOMCindexUN), nomPowerOMC(timeOMCindexUN),time_range,'linear');
% % fissPowerOMCintpl = interp1(timeOMC(timeOMCindexUN),fissPowerOMC(timeOMCindexUN),time_range,'linear');
% % decayPowerOMCintpl = interp1(timeOMC(timeOMCindexUN),decayPowerOMC(timeOMCindexUN),time_range,'linear');
% % 
% % inletTempOMCintpl = interp1(timeOMC(timeOMCindexUN), inletTempOMC(timeOMCindexUN),time_range,'linear');
% % outletTempOMCintpl = interp1(timeOMC(timeOMCindexUN),outletTempOMC(timeOMCindexUN),time_range,'linear');
% % grapTempOMCintpl = interp1(timeOMC(timeOMCindexUN),grapTempOMC(timeOMCindexUN),time_range,'linear');
% % 
% % totalFBOMCintpl = interp1(timeOMC(timeOMCindexUN), totalFBOMC(timeOMCindexUN),time_range,'linear');
% % fuelFBOMCintpl = interp1(timeOMC(timeOMCindexUN), fuelFBOMC(timeOMCindexUN),time_range,'linear');
% % grapFBOMCintpl = interp1(timeOMC(timeOMCindexUN), grapFBOMC(timeOMCindexUN),time_range,'linear');
% 
% 
% 
% % nomPowerOMCintpl = interp1(timeOMC, nomPowerOMC,time_range,'linear');
% % fissPowerOMCintpl = interp1(timeOMCfissPowerOMC,time_range,'linear');
% % decayPowerOMCintpl = interp1(timeOMC,decayPowerOMC,time_range,'linear');
% % 
% % inletTempOMCintpl = interp1(timeOMC, inletTempOMC,time_range,'linear');
% % outletTempOMCintpl = interp1(timeOMC,outletTempOMC,time_range,'linear');
% % grapTempOMCintpl = interp1(timeOMC,grapTempOMC,time_range,'linear');
% % 
% % totalFBOMCintpl = interp1(timeOMC, totalFBOMC,time_range,'linear');
% % fuelFBOMCintpl = interp1(timeOMC, fuelFBOMC,time_range,'linear');
% % grapFBOMCintpl = interp1(timeOMC, grapFBOMC,time_range,'linear');
% 
% nomPowerEr = ((nomPowerOMCintpl-nomPowerSlinkIntpl)./nomPowerSlinkIntpl);
% fissPowerEr = ((fissPowerOMCintpl-fissPowerSlinkIntpl)./fissPowerSlinkIntpl);
% decayPowerEr = ((decayPowerOMCintpl-decayPowerSlinkIntpl)./decayPowerSlinkIntpl);
% 
% inletTempEr = ((inletTempOMCintpl-inletTempSlinkIntpl)./inletTempSlinkIntpl);
% outletTempEr = ((outletTempOMCintpl-outletTempSlinkIntpl)./outletTempSlinkIntpl);
% grapTempEr = ((grapTempOMCintpl-grapTempSlinkIntpl)./grapTempSlinkIntpl);
% 
% totalFbEr = ((totalFBOMCintpl-totalFBslinkIntpl)./totalFBslinkIntpl);
% fuelFbEr = ((fuelFBOMCintpl-fuelFBslinkIntpl)./fuelFBslinkIntpl);
% grapFbEr = ((grapFBOMCintpl-grapFBslinkIntpl)./grapFBslinkIntpl);
% 
% figure(4)
% subplot(2,1,1)
% grid on
% box on 
% hold on
% plot(time_range - zeroStamp,nomPowerSlinkIntpl,'color','#ff0000','LineWidth',1)
% plot(time_range - zeroStamp,nomPowerOMCintpl,'color','#0000ff','LineWidth',1)
% title('Normalized Total Power')
% ylabel('Relative Power')
% legend('Simulink', 'Modelica')
% xlim([start_plot plot_width]) 
% 
% subplot(2,1,2)
% grid on 
% box on
% hold on
% plot(time_range - zeroStamp,nomPowerEr,'color','#0000ff','LineWidth',1)
% title('Modelica Error based on Simulink')
% ylabel('Abs Error')
% xlim([start_plot plot_width]) 
% xlabel('Time [h]')
% 
% x0=10;
% y0=10;
% width=1100;
% height=1050;
% set(gcf,'position',[x0,y0,width,height])
% 
% % Save plot as fig and png
% saveas(gcf,'ERnomPower300pcm.png')
% savefig('ERnomPower300pcm.fig')
% 
% figure(5)
% subplot(2,1,1)
% grid on
% box on 
% hold on
% plot(time_range - zeroStamp,fissPowerSlinkIntpl,'color','#ff0000','LineWidth',1)
% plot(time_range - zeroStamp,fissPowerOMCintpl,'color','#0000ff','LineWidth',1)
% title('Normalized Fission Power')
% ylabel('Relative Power')
% legend('Simulink', 'Modelica')
% xlim([start_plot plot_width]) 
% 
% subplot(2,1,2)
% grid on 
% box on
% hold on
% plot(time_range - zeroStamp,fissPowerEr,'color','#0000ff','LineWidth',1)
% title('Modelica Error based on Simulink')
% ylabel('Abs Error')
% xlim([start_plot plot_width]) 
% xlabel('Time [h]')
% 
% x0=10;
% y0=10;
% width=1100;
% height=1050;
% set(gcf,'position',[x0,y0,width,height])
% 
% % Save plot as fig and png
% saveas(gcf,'ERfisPower300pcm.png')
% savefig('ERfisPower300pcm.fig')
% 
% figure(6)
% subplot(2,1,1)
% grid on
% box on 
% hold on
% plot(time_range - zeroStamp,decayPowerSlinkIntpl,'color','#ff0000','LineWidth',1)
% plot(time_range - zeroStamp,decayPowerOMCintpl,'color','#0000ff','LineWidth',1)
% title('Normalized Decay Power')
% ylabel('Relative Power')
% legend('Simulink', 'Modelica')
% xlim([start_plot plot_width]) 
% 
% subplot(2,1,2)
% grid on 
% box on
% hold on
% plot(time_range - zeroStamp,decayPowerEr,'color','#0000ff','LineWidth',1)
% title('Modelica Error based on Simulink')
% ylabel('Abs Error')
% xlim([start_plot plot_width]) 
% xlabel('Time [h]')
% 
% x0=10;
% y0=10;
% width=1100;
% height=1050;
% set(gcf,'position',[x0,y0,width,height])
% 
% % Save plot as fig and png
% saveas(gcf,'ERdecayPower300pcm.png')
% savefig('ERdecayPower300pcm.fig')
% 
% figure(7)
% subplot(2,1,1)
% grid on
% box on 
% hold on
% plot(time_range - zeroStamp,inletTempSlinkIntpl,'color','#ff0000','LineWidth',1)
% plot(time_range - zeroStamp,inletTempOMCintpl,'color','#0000ff','LineWidth',1)
% title('Core Inlet Temperature')
% ylabel('Temperature [^{\circ}C]')
% legend('Simulink', 'Modelica')
% xlim([start_plot plot_width]) 
% 
% subplot(2,1,2)
% grid on 
% box on
% hold on
% plot(time_range - zeroStamp,inletTempEr,'color','#0000ff','LineWidth',1)
% title('Modelica Error based on Simulink')
% ylabel('Abs Error')
% xlim([start_plot plot_width]) 
% xlabel('Time [h]')
% 
% x0=10;
% y0=10;
% width=1100;
% height=1050;
% set(gcf,'position',[x0,y0,width,height])
% 
% % Save plot as fig and png
% saveas(gcf,'ERtempIn300pcm.png')
% savefig('ERtempIn300pcm.fig')
% 
% figure(8)
% subplot(2,1,1)
% grid on
% box on 
% hold on
% plot(time_range - zeroStamp,outletTempSlinkIntpl,'color','#ff0000','LineWidth',1)
% plot(time_range - zeroStamp,outletTempOMCintpl,'color','#0000ff','LineWidth',1)
% title('Core Outlet Temperature')
% ylabel('Temperature [^{\circ}C]')
% legend('Simulink', 'Modelica')
% xlim([start_plot plot_width]) 
% 
% subplot(2,1,2)
% grid on 
% box on
% hold on
% plot(time_range - zeroStamp,outletTempEr,'color','#0000ff','LineWidth',1)
% title('Modelica Error based on Simulink')
% ylabel('Abs Error')
% xlim([start_plot plot_width]) 
% xlabel('Time [h]')
% 
% x0=10;
% y0=10;
% width=1100;
% height=1050;
% set(gcf,'position',[x0,y0,width,height])
% 
% % Save plot as fig and png
% saveas(gcf,'ERtempOut300pcm.png')
% savefig('ERtempOut300pcm.fig')
% 
% figure(9)
% subplot(2,1,1)
% grid on
% box on 
% hold on
% plot(time_range - zeroStamp,grapTempSlinkIntpl,'color','#ff0000','LineWidth',1)
% plot(time_range - zeroStamp,grapTempOMCintpl,'color','#0000ff','LineWidth',1)
% title('Graphite Temperature')
% ylabel('Temperature [^{\circ}C]')
% legend('Simulink', 'Modelica')
% xlim([start_plot plot_width]) 
% 
% subplot(2,1,2)
% grid on 
% box on
% hold on
% plot(time_range - zeroStamp,grapTempEr,'color','#0000ff','LineWidth',1)
% title('Modelica Error based on Simulink')
% ylabel('Abs Error')
% xlim([start_plot plot_width]) 
% xlabel('Time [h]')
% 
% x0=10;
% y0=10;
% width=1100;
% height=1050;
% set(gcf,'position',[x0,y0,width,height])
% 
% % Save plot as fig and png
% saveas(gcf,'ERtempGrap300pcm.png')
% savefig('ERtempGrap300pcm.fig')
% 
% figure(10)
% subplot(2,1,1)
% grid on
% box on 
% hold on
% plot(time_range - zeroStamp,totalFBslinkIntpl,'color','#ff0000','LineWidth',1)
% plot(time_range - zeroStamp,totalFBOMCintpl,'color','#0000ff','LineWidth',1)
% title('Total Temperature Feedback')
% ylabel('Reativity [pcm]')
% legend('Simulink', 'Modelica')
% xlim([start_plot plot_width]) 
% 
% subplot(2,1,2)
% grid on 
% box on
% hold on
% plot(time_range - zeroStamp,totalFbEr,'color','#0000ff','LineWidth',1)
% title('Modelica Error based on Simulink')
% ylabel('Abs Error')
% xlim([start_plot plot_width]) 
% xlabel('Time [h]')
% 
% x0=10;
% y0=10;
% width=1100;
% height=1050;
% set(gcf,'position',[x0,y0,width,height])
% 
% % Save plot as fig and png
% saveas(gcf,'ERtotalFB300pcm.png')
% savefig('ERtotalFBOut300pcm.fig')
% 
% figure(11)
% subplot(2,1,1)
% grid on
% box on 
% hold on
% plot(time_range - zeroStamp,fuelFBslinkIntpl,'color','#ff0000','LineWidth',1)
% plot(time_range - zeroStamp,fuelFBOMCintpl,'color','#0000ff','LineWidth',1)
% title('Fuel Temperature Feedback')
% ylabel('Reativity [pcm]')
% legend('Simulink', 'Modelica')
% xlim([start_plot plot_width]) 
% 
% subplot(2,1,2)
% grid on 
% box on
% hold on
% plot(time_range - zeroStamp,fuelFbEr,'color','#0000ff','LineWidth',1)
% title('Modelica Error based on Simulink')
% ylabel('Abs Error')
% xlim([start_plot plot_width]) 
% xlabel('Time [h]')
% 
% x0=10;
% y0=10;
% width=1100;
% height=1050;
% set(gcf,'position',[x0,y0,width,height])
% 
% % Save plot as fig and png
% saveas(gcf,'ERfuelFB300pcm.png')
% savefig('ERfuelFBOut300pcm.fig')
% 
% figure(12)
% subplot(2,1,1)
% grid on
% box on 
% hold on
% plot(time_range - zeroStamp,grapFBslinkIntpl,'color','#ff0000','LineWidth',1)
% plot(time_range - zeroStamp,grapFBOMCintpl,'color','#0000ff','LineWidth',1)
% title('Graphite Temperature Feedback')
% ylabel('Reativity [pcm]')
% legend('Simulink', 'Modelica')
% xlim([start_plot plot_width]) 
% 
% subplot(2,1,2)
% grid on 
% box on
% hold on
% plot(time_range - zeroStamp,grapFbEr,'color','#0000ff','LineWidth',1)
% title('Modelica Error based on Simulink')
% ylabel('Abs Error')
% xlim([start_plot plot_width]) 
% xlabel('Time [h]')
% 
% x0=10;
% y0=10;
% width=1100;
% height=1050;
% set(gcf,'position',[x0,y0,width,height])
% 
% % Save plot as fig and png
% saveas(gcf,'ERgrapFB300pcm.png')
% savefig('ERgrapFBOut300pcm.fig')
% 
% 
% 
% 
% 
% 
% 
