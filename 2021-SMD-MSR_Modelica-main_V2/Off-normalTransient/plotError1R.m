%time 1
%Nom power 2
%Nom fission power 3
%Nom decay power 4
%Core inlet temp 5
%core outlet temp 6
%core graphite 7
%fuel feedback 8
%grap feedback 9 
%total feedback 10

zeroStamp = 2000;
start_plot = -500;
stop_plot = 5000;
plot_width = stop_plot - start_plot; 

%%Load Results
omc = readmatrix('omc1R.csv');
load('tran4R1Slink.mat')

%%Find SS time index
timeOmcSSindex= find(omc(:,1)==zeroStamp);
timeSlinkSSindex = find(tout==zeroStamp);

%%Find SS power
powerSlinkSS = powN(timeSlinkSSindex);
powerOmcRSS = omc(timeOmcSSindex,2);

%%Zero out time vectors
timeOmcOg = omc(:,1)-zeroStamp;
timeSlinkOg = tout-zeroStamp;

%Indexing start and end of vectors
timeOmcStartIndex = find(timeOmcOg==start_plot);
timeOmcEndIndex = find(timeOmcOg==stop_plot);

timeSlinkStartIndex = find(timeSlinkOg==start_plot);
timeSlinkEndIndex = find(timeSlinkOg==stop_plot);

%Slice time vectors
timeSlinkSliced = timeSlinkOg(timeSlinkStartIndex:timeSlinkEndIndex);
timeOmc = timeOmcOg(timeOmcStartIndex:timeOmcEndIndex);

%Interpolate SLink time
xdataResamp = linspace(timeOmc(1),timeOmc(end),numel(timeSlinkSliced))'; 
timeSlink = interp1(xdataResamp,timeSlinkSliced,timeOmc);

%Power Results
nomPowerOmc = omc(timeOmcStartIndex:timeOmcEndIndex,2);
nomPowerSlinkSliced = powN(timeSlinkStartIndex:timeSlinkEndIndex);

fisPowerOmc = omc(timeOmcStartIndex:timeOmcEndIndex,3);
fisPowerSlinkSliced = powF(timeSlinkStartIndex:timeSlinkEndIndex);

decPowerOmc = omc(timeOmcStartIndex:timeOmcEndIndex,4);
decPowerSlinkSliced = powD(timeSlinkStartIndex:timeSlinkEndIndex);

nomPowerSlink = interp1(xdataResamp,nomPowerSlinkSliced,timeOmc,"spline");
fisPowerSlink = interp1(xdataResamp,fisPowerSlinkSliced,timeOmc,"spline");
decPowerSlink = interp1(xdataResamp,decPowerSlinkSliced,timeOmc,"spline");

errorPowerNom = ((nomPowerOmc-nomPowerSlink)./nomPowerSlink);
errorPowerFis = ((fisPowerOmc-fisPowerSlink)./fisPowerSlink);
errorPowerDec = ((decPowerOmc-decPowerSlink)./decPowerSlink);

errorPowerNomSm = smoothdata(errorPowerNom,'movmean',1000);
errorPowerFisSm = smoothdata(errorPowerFis,'movmean',1000);
errorPowerDecSm = smoothdata(errorPowerDec,'movmean',1000);

%Temp Results
tempInOmc = omc(timeOmcStartIndex:timeOmcEndIndex,5);
tempInSlinkSliced = tIn(timeSlinkStartIndex:timeSlinkEndIndex);

tempOutOmc = omc(timeOmcStartIndex:timeOmcEndIndex,6);
tempOutSlinkSliced = tOut(timeSlinkStartIndex:timeSlinkEndIndex);

tempGrapOmc = omc(timeOmcStartIndex:timeOmcEndIndex,7);
tempGrapSlinkSliced = tGrap(timeSlinkStartIndex:timeSlinkEndIndex);

tempInSlink = interp1(xdataResamp,tempInSlinkSliced,timeOmc,"spline");
tempOutSlink = interp1(xdataResamp,tempOutSlinkSliced,timeOmc,"spline");
tempGrapSlink = interp1(xdataResamp,tempGrapSlinkSliced,timeOmc,"spline");

errorTempIn = ((tempInOmc-tempInSlink)./tempInSlink);
errorTempOut = ((tempOutOmc-tempOutSlink)./tempOutSlink);
errorTempGrap = ((tempGrapOmc-tempGrapSlink)./tempGrapSlink);

errorTempInSm = smoothdata(errorTempIn,'movmean',10000);
errorTempOutSm = smoothdata(errorTempOut,'movmean',10000);
errorTempGrapSm = smoothdata(errorTempGrap,'movmean',10000);

%Feedback Results
fuelFbOmc = omc(timeOmcStartIndex:timeOmcEndIndex,8)+omc(timeOmcStartIndex:timeOmcEndIndex,9);
fuelFbSlinkSliced = rho_fb_f(timeSlinkStartIndex:timeSlinkEndIndex);

grapFbOmc = omc(timeOmcStartIndex:timeOmcEndIndex,10);
grapFbSlinkSliced = rho_fb_g(timeSlinkStartIndex:timeSlinkEndIndex);

totalFbOmc = omc(timeOmcStartIndex:timeOmcEndIndex,11);
totalFbSlinkSliced = rho_fb_tot(timeSlinkStartIndex:timeSlinkEndIndex);

fuelFbSlink = interp1(xdataResamp,fuelFbSlinkSliced,timeOmc,"spline");
grapFbSlink = interp1(xdataResamp,grapFbSlinkSliced,timeOmc,"spline");
totalFbSlink = interp1(xdataResamp,totalFbSlinkSliced,timeOmc,"spline");

errorFbFuel = ((fuelFbOmc-fuelFbSlink)./fuelFbSlink);
errorFbGrap = ((grapFbOmc-grapFbSlink)./grapFbSlink);
errorFbTot = ((totalFbOmc-totalFbSlink)./totalFbSlink);

errorFbFuelSm = smoothdata(errorFbFuel,'movmean',100000);
errorFbGrapSm = smoothdata(errorFbGrap,'movmean',10000);
errorFbTotSm = smoothdata(errorFbTot,'movmean',100000);

figure(1)
subplot(3,1,1)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink/3600,nomPowerSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc/3600,nomPowerOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Relative Power')
yyaxis right
hold on
plot(timeOmc/3600,errorPowerNomSm,'color','#37783a','LineWidth',2);
grid on
title('Normalized Total Power')
ylabel('Absolute Error')
legend('Simulink','Modelica')
xlim([start_plot/3600 stop_plot/3600])

subplot(3,1,2)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink/3600,fisPowerSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc/3600,fisPowerOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Relative Power')
yyaxis right
hold on
plot(timeOmc/3600,errorPowerFisSm,'color','#37783a','LineWidth',2);
grid on
title('Normalized Fission Power')
ylabel('Absolute Error')
legend('Simulink','Modelica')
xlim([start_plot/3600 stop_plot/3600])

subplot(3,1,3)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink/3600,decPowerSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc/3600,decPowerOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Relative Power')
yyaxis right
hold on
plot(timeOmc/3600,errorPowerDecSm,'color','#37783a','LineWidth',2);
grid on
title('Normalized Decay Power')
ylabel('Absolute Error')
xlabel('Time [h]')
legend('Simulink','Modelica')
xlim([start_plot/3600 stop_plot/3600])

x0=0;
y0=0;
width=1100;
height=1050;
set(gcf,'position',[x0,y0,width,height])

saveas(gcf,'tran4Power1R.png')
savefig('tran4Power1R.fig')

figure(2)
subplot(3,1,1)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink/3600,tempInSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc/3600,tempInOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Temperature [^\circC]')
yyaxis right
hold on
plot(timeOmc/3600,errorTempInSm,'color','#37783a','LineWidth',2) 
grid on
title('Core Fuel Inlet Temperature')
ylabel('Absolute Error')
legend('Simulink','Modelica')
xlim([start_plot/3600 stop_plot/3600])

subplot(3,1,2)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink/3600,tempOutSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc/3600,tempOutOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Temperature [^\circC]')
yyaxis right
hold on
plot(timeOmc/3600,errorTempOutSm,'color','#37783a','LineWidth',2)
grid on
title('Core Fuel Outlet Temperature')
ylabel('Absolute Error')
legend('Simulink','Modelica')
xlim([start_plot/3600 stop_plot/3600])

subplot(3,1,3)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink/3600,tempGrapSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc/3600,tempGrapOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Temperature [^\circC]')
yyaxis right
hold on
plot(timeOmc/3600,errorTempGrapSm,'color','#37783a','LineWidth',2)
grid on
title('Core Graphite Temperature')
ylabel('Absolute Error')
xlabel('Time [h]')
legend('Simulink','Modelica')
xlim([start_plot/3600 stop_plot/3600])

x0=0;
y0=0;
width=1100;
height=1050;
set(gcf,'position',[x0,y0,width,height])

saveas(gcf,'tran4Temp1R.png')
savefig('tran4Temp1R.fig')

figure(3)
subplot(3,1,1)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink/3600,(totalFbSlink)*1E5,'color','#FF0000','LineWidth',2)
plot(timeOmc/3600,(totalFbOmc)*1E5,'--','color','#0000ff','LineWidth',2)
ylabel('Feedback [pcm]')
yyaxis right
hold on
plot(timeOmc/3600,errorFbTotSm,'color','#37783a','LineWidth',2)
grid on
title('Total Temperature Feedback')
ylabel('Absolute Error')
legend('Simulink','Modelica')
xlim([start_plot/3600 stop_plot/3600])

subplot(3,1,2)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink/3600,(fuelFbSlink)*1E5,'color','#FF0000','LineWidth',2)
plot(timeOmc/3600,(fuelFbOmc)*1E5,'--','color','#0000ff','LineWidth',2)
ylabel('Feedback [pcm]')
yyaxis right
hold on
plot(timeOmc/3600,errorFbFuelSm,'color','#37783a','LineWidth',2)
grid on
title('Fuel Temperature Feedback')
ylabel('Absolute Error')
legend('Simulink','Modelica')
xlim([start_plot/3600 stop_plot/3600])

subplot(3,1,3)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink/3600,(grapFbSlink)*1E5,'color','#FF0000','LineWidth',2)
plot(timeOmc/3600,(grapFbOmc)*1E5,'--','color','#0000ff','LineWidth',2)
ylabel('Feedback [pcm]')
yyaxis right
hold on
plot(timeOmc/3600,errorFbGrapSm,'color','#37783a','LineWidth',2)
grid on
title('Graphite Temperature Feedback')
ylabel('Absolute Error')
xlabel('Time [h]')
legend('Simulink','Modelica')
xlim([start_plot/3600 stop_plot/3600])

x0=0;
y0=0;
width=1100;
height=1050;
set(gcf,'position',[x0,y0,width,height])

saveas(gcf,'tran4React1R.png')
savefig('tran4React1R.fig')