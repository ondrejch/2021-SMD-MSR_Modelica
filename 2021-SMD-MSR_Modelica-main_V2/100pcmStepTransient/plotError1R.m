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
start_plot = -100;
stop_plot = 300;
plot_width = stop_plot - start_plot; 

%%Load Results
omc = readmatrix('omc1R.csv');
load('tran1R1Slink.mat')

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

errorPowerNomSm = smoothdata(errorPowerNom,'movmean',10000);
errorPowerFisSm = smoothdata(errorPowerFis,'movmean',10000);
errorPowerDecSm = smoothdata(errorPowerDec,'movmean',10000);

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

errorFbFuelSm = smoothdata(errorFbFuel,'movmean',1000);
errorFbGrapSm = smoothdata(errorFbGrap,'movmean',10000);
errorFbTotSm = smoothdata(errorFbTot,'movmean',1000);

figure(1)
subplot(3,1,1)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink,nomPowerSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc,nomPowerOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Relative Power')
yyaxis right
hold on
plot(timeOmc,errorPowerNomSm,'color','#37783a','LineWidth',2);
grid on
title('Normalized Total Power')
ylabel('Absolute Error')
legend('Simulink','Modelica')
xlim([start_plot stop_plot])

subplot(3,1,2)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink,fisPowerSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc,fisPowerOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Relative Power')
yyaxis right
hold on
plot(timeOmc,errorPowerFisSm,'color','#37783a','LineWidth',2);
grid on
title('Normalized Fission Power')
ylabel('Absolute Error')
legend('Simulink','Modelica')
xlim([start_plot stop_plot])

subplot(3,1,3)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink,decPowerSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc,decPowerOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Relative Power')
yyaxis right
hold on
plot(timeOmc,errorPowerDecSm,'color','#37783a','LineWidth',2);
grid on
title('Normalized Decay Power')
ylabel('Absolute Error')
xlabel('Time [s]')
legend('Simulink','Modelica')
xlim([start_plot stop_plot])

x0=0;
y0=0;
width=1100;
height=1050;
set(gcf,'position',[x0,y0,width,height])

saveas(gcf,'tran1Power1R.png')
savefig('tran1Power1R.fig')

figure(2)
subplot(3,1,1)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink,tempInSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc,tempInOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Temperature [^\circC]')
yyaxis right
hold on
plot(timeOmc,errorTempInSm,'color','#37783a','LineWidth',2) 
grid on
title('Core Fuel Inlet Temperature')
ylabel('Absolute Error')
legend('Simulink','Modelica')
xlim([start_plot stop_plot])

subplot(3,1,2)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink,tempOutSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc,tempOutOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Temperature [^\circC]')
yyaxis right
hold on
plot(timeOmc,errorTempOutSm,'color','#37783a','LineWidth',2)
grid on
title('Core Fuel Outlet Temperature')
ylabel('Absolute Error')
legend('Simulink','Modelica')
xlim([start_plot stop_plot])

subplot(3,1,3)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink,tempGrapSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc,tempGrapOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Temperature [^\circC]')
yyaxis right
hold on
plot(timeOmc,errorTempGrapSm,'color','#37783a','LineWidth',2)
grid on
title('Core Graphite Temperature')
ylabel('Absolute Error')
xlabel('Time [s]')
legend('Simulink','Modelica')
xlim([start_plot stop_plot])

x0=0;
y0=0;
width=1100;
height=1050;
set(gcf,'position',[x0,y0,width,height])

saveas(gcf,'tran1Temp1R.png')
savefig('tran1Temp1R.fig')

figure(3)
subplot(3,1,1)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink,(totalFbSlink)*1E5,'color','#FF0000','LineWidth',2)
plot(timeOmc,(totalFbOmc)*1E5,'--','color','#0000ff','LineWidth',2)
ylabel('Feedback [pcm]')
yyaxis right
hold on
plot(timeOmc,errorFbTotSm,'color','#37783a','LineWidth',2)
grid on
title('Total Temperature Feedback')
ylabel('Absolute Error')
legend('Simulink','Modelica')
xlim([start_plot stop_plot])

subplot(3,1,2)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink,(fuelFbSlink)*1E5,'color','#FF0000','LineWidth',2)
plot(timeOmc,(fuelFbOmc)*1E5,'--','color','#0000ff','LineWidth',2)
ylabel('Feedback [pcm]')
yyaxis right
hold on
plot(timeOmc,errorFbFuelSm,'color','#37783a','LineWidth',2)
grid on
title('Fuel Temperature Feedback')
ylabel('Absolute Error')
legend('Simulink','Modelica')
xlim([start_plot stop_plot])

subplot(3,1,3)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink,(grapFbSlink)*1E5,'color','#FF0000','LineWidth',2)
plot(timeOmc,(grapFbOmc)*1E5,'--','color','#0000ff','LineWidth',2)
ylabel('Feedback [pcm]')
yyaxis right
hold on
plot(timeOmc,errorFbGrapSm,'color','#37783a','LineWidth',2)
grid on
title('Graphite Temperature Feedback')
ylabel('Absolute Error')
xlabel('Time [s]')
legend('Simulink','Modelica')
xlim([start_plot stop_plot])

x0=0;
y0=0;
width=1100;
height=1050;
set(gcf,'position',[x0,y0,width,height])

saveas(gcf,'tran1React1R.png')
savefig('tran1React1R.fig')