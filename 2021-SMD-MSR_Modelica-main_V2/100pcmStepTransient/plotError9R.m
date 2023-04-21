zeroStamp = 2000;
start_plot = -100;
stop_plot = 500;
plot_width = stop_plot - start_plot; 

%%Load Results
%omc = readmatrix('omc9RwithNoDH.csv');
omc = readmatrix('omc9RwithDH.csv');
load('tran1R9Slink.mat')

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

nomPowerSlink = interp1(xdataResamp,nomPowerSlinkSliced,timeOmc,"spline");

errorPowerNom = ((nomPowerOmc-nomPowerSlink)./nomPowerSlink);

errorPowerNomSm = smoothdata(errorPowerNom,'movmean',10000);

%Temp Results
tempInOmc = omc(timeOmcStartIndex:timeOmcEndIndex,3);
tempInSlinkSliced = tIn(timeSlinkStartIndex:timeSlinkEndIndex);

tempOutOmc = omc(timeOmcStartIndex:timeOmcEndIndex,4);
tempOutSlinkSliced = tOut(timeSlinkStartIndex:timeSlinkEndIndex);

tempInSlink = interp1(xdataResamp,tempInSlinkSliced,timeOmc,"spline");
tempOutSlink = interp1(xdataResamp,tempOutSlinkSliced,timeOmc,"spline");

errorTempIn = ((tempInOmc-tempInSlink)./tempInSlink);
errorTempOut = ((tempOutOmc-tempOutSlink)./tempOutSlink);

errorTempInSm = smoothdata(errorTempIn,'movmean',10000);
errorTempOutSm = smoothdata(errorTempOut,'movmean',10000);

%Feedback Results
totalFbOmc = omc(timeOmcStartIndex:timeOmcEndIndex,5);
totalFbSlinkSliced = rho_fb_tot(timeSlinkStartIndex:timeSlinkEndIndex);

totalFbSlink = interp1(xdataResamp,totalFbSlinkSliced,timeOmc,"spline");

errorFbTot = ((totalFbOmc-totalFbSlink)./totalFbSlink);

%errorFbTotSm = smoothdata(errorFbTot,'movmean',10000);
errorFbTotSm = smoothdata(errorFbTot,'movmean',100);

figure(1)
subplot(4,1,1)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink,nomPowerSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc,nomPowerOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Relative Power')
yyaxis right
plot(timeOmc,errorPowerNomSm,'color','#37783a','LineWidth',2);          
grid on
title('Normalized Total Power')
ylabel('Absolute Error')
%legend('9R Simulink','9R Modelica (No DH)')
legend('9R Simulink','9R Modelica (With DH)')
xlim([start_plot stop_plot])

subplot(4,1,2)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink,tempInSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc,tempInOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Temperature [^\circC]')
yyaxis right
plot(timeOmc,errorTempInSm,'color','#37783a','LineWidth',2);          
grid on
title('Core Fuel Inlet Temperature')
ylabel('Absolute Error')
%legend('9R Simulink','9R Modelica (No DH)')
legend('9R Simulink','9R Modelica (With DH)')
xlim([start_plot stop_plot])

subplot(4,1,3)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink,tempOutSlink,'color','#FF0000','LineWidth',2)
plot(timeOmc,tempOutOmc,'--','color','#0000ff','LineWidth',2)
ylabel('Temperature [^\circC]')
yyaxis right
plot(timeOmc,errorTempOutSm,'color','#37783a','LineWidth',2);
grid on
title('Core Fuel Outlet Temperature')
ylabel('Absolute Error')
% legend('9R Simulink','9R Modelica (No DH)')
legend('9R Simulink','9R Modelica (With DH)')
xlim([start_plot stop_plot])

subplot(4,1,4)
colororder({'b','#37783a'})
yyaxis left
box on 
hold on
plot(timeSlink,(totalFbSlink)*1E5,'color','#FF0000','LineWidth',2)
plot(timeOmc,(totalFbOmc)*1E5,'--','color','#0000ff','LineWidth',2)
ylabel('Feedback [pcm]')
yyaxis right
plot(timeOmc,errorFbTotSm,'color','#37783a','LineWidth',2);          
grid on
title('Total Temperature Feedback')
ylabel('Absolute Error')
%legend('9R Simulink','9R Modelica (No DH)')
legend('9R Simulink','9R Modelica (With DH)')
xlabel('Time [s]')
xlim([start_plot stop_plot])

x0=0;
y0=0;
width=1100;
height=1050;
set(gcf,'position',[x0,y0,width,height])

% Save plot as fig and png
% saveas(gcf,'tran19RnoDH.png')
% savefig('tran19RnoDH.fig')
saveas(gcf,'tran19RwithDH.png')
savefig('tran19RwithDH.fig')