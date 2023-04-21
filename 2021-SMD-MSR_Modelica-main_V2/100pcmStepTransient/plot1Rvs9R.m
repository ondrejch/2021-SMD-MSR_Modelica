zeroStamp = 2000;
start_plot = -100;
stop_plot = 300;
plot_width = stop_plot - start_plot; 

%%Load Results
omc1R = readmatrix('omc1R.csv');
omc9R = readmatrix('omc9RwithDH.csv');

%%Find SS time index
timeOmc1RSSindex = find(omc1R(:,1)==zeroStamp);
timeOmc9RSSindex = find(omc9R(:,1)==zeroStamp);

%%Find SS power
powerOmc1RSS = omc1R(timeOmc1RSSindex,2);
powerOmc9RSS = omc9R(timeOmc9RSSindex,2);

%%Zero out time vectors
timeOmc1ROg = omc1R(:,1)-zeroStamp;
timeOmc9ROg = omc9R(:,1)-zeroStamp;

%Indexing start and end of vectors
timeOmc1RStartIndex = find(timeOmc1ROg==start_plot);
timeOmc1REndIndex = find(timeOmc1ROg==stop_plot);

timeOmc9RStartIndex = find(timeOmc9ROg==start_plot);
timeOmc9REndIndex = find(timeOmc9ROg==stop_plot);

%Slice time vectors
timeOmc1R = timeOmc1ROg(timeOmc1RStartIndex:timeOmc1REndIndex);
timeOmc9R = timeOmc9ROg(timeOmc9RStartIndex:timeOmc9REndIndex);

%Power Results
nomPowerOmc1R = omc1R(timeOmc1RStartIndex:timeOmc1REndIndex,2);
nomPowerOmc9R = omc9R(timeOmc9RStartIndex:timeOmc9REndIndex,2);

%Temp Results
tempInOmc1R = omc1R(timeOmc1RStartIndex:timeOmc1REndIndex,5);
tempInOmc9R = omc9R(timeOmc9RStartIndex:timeOmc9REndIndex,3);

tempOutOmc1R = omc1R(timeOmc1RStartIndex:timeOmc1REndIndex,6);
tempOutOmc9R = omc9R(timeOmc9RStartIndex:timeOmc9REndIndex,4);

%Feedback Results
totalFbOmc1R = omc1R(timeOmc1RStartIndex:timeOmc1REndIndex,11);
totalFbOmc9R = omc9R(timeOmc9RStartIndex:timeOmc9REndIndex,5);

figure(1)
subplot(4,1,1)
box on 
hold on
plot(timeOmc1R,nomPowerOmc1R,'color','#FF0000','LineWidth',2)
plot(timeOmc9R,nomPowerOmc9R,'--','color','#0000ff','LineWidth',2)
ylabel('Relative Power')        
grid on
title('Normalized Total Power')
legend('1R Modelica','9R Modelica')
xlim([start_plot stop_plot])

subplot(4,1,2)
box on 
hold on
plot(timeOmc1R,tempInOmc1R,'color','#FF0000','LineWidth',2)
plot(timeOmc9R,tempInOmc9R,'--','color','#0000ff','LineWidth',2)
ylabel('Temperature [^\circC]')        
grid on
title('Core Fuel Inlet Temperature')
legend('1R Modelica','9R Modelica')
xlim([start_plot stop_plot])

subplot(4,1,3)
box on 
hold on
plot(timeOmc1R,tempOutOmc1R,'color','#FF0000','LineWidth',2)
plot(timeOmc9R,tempOutOmc9R,'--','color','#0000ff','LineWidth',2)
ylabel('Temperature [^\circC]')
grid on
title('Core Fuel Outlet Temperature')
legend('1R Modelica','9R Modelica')
xlim([start_plot stop_plot])

subplot(4,1,4)
box on 
hold on
plot(timeOmc1R,(totalFbOmc1R)*1E5,'color','#FF0000','LineWidth',2)
plot(timeOmc9R,(totalFbOmc9R)*1E5,'--','color','#0000ff','LineWidth',2)
ylabel('Feedback [pcm]') 
xlabel('Time [s]')
grid on
title('Total Temperature Feedback')
legend('1R Modelica','9R Modelica')
xlim([start_plot stop_plot])

x0=0;
y0=0;
width=1100;
height=1050;
set(gcf,'position',[x0,y0,width,height])

% Save plot as fig and png
saveas(gcf,'tran1_1Rvs9R.png')
savefig('tran1_1Rvs9R.fig')