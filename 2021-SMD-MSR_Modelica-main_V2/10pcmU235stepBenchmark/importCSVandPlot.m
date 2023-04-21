%% 10 pcm Step @ 3 power levels
% Pannel 1 -> 1MW
% Pannet 2 -> 5MW
% Pannel 3 -> 8MW

%Plot setup
zeroStamp = 5000;
start_plot = -10;
stop_plot = 350;
P = 8;

plot_width = stop_plot - start_plot; 

%%Load Results 1MW
omc1R = readmatrix('U235_1R_MSRE_1MW.csv');
omc9R = readmatrix('U235_9R_MSRE_1MW.csv');

timeOmc1R = omc1R(:,1);
timeOmc9R = omc9R(:,1);

timeOmc1RSSindex = find(timeOmc1R==zeroStamp,1);
timeOmc9RSSindex= find(timeOmc9R==zeroStamp,1);

powerOmc1RSS = omc1R(timeOmc1RSSindex,2);
powerOmc9RSS = omc9R(timeOmc9RSSindex,2);

figure(1)
subplot(3,1,1)
grid on
box on 
hold on
plot(timeOmc1R-zeroStamp,(omc1R(:,2)-powerOmc1RSS)*P,'color','#ff0000','LineWidth',2)
plot(timeOmc9R-zeroStamp,(omc9R(:,2)-powerOmc9RSS)*P,'--','color','#0000ff','LineWidth',2)
title('+10[pcm] step insertion @ 1[MW_t]MSRE-U235')
ylabel('\Delta Power [MW]')
legend('1R Modelica','9R Modelica')
xlim([start_plot plot_width]) 


%%Load Results 5MW
omc1R = readmatrix('U235_1R_MSRE_5MW.csv');
omc9R = readmatrix('U235_9R_MSRE_5MW.csv');

timeOmc1R = omc1R(:,1);
timeOmc9R = omc9R(:,1);

timeOmc1RSSindex = find(timeOmc1R==zeroStamp,1);
timeOmc9RSSindex= find(timeOmc9R==zeroStamp,1);

powerOmc1RSS = omc1R(timeOmc1RSSindex,2);
powerOmc9RSS = omc9R(timeOmc9RSSindex,2);

subplot(3,1,2)
grid on
box on 
hold on
plot(timeOmc1R-zeroStamp,(omc1R(:,2)-powerOmc1RSS)*P,'color','#ff0000','LineWidth',2)
plot(timeOmc9R-zeroStamp,(omc9R(:,2)-powerOmc9RSS)*P,'--','color','#0000ff','LineWidth',2)
title('+10[pcm] step insertion @ 5[MW_t] MSRE-U235')
ylabel('\Delta Power [MW]')
legend('1R Modelica','9R Modelica')
xlim([start_plot plot_width])  

%%Load Results
omc1R = readmatrix('U235_1R_MSRE_8MW.csv');
omc9R = readmatrix('U235_9R_MSRE_8MW.csv');

timeOmc1R = omc1R(:,1);
timeOmc9R = omc9R(:,1);

timeOmc1RSSindex = find(timeOmc1R==zeroStamp,1);
timeOmc9RSSindex= find(timeOmc9R==zeroStamp,1);

powerOmc1RSS = omc1R(timeOmc1RSSindex,2);
powerOmc9RSS = omc9R(timeOmc9RSSindex,2);

subplot(3,1,3)
grid on
box on 
hold on
plot(timeOmc1R-zeroStamp,(omc1R(:,2)-powerOmc1RSS)*P,'color','#ff0000','LineWidth',2)
plot(timeOmc9R-zeroStamp,(omc9R(:,2)-powerOmc9RSS)*P,'--','color','#0000ff','LineWidth',2)
title('+10[pcm] step insertion @ 8[MW_t] MSRE-U235')
ylabel('\Delta Power [MW]')
xlabel('Time [s]')
legend('1R Modelica','9R Modelica')
xlim([start_plot plot_width]) 

x0=10;
y0=10;
width=1100;
height=1050;
set(gcf,'position',[x0,y0,width,height])

% Save plot as fig and png
saveas(gcf,'U235Benchmark.png')
savefig('U235Benchmark.fig')