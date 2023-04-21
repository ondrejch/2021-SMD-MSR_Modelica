zeroStamp = 5000;
start_plot = -10;
stop_plot = 350;
P = 8;

plot_width = stop_plot - start_plot; 


%%Load Results 1MW
Results = readmatrix('U235step1MW.csv');
load('U235step1MW.mat')

timeSlink = tout;
timeModelica = Results(:,1);

timeSlinkSSindex = find(timeSlink==zeroStamp);
timeModSSindex = find(timeModelica==zeroStamp,1);

powerModSS = Results(timeModSSindex,2);

powerSlink = Temp_mux(:,1) + Temp_mux(:,2);
powerSlinkSS = Temp_mux(timeSlinkSSindex,1) + Temp_mux(timeSlinkSSindex,2);

figure(1)
subplot(3,1,1)
grid on
box on 
hold on
plot(timeSlink-zeroStamp,(powerSlink-powerSlinkSS)*P,'color','#ff0000','LineWidth',2)
plot(timeModelica-zeroStamp,(Results(:,2)-powerModSS)*P,'--','color','#0000ff','LineWidth',2)
title('+10[pcm] step insertion @ 1MW MSRE-U233')
ylabel('\Delta Power [MW]')
legend('Simulink', 'Modelica')
xlim([start_plot plot_width]) 

%%Load Results 5MW
Results = readmatrix('U235step5MW.csv');
load('U235step5MW.mat')

timeSlink = tout;
timeModelica = Results(:,1);

timeSlinkSSindex = find(timeSlink==zeroStamp);
timeModSSindex = find(timeModelica==zeroStamp,1);

powerModSS = Results(timeModSSindex,2);

powerSlink = Temp_mux(:,1) + Temp_mux(:,2);
powerSlinkSS = Temp_mux(timeSlinkSSindex,1) + Temp_mux(timeSlinkSSindex,2);

subplot(3,1,2)
grid on
box on 
hold on
plot(timeSlink-zeroStamp,(powerSlink-powerSlinkSS)*P,'color','#ff0000','LineWidth',2)
plot(timeModelica-zeroStamp,(Results(:,2)-powerModSS)*P,'--','color','#0000ff','LineWidth',2)
title('+10[pcm] step insertion @ 5MW MSRE-U233')
ylabel('\Delta Power [MW]')
legend('Simulink', 'Modelica')
xlim([start_plot plot_width])  

%%Load Results
Results = readmatrix('U235step8MW.csv');
load('U235step8MW.mat')

timeSlink = tout;
timeModelica = Results(:,1);

timeSlinkSSindex = find(timeSlink==zeroStamp);
timeModSSindex = find(timeModelica==zeroStamp,1);

powerModSS = Results(timeModSSindex,2);

powerSlink = Temp_mux(:,1) + Temp_mux(:,2);
powerSlinkSS = Temp_mux(timeSlinkSSindex,1) + Temp_mux(timeSlinkSSindex,2);

subplot(3,1,3)
grid on
box on 
hold on
plot(timeSlink-zeroStamp,(powerSlink-powerSlinkSS)*P,'color','#ff0000','LineWidth',2)
plot(timeModelica-zeroStamp,(Results(:,2)-powerModSS)*P,'--','color','#0000ff','LineWidth',2)
title('+10[pcm] step insertion @ 8MW MSRE-U235')
ylabel('\Delta Power [MW]')
xlabel('[s]')
legend('Simulink', 'Modelica')
xlim([start_plot plot_width]) 

x0=10;
y0=10;
width=1100;
height=1050;
set(gcf,'position',[x0,y0,width,height])

% Save plot as fig and png
saveas(gcf,'U235stepBenchmark.png')
savefig('U235stepBenchmark.fig')