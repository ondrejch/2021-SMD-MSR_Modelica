%% Read results csv
T = readtable('startUp_res.csv','VariableNamingRule','preserve');

%% Analysis
l = 0.017;
s = 1E6;
phi0 = 1.58E20;
t = T.time/3600; %turn into hours
n = T.("pke.n_population.n")*phi0;

t0Indx = find(t==1.9);
n0 = n(t0Indx);
M = n./n0;
k_eff1 = 1 - (1 ./M);

xStamp = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30];

[tf,indxt] = ismember(t,xStamp+0.9);
idx=1:length(t);
idx=idx(tf);
idx=idx(indxt(tf));

k_eff = 1 - ((l*s)./n);

%% Plot neutron population
figure(1)
box on 
grid on
hold on
plot(t,n,'color','#0000FF','LineWidth',1)
xline(13.5,'--','LineWidth',2)
xline(19.5,'--','LineWidth',2)
xline(24.5,'--','LineWidth',2)
ylabel('Neutron Population')
text(5,3.7E6,'Phase 1','FontSize',14,'color','#FF0000')
text(16,3.7E6,'Phase 2','FontSize',14,'color','#FF0000')
text(21,3.7E6,'Phase 3','FontSize',14,'color','#FF0000')
text(26,3.7E6,'Phase 4','FontSize',14,'color','#FF0000')
xticks(xStamp)
xlabel('Time [h]')
xlim([0 29]) 

x0=0;
y0=0;
width=2000;
height=1000;
set(gcf,'position',[x0,y0,width,height])

figure(2)
box on
grid on 
hold on
plot(t,n,'color','#0000FF','LineWidth',1)
ylabel('Neutron Population')
xticks(xStamp)
xlabel('Time [h]')
xlim([0 13.5]) 

x0=0;
y0=0;
width=2000;
height=1000;
set(gcf,'position',[x0,y0,width,height])

figure(3)
box on
grid on 
hold on
plot(t,n,'color','#0000FF','LineWidth',1)
ylabel('Neutron Population')
xticks(xStamp)
xlabel('Time [h]')
xlim([13.5 19.5]) 

x0=0;
y0=0;
width=2000;
height=1000;
set(gcf,'position',[x0,y0,width,height])

figure(4)
box on
grid on 
hold on
plot(t,n,'color','#0000FF','LineWidth',1)
ylabel('Neutron Population')
xticks(xStamp)
xlabel('Time [h]')
xlim([19.5 24.5]) 

x0=0;
y0=0;
width=2000;
height=1000;
set(gcf,'position',[x0,y0,width,height])

figure(5)
box on
grid on 
hold on
plot(t,n,'color','#0000FF','LineWidth',1)
ylabel('Neutron Population')
xticks(xStamp)
xlabel('Time [h]')
xlim([24.5 29])

x0=0;
y0=0;
width=2000;
height=1000;
set(gcf,'position',[x0,y0,width,height])

%% Plot external reactivity

figure(6)
box on 
grid on
hold on
plot(t,T.("pke.externalReactivityIn")*1E5,'color','#0000FF','LineWidth',1)
xline(13.5,'--','LineWidth',2)
xline(19.5,'--','LineWidth',2)
xline(24.5,'--','LineWidth',2)
ylabel('External Reactivity [pcm]')
text(5,300,'Phase 1','FontSize',14,'color','#FF0000')
text(16,300,'Phase 2','FontSize',14,'color','#FF0000')
text(21,300,'Phase 3','FontSize',14,'color','#FF0000')
text(26,300,'Phase 4','FontSize',14,'color','#FF0000')
xticks(xStamp)
xlabel('Time [h]')
xlim([0 29]) 

x0=0;
y0=0;
width=2000;
height=1000;
set(gcf,'position',[x0,y0,width,height])

%% Plot keff

figure(7)
box on 
grid on
hold on
plot(t,k_eff1,'color','#0000FF','LineWidth',1)
plot(t,k_eff)
xline(13.5,'--','LineWidth',2)
xline(19.5,'--','LineWidth',2)
xline(24.5,'--','LineWidth',2)
ylabel('Multipication Factor')
text(5,0.2,'Phase 1','FontSize',14,'color','#FF0000')
text(16,0.2,'Phase 2','FontSize',14,'color','#FF0000')
text(21,0.2,'Phase 3','FontSize',14,'color','#FF0000')
text(26,0.2,'Phase 4','FontSize',14,'color','#FF0000')
xticks(xStamp)
xlabel('Time [h]')
xlim([1 14]) 

x0=0;
y0=0;
width=2000;
height=1000;
set(gcf,'position',[x0,y0,width,height])