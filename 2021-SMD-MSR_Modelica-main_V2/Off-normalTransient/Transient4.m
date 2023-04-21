%%% Scalable Modular Dynamic MSR Model
%%% Authors - Visura Pathirana, Alex Wheeler
%%% Building on work done by Vikram Sinha and Alex Wheeler
%%% Project advisor - Dr. Ondrej Chvala

%% Transient - 2
%%% Runs a simulation of UHX failure with DHRS open & SCRAM simultaniously
%%% Results are plotted in three panels; Reactor power, Core tempratures & feedbacks
%%% Simulation done in five steps

%%% Step - 1; Simulation is run for 2000[s] at 8[Mw_t]
%%% Step - 2; Using UHX_MODE = 2, UHX tripped at 2000[s]
%%% Step - 3; Using DHRS_MODE = 1, 8% Normal DHRS open at 2000[s]
%%% Step - 4; -1400[pcm] worth reactivity inserted as a step at 2000[s]
%%% Step - 5; Simulation continued till 12000[s]

%% User Inputs

%%% Basic Simulation Parameters
simtime = 8000;                                                            %Simulation time [s]
ts_max = 1e-2;                                                             %Maximum timestep [s] 
P=8;                                                                       %Operational thermal power [MW]

%%% Fuel Type
%%% fuel_type = 235; for FLibe with U235
%%% fuel_type = 233; for FLiBe with U233
fuel_type = 235;                                                           

%%% Source Step Reactivity Insertions & Sinusoidal Reactvity Insertions 
sourcedata = [0 0 0];                                                      %Neutron source insertions [abs]
sourcetime = [0 1000 2500];                                                %Neutron source insertion time [s]
source = timeseries(sourcedata,sourcetime);                                %Defining source timeseries  

reactdata = [0 0 -1400E-5 -1400E-5 -1400E-5 -1400E-5];                        %Reactivity insertions [abs]
reacttime = [0 2000 2000 5000 7500 15000];                                 %Reactivity insertion time [s]
react = timeseries(reactdata,reacttime);                                   %Defining source timeseries

omega          = 10.00000;                                                 %Frequncy of the sine wave [rad]
sin_mag        = 0;                                                        %Amplitude of the sine wave [abs]
dx             = round((2*pi/omega)/25, 2, 'significant');                 %Size of the time step [s]

%%% Pump Trips
Trip_P_pump=2000;                                                      %Time at which primary pump is tripped [s]
Trip_S_pump=2000;                                                      %Time at which secondary pump is tripped [s]

%%% UHX Parameters
%%% UHX_MODE = 1; uses a radiator
%%% UHX_MODE = 2; uses a constant power removal block
%%% Demand following is allowed only in UHX_MODE = 2 (constant power removal block)
UHX_MODE = 2;
Block_UHE=2000;                                                            %Time at which ultimate heat exchanger will be cut off [s]

%%% Only for UHX_MODE = 2
demanddata = [1 1 1 1 1];                                                  %Reactivity insertions [abs]
demandtime = [0 1000 2000 3000 5000];                                      %Reactivity insertion time [s]
demand = timeseries(demanddata,demandtime);                                %Defining source timeseries                               %Defining source timeseries

%%% DHRS Parameters
%%% DHRS_MODE = 1; a sigmoid based DHRS (Normal DHRS)
%%% DHRS_MODE = 2; a square pulse based DHRS (Broken DHRS)
%%% DHRS_MODE = 1 all  ows modifications to sigmoid behavior. Listed under Normal DHRS below
%%% DHRS_MODE = 2 allows modifications through matlab function block from sim
%%% DHRS_MODE = 2 allows cold slug insertions
DHRS_MODE = 1; 
DHRS_time=2000;                                                            %Time at which DRACS will be activated. use simtime to keep it off [s]

%%% Only for DHRS_MODE = 1
DHRS_Power= P*(0.08);                                                      %Maximum power that can be removed by DHRS
Power_Bleed= P*(0.00);                                                     %Some power will removed from DRACS even when its not used 

%%% Only for DHRS_MODE = 2
deltaTf_DHRS = 30;                                                         %Temperature drop by broken DHRS [deg. C]
slug_time = 8.46;                                                          %Duration of slug [s]


run('SMD_MSR_Para_V1')
sim('SMD_MSR_Sim_V1.slx');

powN = Temp_mux(:,1) + Temp_mux(:,2);
powF = Temp_mux(:,1);
powD = Temp_mux(:,2);
tIn = Temp_mux(:,3);
tGrap = Temp_mux(:,4);
tOut = Temp_mux(:,6);

save('tran4R1Slink.mat','tout','powN','powF','powD','tIn','tGrap','tOut','rho_fb_tot','rho_fb_f','rho_fb_g')