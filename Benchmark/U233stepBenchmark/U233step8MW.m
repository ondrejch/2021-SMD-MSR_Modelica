%%% Scalable Modular Dynamic MSR Model
%%% Authors - Visura Pathirana, Alex Wheeler
%%% Building on work done by Vikram Sinha and Alex Wheeler
%%% Project advisor - Dr. Ondrej Chvala

%% Benchmark Transient - 1 --> +10[pcm] for U235 MSRE @ 1MW, 5MW & 8MW
%%% Runs step insertions of +10[pcm] for U235 MSRE @ 1MW, 5MW & 8MW
%%% Results are ploted along with MSRE one region, nine region in three panels for each power
%%% Simulations done in four steps in a loop

%%% Step - 1; Simulation is run for 2000[s] at 8[Mw_t]
%%% Step - 2; Using UHX_MODE = 2, simulation is brought to 8[MW_t]
%%% Step - 3; Allow 3000[s] to come to steady state at 8[MW_t]
%%% Step - 4; +10.0[pcm] worth reactivity is inserted as a step at 5000[s]

%% Benchmark 1 - Inputs

%%% Basic Simulation Parameters
simtime = 5500;                                                            %Simulation time [s]
ts_max = 1e-2;                                                             %Maximum timestep [s] 
P=8;                                                                       %Operational thermal power [MW]

%%% Fuel Type
%%% fuel_type = 235; for FLibe with U235
%%% fuel_type = 233; for FLiBe with U233
fuel_type = 233;                                                           

%%% Source Step Reactivity Insertions & Sinusoidal Reactvity Insertions 
sourcedata = [0 0 0];                                                      %Neutron source insertions [abs]
sourcetime = [0 1000 2500];                                                %Neutron source insertion time [s]
source = timeseries(sourcedata,sourcetime);                                %Defining source timeseries  

reactdata = [0 10.0E-5 10.0E-5];                                           %Reactivity insertions [abs]
reacttime = [0 5000 6000];                                                 %Reactivity insertion time [s]
react = timeseries(reactdata,reacttime);                                   %Defining source timeseries

omega          = 10.00000;                                                 %Frequncy of the sine wave [rad]
sin_mag        = 0;                                                        %Amplitude of the sine wave [abs]
dx             = round((2*pi/omega)/25, 2, 'significant');                 %Size of the time step [s]

%%% Pump Trips
Trip_P_pump=20000000;                                                      %Time at which primary pump is tripped [s]
Trip_S_pump=20000000;                                                      %Time at which secondary pump is tripped [s]

%%% UHX Parameters
%%% UHX_MODE = 1; uses a radiator
%%% UHX_MODE = 2; uses a constant power removal block
%%% Demand following is allowed only in UHX_MODE = 2 (constant power removal block)
UHX_MODE = 2;
Block_UHE=20000000;                                                        %Time at which ultimate heat exchanger will be cut off [s]

DHRS_Power = 0;
Power_Bleed = 0;

%%% Only for UHX_MODE = 2
demanddata = [1 1 1 1 1];                                            %Reactivity insertions [abs]
demandtime = [0 1000 2000 3000 6000];                                      %Reactivity insertion time [s]
demand = timeseries(demanddata,demandtime);                                %Defining source timeseries

%%% DHRS Parameters
%%% DHRS_MODE = 1; a sigmoid based DHRS (Normal DHRS)
%%% DHRS_MODE = 2; a square pulse based DHRS (Broken DHRS)
%%% DHRS_MODE = 1 allows modifications to sigmoid behavior. Listed under Normal DHRS below
%%% DHRS_MODE = 2 allows modifications through matlab function block from sim
%%% DHRS_MODE = 2 allows cold slug insertions
DHRS_MODE = 1; 
DHRS_time=20000000;                                                        %Time at which DRACS will be activated. use simtime to keep it off [s]

%%% Only for DHRS_MODE = 2
deltaTf_DHRS = 30;                                                         %Temperature drop by broken DHRS [deg. C]
slug_time = 100000;                                                        %Duration of slug [s]

run('SMD_MSR_Para_V1')
sim('SMD_MSR_Sim_V1.slx');

save('U233step8MW.mat','tout','Temp_mux','rho_fb_tot','rho_fb_f','rho_fb_g')