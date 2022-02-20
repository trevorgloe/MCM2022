% run modeled for sampled parameters

clear all
close all
%% Define randome parameters
g = 9.81;   % acceleration due to gravity [m/s^2]

%% Define biker parameters
% m = 66.25;      % rider mass [kg]
m = 66.25;
Cr = 0.005;   % wheel resistance coefficient
% biker.CdA = 0.194*0.1;   % drag coeff of area [m^2]   (drag coefficient of 0.1 is also in there)
% CP = 142.8;   % rider critical power [Watts]
% CP = 180;
% Wcap = 24235;   % rider anaerobic work capacity [J]
% tau_w = 500;   % W' recovery time constant
% Pm = 350; % max power [Watts]

biker.m = m;
biker.Cr = Cr;
% biker.CdA = A;
% biker.CP = CP;
% biker.Wcap = Wcap;
% biker.tau_w = tau_w;
% biker.Pm = Pm;
biker.CdA = 0.194*0.1;   % drag coeff of area [m^2]   (drag coefficient of 0.1 is also in there)

%% Define course parameters
L = 13e3;    %total course length [m]
% phi = [1 2 3 4 6 7 8 7 6 5 4 3 2 1];   %angle of the slope of the course over the length of the course
% phi = [2 7 6 6 6 6 5 6 8 8 7 8 7.9];
phi = [0 1 0];
rho = 1.1455; % density at location [kg/m^3]

%structure to store course parameters
course.L = L;
course.phi = phi;
course.rho = rho;

%% Discretize course into lil chunky bits

disc.N = 100;    %number of chunks in discretization
N = disc.N;

sample_biker_params
% sample parameters are stored in param_mesh


for point = 1:length(param_mesh)
    biker.CP = param_mesh(point).CP;
    biker.Pm = param_mesh(point).Pm;
    biker.Wcap = param_mesh(point).Wcap;
    biker.tau_w = param_mesh(point).tau_w
    
    
    sqp_solve;
    convert_v;
    
    data_run(point).time = Tf(end);
    data_run(point).P = P;
    data_run(point).v = v;
    
    data_run(point).CP = param_mesh(point).CP;
    data_run(point).Pm = param_mesh(point).Pm;
    data_run(point).Wcap = param_mesh(point).Wcap;
    data_run(point).tau_w = param_mesh(point).tau_w;
    
    all_params = {biker, course, disc};
    save_data_params(all_params,v,P);
    
end
    
    
    