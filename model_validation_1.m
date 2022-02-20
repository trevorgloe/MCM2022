clear
close all

%% Load Data
distcalcs_verif

%% Set Params
g = 9.81;   % acceleration due to gravity [m/s^2]

% Define biker parameters
biker.m = 66.25;      % rider mass [kg]
biker.Cr0 = 0.002;   % wheel resistance coefficient
biker.CdA = 0.172;   % drag coeff of area [m^2]   (drag coefficient of 0.1 is also in there)
biker.CP = 180;   % rider critical power [Watts]
biker.Wcap = 128e3;   % rider anaerobic work capacity [J]
biker.tau_w = 500;   % W' recovery time constant
biker.Pm = 300; % max power [Watts]

% Define course parameters
course.L = 48.7e3;    %total course length [m]
course.phi = [0 1 0];%   %angle of the slope of the course over the length of the course
course.rho = 1.1455; % density at location [kg/m^3]
course.r_c = R;% curvature from ma boi callan***

% Discretize course into lil chunky bits
disc.N = 200;    %number of chunks in discretization

%% Run Model
[v,P,x] = sqp_run_5(course, biker, disc);
convert_v2;
time_values = Tf(end);
all_params = {biker, course, disc};

%% Plotting
plotting;

%% Model Comparison