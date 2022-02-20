clear
close all

%% Load Data

%% Set Params
g = 9.81;   % acceleration due to gravity [m/s^2]

% Define biker parameters
biker.m = 66.25;      % rider mass [kg]
biker.Cr = 0.005;   % wheel resistance coefficient
biker.CdA = 0.194*0.1;   % drag coeff of area [m^2]   (drag coefficient of 0.1 is also in there)
biker.CP = 180;   % rider critical power [Watts]
biker.Wcap = 24235;   % rider anaerobic work capacity [J]
biker.tau_w = 500;   % W' recovery time constant
biker.Pm = 350; % max power [Watts]

% Define course parameters
course.L = 13e3;    %total course length [m]
course.phi = [1 2 3 4 6 7 8 7 6 5 4 3 2 1];   %angle of the slope of the course over the length of the course
course.rho = 1.1455; % density at location [kg/m^3]

% Discretize course into lil chunky bits
disc.N = 100;    %number of chunks in discretization

%% Run Model
[v,P,x] = sqp_run_new(course, biker, disc);
convert_v;
time_values = Tf(end);
all_params = {biker, course, disc};

%% Plotting
plotting;

%% Model Comparison