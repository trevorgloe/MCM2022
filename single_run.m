% run a single time
%% Define randome parameters
g = 9.81;   % acceleration due to gravity [m/s^2]

%% Define biker parameters
% m = 66.25;      % rider mass [kg]
m = 66.25;
Cr = 0.005;   % wheel resistance coefficient
% biker.CdA = 0.194*0.1;   % drag coeff of area [m^2]   (drag coefficient of 0.1 is also in there)
% CP = 142.8;   % rider critical power [Watts]
CP = 180;
Wcap = 24235;   % rider anaerobic work capacity [J]
tau_w = 500;   % W' recovery time constant
Pm = 225; % max power [Watts]

biker.m = m;
biker.Cr = Cr;
% biker.CdA = A;
biker.CP = CP;
biker.Wcap = Wcap;
biker.tau_w = tau_w;
biker.Pm = Pm;
biker.CdA = 0.194*0.1;   % drag coeff of area [m^2]   (drag coefficient of 0.1 is also in there)

%% Define course parameters
L = 13e3;    %total course length [m]
% phi = [1 2 3 4 6 7 8 7 6 5 4 3 2 1];   %angle of the slope of the course over the length of the course
% phi = [2 7 6 6 6 6 5 6 8 8 7 8 7.9];
phi = [0 1 0];
rho = 1.1455; % density at location [kg/m^3]
beta = linspace(1,2*pi,N);      % angle of current velocity and headwind
headwind = 5;   % [m/s]

%structure to store course parameters
course.L = L;
course.phi = phi;
course.rho = rho;
course.beta = beta;
course.headwind = headwind;

%% Discretize course into lil chunky bits

disc.N = 100;    %number of chunks in discretization
N = disc.N;

% sqp_solve;
[v,P,x] = sqp_run_new_wind(course, biker, disc);
convert_v;
plotting;