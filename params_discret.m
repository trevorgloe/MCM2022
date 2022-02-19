%% Define parameters and discretize

%% Define randome parameters
g = 9.81;   % acceleration due to gravity [m/s^2]

%% Define biker parameters
m = 66.25;      % rider mass [kg]
Cr = 0.001;   % wheel resistance coefficient
A = 0.194;   % frontal area [m^2]
CP = 142.8;   % rider critical power [Watts]
Wcap = 24235;   % rider anaerobic work capacity [J]
tau_w = 580;   % W' recovery time constant
Pm = 350; % max power [Watts]

biker.m = m;
biker.Cr = Cr;
biker.A = A;
biker.CP = CP;
biker.Wcap = Wcap;
biker.tau_w = tau_w;
biker.Pm = Pm;

%% Define course parameters
L = 44.2e3;    %total course length [m]
% phi = [1 2 3 4 6 7 8 7 6 5 4 3 2 1];   %angle of the slope of the course over the length of the course
% phi = [0 0 0 0 0 0 0];
phi = [0 0 -1 -2 -3 -5 -8 -9 -7 -4 -3 -2 1 3 7 9];
% phi = totaldist.m;
rho = 1.1455; % density at location [kg/m^3]

% course.data = 'TokyoWomens_indivTT.mat';
% get_phi

%structure to store course parameters
course.L = L;
course.phi = phi;
course.rho = rho;



%% Discretize course into lil chunky bits

disc.N = 120;    %number of chunks in discretization

