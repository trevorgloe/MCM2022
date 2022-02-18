%% Define parameters and discretize

%% Define randome parameters
g = 9.81;   % acceleration due to gravity [m/s^2]

%% Define biker parameters
m = 1;      % rider mass [kg]
Cr = 0.1;   % wheel resistance coefficient
Cd = 0.1;   % air drag coefficient
CP = 100;   % rider critical power
Wcap = 100;   % rider anaerobic work capacity
tau_w = 0.01;   % W' recovery time constant

biker.m = m;
biker.Cr = Cr;
biker.Cd = Cd;
biker.CP = CP;
biker.Wp = Wp;
biker.tau_w = tau_w;

%% Define course parameters
L = 100;    %total course length [m]
phi = [1 2 3 4 6 7 8 7 6 5 4 3 2 1];   %angle of the slope of the course over the length of the course

%structure to store course parameters
course.L = L;
course.phi = phi;

%% Discretize course into lil chunky bits

disc.N = 100;    %number of chunks in discretization

disc.x = linspace(0,L,N);
disc.phi_dis = interp1(linspace(0,L,length(phi)),phi,x);

disc.c2_array = m*g*(phi_dis + Cr);

