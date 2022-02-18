%% Define parameters and discretize

%% Define randome parameters
g = 9.81;   % acceleration due to gravity [m/s^2]

%% Define biker parameters
m = 1;      % rider mass [kg]
Cr = 0.1;   % wheel resistance coefficient
A = 10;   % frontal area [TBD units]***
CP = 100;   % rider critical power
Wcap = 100;   % rider anaerobic work capacity
tau_w = 0.01;   % W' recovery time constant
Pm = 100; % max power [TBD units]***

biker.m = m;
biker.Cr = Cr;
biker.A = A;
biker.CP = CP;
biker.Wcap = Wcap;
biker.tau_w = tau_w;
biker.Pm = Pm;

%% Define course parameters
L = 100;    %total course length [m]
phi = [1 2 3 4 6 7 8 7 6 5 4 3 2 1];   %angle of the slope of the course over the length of the course
rho = 0.1; % density at location [TBD units]***

%structure to store course parameters
course.L = L;
course.phi = phi;
course.rho = rho;

%% Discretize course into lil chunky bits

disc.N = 100;    %number of chunks in discretization

