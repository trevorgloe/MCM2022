%% Define parameters and discretize

%% Define randome parameters
g = 9.81;   % acceleration due to gravity [m/s^2]

%% Define biker parameters
m = 1;      % rider mass [kg]
Cr = 0.1;   % wheel resistance coefficient
Cd = 0.1;   % air drag coefficient
CP = 100;   % rider critical power
Wp = 100;   % rider anaerobic work capacity

biker.m = m;
biker.Cr = Cr;
biker.Cd = Cd;
biker.CP = CP;
biker.Wp = Wp;


%% Define course parameters
L = 100;    %total course length [m]
phi = [1 2 3 4 6 7 8 7 6 5 4 3 2 1];   %angle of the slope of the course over the length of the course

%structure to store course parameters
course.L = L;
course.phi = phi;

%% Discretize course into lil chunky bits
N = 100;    %number of chunks in discretization

x_vals = linspace(0,L,N);
phi_dis = interp1(linspace(0,L,length(phi)),phi,x_vals);

dx = x_vals(end)-x_vals(end-1);

c2_array = m*g*(phi_dis + Cr);

