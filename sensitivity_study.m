%% Define parameters and discretize for sensitivity study

%% Define randome parameters
g = 9.81;   % acceleration due to gravity [m/s^2]

%% Define biker parameters
% m = 66.25;      % rider mass [kg]
m = 66.25;
Cr = 0.005;   % wheel resistance coefficient
A = 0.194*0.1;   % frontal area [m^2]   (drag coefficient of 0.1 is also in there)
% CP = 142.8;   % rider critical power [Watts]
% CP = 180;
Wcap = 24235;   % rider anaerobic work capacity [J]
tau_w = 500;   % W' recovery time constant
Pm = 350; % max power [Watts]

% biker.m = m;
biker.Cr = Cr;
biker.A = A;
biker.CP = CP;
biker.Wcap = Wcap;
biker.tau_w = tau_w;
biker.Pm = Pm;

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

CP_vec = linspace(95,220,9);
for i=1:length(CP_vec)
    % run the model 3 times with random initial conditions and take the
    % best run
    CP = CP_vec(i);
    biker.CP = CP;
    time_values = zeros(1,3);
    Ptot = zeros(3,N);
    vtot = zeros(3,N);
    for j=1:3
        sqp_solve;
    	convert_v;
        time_values(j)=Tf(end);
        Ptot(j,:)=P;
        vtot(j,:)=v;
    end
    best = find(time_values==min(time_values));
    best_P = Ptot(best,:);
    best_v = vtot(best,:);
    all_params = {biker, course, disc};
    save_data_params(all_params,best_v,best_P);
end
    

