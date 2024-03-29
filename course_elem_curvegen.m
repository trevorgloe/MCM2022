clear all
close all
clc

%% Load Data
course_elem_structs;

%% Set Params
g = 9.81;   % acceleration due to gravity [m/s^2]

% Define biker parameters
biker.m = 66.25;      % rider mass [kg]
biker.Cr0 = 0.002;   % wheel resistance coefficient
biker.CdA = 0.172;   % drag coeff of area [m^2]   (drag coefficient of 0.1 is also in there)
biker.CP = 200;   % rider critical power [Watts]
biker.Wcap = -20e3;   % rider anaerobic work capacity [J]
biker.tau_w = 500;   % W' recovery time constant
biker.Pm = 300; % max power [Watts]

% Discretize course into lil chunky bits
disc.N = 150;    %number of chunks in discretization

%% Run Model
for qq = 3
    disp(qq)
    c(qq).beta = 0;
    [v,P,x] = course_elem_sqp_run_6(c(qq), biker, disc, qq);
    course_elem_convert_v;
    c(qq).time_values = Tf(end);
    all_params(:,qq) = {biker, c(qq), disc};
    
    course_elem_plotting;
end

