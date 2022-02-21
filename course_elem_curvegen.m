clear all
clc

%% Load Data
course_elem_structs;

%% Set Params
g = 9.81;   % acceleration due to gravity [m/s^2]

% Define biker parameters
rider_structs;
biker = c_m;
biker.Cr0 = 0.002;   % wheel resistance coefficient
% biker.m = 66.25;      % rider mass [kg]
% biker.CdA = 0.172;   % drag coeff of area [m^2]   (drag coefficient of 0.1 is also in there)
% biker.CP = 200;   % rider critical power [Watts]
% biker.Wcap = -20e3;   % rider anaerobic work capacity [J]
% biker.tau_w = 500;   % W' recovery time constant
% biker.Pm = 300; % max power [Watts]


% Discretize course into lil chunky bits
disc.N = 100;    %number of chunks in discretization

%% Run Model
for qq = 1
    disp(qq)
    
    biker = tt_m;
    biker.Cr0 = 0.002;
    c(qq).beta = [0*ones(1,25) 90*ones(1,25) 180*ones(1,25) 270*ones(1,25)];
    [vtt(qq,:),Ptt(qq,:),xtt(qq,:)] = course_elem_sqp_run_6(c(qq), biker, disc, qq);
    Ptt(find(Ptt<0)) = 0;
    v = vtt(qq,:);
    P = Ptt(qq,:);
    x = xtt(qq,:);
    course_elem_convert_v;
    c(qq).time_values = Tf(end);
    all_params(:,qq) = {biker, c(qq), disc};
    
    biker = c_m;
    biker.Cr0 = 0.002;  
    
    [vc(qq,:),Pc(qq,:),xc(qq,:)] = course_elem_sqp_run_6(c(qq), biker, disc, qq);
    Pc(find(Pc<0)) = 0;
    v4 = vc(qq,:);
    P4 = Pc(qq,:);
    x4 = xc(qq,:);
    course_elem_convert_v;
    c(qq).time_values = Tf(end);
    all_params(:,qq) = {biker, c(qq), disc};
end



