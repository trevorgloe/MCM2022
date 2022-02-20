% compare three bikers

clear all
close all

%% Set Params
g = 9.81;   % acceleration due to gravity [m/s^2]

% Define biker parameters
% biker.m = 66.25;      % rider mass [kg]
% biker.Cr0 = 0.002;   % wheel resistance coefficient
% biker.CdA = 0.172;   % drag coeff of area [m^2]   (drag coefficient of 0.1 is also in there)
% biker.CP = 180;   % rider critical power [Watts]
% biker.Wcap = 128e3;   % rider anaerobic work capacity [J]
% biker.tau_w = 500;   % W' recovery time constant
% biker.Pm = 300; % max power [Watts]

%% Biker 1 (TT specialist)
%male rider
% biker1.m = 71;      %[kg]
% biker1.Cr0 = 0.002;
% biker1.CdA = 0.35;
% biker1.CP = 357;
% biker1.Wcap = 22000;
% biker1.tau_w = 2000;
% biker1.Pm = 450;

%female rider
biker1.m = 61.6;      %[kg]
biker1.Cr0 = 0.002;
biker1.CdA = 0.318;
biker1.CP = 282;
biker1.Wcap = 17500;
biker1.tau_w = 2000;
biker1.Pm = 363;

%% Biker 2 (Climber)
% biker2.m = 62;
% biker2.Cr0 = 0.002;
% biker2.CdA = 0.33;
% biker2.CP = 320;
% biker2.Wcap = 45000;
% biker2.tau_w = 10;
% biker2.Pm = 398;

biker2.m = 54;
biker2.Cr0 = 0.002;
biker2.CdA = 0.30;
biker2.CP = 250;
biker2.Wcap = 31000;
biker2.tau_w = 10;
biker2.Pm = 310;

% %% Biker 3 (Sprinter)
% biker3.m = 76;
% biker3.Cr0 = 0.002;
% biker3.CdA = 0.37;
% biker3.CP = 290;
% biker3.Wcap = 35000;
% biker3.tau_w = 500;
% biker3.Pm = 800;


% Define course parameters
course.L = 1e3;    %total course length [m]

phi = [6 6 6 6];
% phi = [ 0 0 0 0];

course.rho = 1.1455; % density at location [kg/m^3]
% for ii = 1:length(R)
%     if mod(ii,10) == 0
%         r_c(ii/10) = R(ii);% curvature
%         course.beta(ii/10) = beta(ii);
%     end
% end
% course.r_c = r_c;
course.headwind = 0;
course.phi = phi;

disc.N = 100;    %number of chunks in discretization

% make all r_c values 1 (cause im ignoring them)
course.r_c = ones(1,disc.N);
course.beta = zeros(1,disc.N);

figure
plot(linspace(0,course.L,length(phi)),phi)
title('Course')
% Discretize course into lil chunky bits

N = disc.N;

%% Run Model
all_Ps = zeros(2,disc.N);
all_vs = zeros(2,disc.N);
all_Tf = zeros(1,2);

all_bikers = {biker1, biker2}

for k=1:2
    biker = all_bikers{k};
    
    % run it 3 times to make sure its actually optimal
    test_Ts = zeros(1,2);
    test_Ps = zeros(2,N);
    test_vs = zeros(2,N);
    for j = 1:2
        [v,P,x] = sqp_run_6(course, all_bikers{k}, disc);
        convert_v2_multiple;
        test_Ts(j) = Tf(end);
        all_params = {all_bikers{k}, course, disc};
        test_Ps(j,:) = P;
        test_vs(j,:) = v;
    end
    ind = find(test_Ts==min(test_Ts));
    all_Ps(k,:) = test_Ps(ind,:);
    all_vs(k,:) = test_vs(ind,:);
    all_Tf(k) = test_Ts(ind);
%     plotting
    save_data_params(all_params,v,P);
    
end

figure
hold on
for i=1:2
    plot(linspace(0,course.L,disc.N),all_Ps(i,:))
end

figure
hold on
for i=1:2
    plot(linspace(0,course.L,disc.N),all_vs(i,:))
end
