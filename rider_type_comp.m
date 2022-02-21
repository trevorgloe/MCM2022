%% Rider Typer Comparisons
clear all
close all
clc

%% Load data
rider_structs;
get_courses;

%% Set Params
g = 9.81;   % acceleration due to gravity [m/s^2]

% Define course parameters
% course.L = 48.7e3;    %total course length [m]
% course.phi = phi;%   %angle of the slope of the course over the length of the course
% course.rho = 1.1455; % density at location [kg/m^3]
% for ii = 1:length(R)
%     if mod(ii,10) == 0
%         r_c(ii/10) = R(ii);% curvature
%         course.beta(ii/10) = beta(ii);
%     end
% end
% course.r_c = r_c;
course.headwind = 0;
course.L = 22e3;    %total course length [m]
% course.phi = [0 1 0];
bad_phi = 4*load('tokyo_inclination.mat').inclination_5sample;
bad_phi(find(isnan(bad_phi)))=0;
course.phi = lowpass(bad_phi,0.05)
course.rho = 1.1455; % density at location [kg/m^3]
course.beta = 0;

% Discretize course into lil chunky bits
disc.N = 100;
%length(course.r_c);    %number of chunks in discretization

%% Run Model
for ww = 1:4
    biker = rider_type(ww);
    biker.Cr0 = 0.002;   % wheel resistance coefficient
    
    [v2(ww,:),P2(ww,:),x2(ww,:)] = rider_type_sqp_run_6(course, biker, disc);
    %[v(ww,:),P(ww,:),x(ww,:)]
    rider_type_convert_v;
    time_values = Tf(end);
    all_params(:,ww) = {biker, course, disc};
end

%% Plotting

% plot male bikers
figure(1)
hold on
% for qq = 1:2
% %     x = linspace(0,all_params{2,qq}.L,all_params{3,qq}.N);
%     plot(x2(qq,:),v2(qq,:),'.-','MarkerSize', 10)
% end
x = linspace(0,all_params{2,1}.L,all_params{3,1}.N)
plot(x2(1,:),v2(1,:),'.-','MarkerSize',10)
plot(x2(3,:),v2(3,:),'.-','MarkerSize',10)
xlabel('x position')
ylabel('v')
grid on
hold off
legend('Time Trial - Male','Climber - Male')
% legend('Time Trial - Male','Time Trial - Female','Climber - Male','Climber - Female')

% plot female bikers
figure(2)
hold on
% for qq = 1:2
% %     x = linspace(0,all_params{2,qq}.L,all_params{3,qq}.N);
%     plot(x2(qq,:),v2(qq,:),'.-','MarkerSize', 10)
% end
% x = linspace(0,all_params{2,1}.L,all_params{3,1}.N)
plot(x2(2,:),v2(2,:),'.-','MarkerSize',10)
plot(x2(4,:),v2(4,:),'.-','MarkerSize',10)
xlabel('x position')
ylabel('v')
grid on
hold off
legend('Time Trial - Female','Climber - Female')

% plot males
figure(3)
hold on
% for qq = 1:4
% %     x = linspace(0,all_params{2,qq}.L,all_params{3,qq}.N);
%     plot(x2(qq,:),P2(qq,:),'.-','MarkerSize', 10)
% end
plot(x2(1,:),P2(1,:),'.-','MarkerSize',10)
plot(x2(3,:),P2(3,:),'.-','MarkerSize',10)
xlabel('x position')
ylabel('P')
grid on
hold off
legend('Time Trial - Male','Climber - Male')

% plot females
figure(4)
hold on
% for qq = 1:4
% %     x = linspace(0,all_params{2,qq}.L,all_params{3,qq}.N);
%     plot(x2(qq,:),P2(qq,:),'.-','MarkerSize', 10)
% end
plot(x2(2,:),P2(2,:),'.-','MarkerSize',10)
plot(x2(4,:),P2(4,:),'.-','MarkerSize',10)
xlabel('x position')
ylabel('P')
grid on
hold off
legend('Time Trial - Female','Climber - Female')