%% Course element structs

%% c1
% geom
c1.L = 5e3;    %total course length [m]
c1.phi = [0 0];%   %angle of the slope of the course over the length of the course
% env
c1.rho = 1.1455; % density at location [kg/m^3]
c1.r_c = 1; %*** needs fixed to run without
c1.headwind = 0;

%% c2
% geom
c2.L = 5e3;    %total course length [m]
c2.phi = [5 5];%   %angle of the slope of the course over the length of the course
% env
c2.rho = 1.1455; % density at location [kg/m^3]
c2.r_c = 1; %*** needs fixed to run without
c2.headwind = 0;

%% c3
% geom
c3.L = 5e3;    %total course length [m]
c3.phi = [-5 -5];
%-secd(linspace(-45,45,20)).^2;
%1./(1+linspace(-45,45,20).^2);%   %angle of the slope of the course over the length of the course
% env
c3.rho = 1.1455; % density at location [kg/m^3]
c3.r_c = 1; %*** needs fixed to run without
c3.headwind = 0;

%% c4
% geom
c4.L = 5e3;    %total course length [m]
c4.phi = linspace(5,-5,20);
%cosd()*360; %[45 -45];%   %angle of the slope of the course over the length of the course
% env
c4.rho = 1.1455; % density at location [kg/m^3]
c4.r_c = 1; %*** needs fixed to run without
c4.headwind = 0;

%% c5
% geom
c5.L = 5e3;    %total course length [m]
c5.phi = linspace(-5,5,20);
%[-cosd(linspace(-40,40,20))*360 cosd(linspace(-40,40,20))*360];%[-45 45];%   %angle of the slope of the course over the length of the course
% env
c5.rho = 1.1455; % density at location [kg/m^3]
c5.r_c = 1; %*** needs fixed to run without
c5.headwind = 0;

%% Combine
c = [c1;c2;c3;c4;c5];