%% Course element structs

%% c1
% geom
c1.L = 5e3;    %total course length [m]
c1.phi = [0 0];%   %angle of the slope of the course over the length of the course
% env
c1.rho = 1.1455; % density at location [kg/m^3]
c1.r_c = 1; %*** needs fixed to run without
c1.headwind = 0;
c1p = zeros(1,20);

%% c2
% geom
c2.L = 5e3;    %total course length [m]
c2.phi = [5 5];%   %angle of the slope of the course over the length of the course
% env
c2.rho = 1.1455; % density at location [kg/m^3]
c2.r_c = 1; %*** needs fixed to run without
c2.headwind = 0;
c2p = 5*ones(1,20);

%% c3
% geom
c3.L = 5e3;    %total course length [m]
c3.phi = [0 -5 0];
%-secd(linspace(-45,45,20)).^2;
%1./(1+linspace(-45,45,20).^2);%   %angle of the slope of the course over the length of the course
% env
c3.rho = 1.1455; % density at location [kg/m^3]
c3.r_c = 1; %*** needs fixed to run without
c3.headwind = 0;
c3p = [zeros(1,7) -5*ones(1,20) zeros(1,7)];
%% c4 here
% geom
c4.L = 5e3;    %total course length [m]
c4.phi = linspace(5,-5,20) ; % [linspace(5,-5,20) linspace(5,-5,20) linspace(5,-5,20) linspace(5,-5,20)];
%cosd()*360; %[45 -45];%   %angle of the slope of the course over the length of the course
% env
c4.rho = 1.1455; % density at location [kg/m^3]
c4.r_c = 1; %*** needs fixed to run without
c4.headwind = 0;
c4p = linspace(5,-5,20);

%% c5
% geom
c5.L = 5e3;    %total course length [m]
c5.phi = linspace(-5,5,20);
% env
c5.rho = 1.1455; % density at location [kg/m^3]
c5.r_c = 1; %*** needs fixed to run without
c5.headwind = 0;
c5p = linspace(-5,5,20);

%% c6
c6.L = 5e3;    %total course length [m]
c6.phi = [linspace(5,-5,20) linspace(-5,5,20)];
% env
c6.rho = 1.1455; % density at location [kg/m^3]
c6.r_c = 1; %*** needs fixed to run without
c6.headwind = 0;
%% c7
c7.L = 5e3;    %total course length [m]
c7.phi = [linspace(-5,5,20) linspace(5,-5,20)];
% env
c7.rho = 1.1455; % density at location [kg/m^3]
c7.r_c = 1; %*** needs fixed to run without
c7.headwind = 0;

%% c8
c8.L = 5e3;    %total course length [m]
c8.phi = [linspace(-5,5,20) zeros(1,20) linspace(-5,5,20)];
% env
c8.rho = 1.1455; % density at location [kg/m^3]
c8.r_c = 1; %*** needs fixed to run without
c8.headwind = 0;

%% c9
c9.L = 5e3;    %total course length [m]
c9.phi = [c4p c4p c1p c1p ];
% env
c9.rho = 1.1455; % density at location [kg/m^3]
c9.r_c = 1; %*** needs fixed to run without
c9.headwind = 0;

%% Combine
c = [c1;c2;c3;c4;c5;c6;c7;c8;c9];