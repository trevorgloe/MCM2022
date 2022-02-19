% converts and output of the nonlinear optimization, a velocity vector to a
% vector of power and calculates the total time

% must have velocity vector v defined

m = biker.m;
Cr = biker.Cr;
%     Cd = biker.Cd;
CP = biker.CP;
Wcap = biker.Wcap;
tau_w = biker.tau_w;
A = biker.A;
Pm = biker.Pm;

L = course.L;
phi = course.phi;
rho = course.rho;

N = disc.N;

dx = L/N;
v0 = ones(1,N);
x = linspace(0,L,N);
phi_dis = interp1(linspace(0,L,length(phi)),phi,x);
g = 9.8; % m/s

c1 = 0.5*rho*A;% air drag coefficient *** get equation
c2 = m.*g.*(phi_dis + Cr);
c3 = m; % get eq ***

% power
P = (c1*v + c2.*v + c3*[diff(v) 0]).*v;

% total time
Tf = cumtrapz(x,1./v)
