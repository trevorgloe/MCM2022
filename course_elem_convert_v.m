%% Converts and output of the nonlinear optimization, a velocity vector to a
% vector of power and calculates the total time
% must have velocity vector v defined
% version 2: includes implementation of curvature consideration
course = c(qq);

N = disc.N;
m = biker.m;
Cr = zeros(1,N);
% for qq = 1:disc.N
%     if course.r_c(qq) == inf %|| course.curvature_flag == 0
%         Cr(qq) = biker.Cr0;
%     else
%         Cr(qq) = biker.Cr0*sqrt(v(qq).^4 +g^2.*r_c(qq).^2)./(g.*r_c(qq));
%     end
% end
Cr0 = biker.Cr0;
Cr = Cr0*ones(1,N);
CP = biker.CP;
Wcap = biker.Wcap;
tau_w = biker.tau_w;
CdA = biker.CdA;
Pm = biker.Pm;

L = course.L;
phi = course.phi;
rho = course.rho;

dx = L/N;
v0 = ones(1,N);
x = linspace(0,L,N);
phi_dis = interp1(linspace(0,L,length(phi)),phi,x);
g = 9.8; % m/s

c1 = 0.5*rho*CdA;% area air drag coefficient
c2 = m.*g.*(sind(phi_dis) + Cr);
c3 = m; % get eq ***

dvdt = zeros(1,N);
dvdt(1) =  (v(1))/(dx*v(1));
for ii = 2:N
    dvdt(ii) = (v(ii) - v(ii-1))/(dx*v(ii));
end
% ineq constraint 1
Pcalc = (c1.*v.^2 + c2 + c3.*dvdt).*v; 

P_diff = Pcalc - P;
delta_v = delta(P,CP);
        
% check to see how much W' was actually used
x2 = repmat(x(1:30),1,ceil(length(x)/30));
shifted_x = x2(1:N);
        
Wptot = Wcap;
%         Wexp = 0;
Wptot = Wcap - sum((1-delta_v+delta_v.*exp(-shifted_x./(v*tau_w))).*(P-CP)*dx);

% total time
Tf = cumtrapz(x,1./v);
disp(['Total time: ' num2str(Tf(end)/60)])

function d = delta(P,CP)
    d = P < CP;
end
