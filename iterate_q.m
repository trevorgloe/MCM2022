%% iterate solving of modified lagrangian until contraint is satisified
q = 10e-4;     % initial guess
tol = 10;

err = 1e9;
prev_err = 0;
while abs(err) > tol
    params_discret
    sqp_solve
    convert_v
    
    % check to see if contraint is fullfilled
    Wprime = W(v,P,biker,course,disc);
    err = (Wcap - Wprime)
    
%     q = q+1e-4
%     diff_err = err - prev_err;
%     q = q + 1e-7*err
    
    err = 0
%     prev_err = err;
end


function [Wprime] = W(v,P,biker,course,disc)
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
    v0 = 20*linspace(0.1,1,N);
%     v0 = 10*ones(1,N);
    x = linspace(0,L,N);
    phi_dis = interp1(linspace(0,L,length(phi)),phi,x);
    g = 9.8; % m/s
    
    c1 = 0.5*rho*A;% air drag coefficient *** get equation
    c2 = m.*g.*(sind(phi_dis) + Cr);
    c3 = m; % get eq ***
    % get x values of changes
    for kk = 1:length(P)
        if P(kk) > CP
            x_flag(kk) = 1;
        else
            x_flag(kk) = 0;
        end
    end
    x_flag = P > CP;

    top_x = x;
    top_v = v;
    for kk = 2:length(P)
        if x_flag(kk) == 0
            ii = kk;
            poop = x_flag(ii);
            while poop == 0 && ii > 0
                poop = x_flag(ii);
                ind = (ii);
                ii = ii-1;
            end
            if ii == 0
                top_x(kk) = 0;
                top_v(kk) = v0(1);
            else
                top_x(kk) = x(ind);
                top_v(kk) = v(ind);
            end
        end  
    end        

    % integral of modified lagrangian with lagrangian multiplier q
    delta = P < CP;
    bigphi = (P-CP).*(1 + delta.*exp(((top_x./top_v)-(x./v))/tau_w));

    G = bigphi./v;
    
    all_vals = cumtrapz(x,G);
    Wprime = all_vals(end);
end

