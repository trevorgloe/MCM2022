function [v,x] = sqp_run2(course, biker, disc,q)
% takes in rider and course struct and runs SQP based model
    
    %% Get params
    
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
    
    %% F function
%     fun = @(v) dx.*sum(1./v);
    function [Tp] = ModLag(v,q,c1,c2,c3,Wcap,tau_w,CP,L)
        
        % Backward difference approx for dvdt
        dvdt(1) =  (v(1))/(dx*v(1));
        for ii = 2:N
            dvdt(ii) = (v(ii) - v(ii-1))/(dx*v(ii));
        end
        
        P = (c1.*v + c2 + c3.*dvdt).*v;
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
        
        Lag = 1./v;
        % modified lagrangian
        mL = Lag + q*(G - Wcap/L);
        
        total_int = cumtrapz(x,mL);
        
        Tp = total_int(end);
    end

    % wrap function
    fun = @(v) ModLag(v,q,c1,c2,c3,Wcap,tau_w,CP,L);
        
    
    %% Constraint funct
    function [c,ceq] = constraint(v)%,x,Pm,N,Wcap,CP,c1,c2,c3,tau_w)
        % Backward difference approx for dvdt
        dvdt(1) =  (v(1))/(dx*v(1));
        for ii = 2:N
            dvdt(ii) = (v(ii) - v(ii-1))/(dx*v(ii));
        end
        % ineq constraint 1
        P = (c1.*v + c2 + c3.*dvdt).*v; 
        % P <= Pm;
        c = P - Pm*ones(1,N);
%         disp(max(P))
        if c > 0
            disp('c not satisfied')
        end
        if isnan(c)
%             disp('AHHHH')
%             disp(v)
        end
        ceq = [];
    end
    
    %% fmincon arguments
    A = eye(N);
    b = ones(1,N)*30;
    Aeq = [];
    beq = [];
    lb = zeros(1,N);
    ub = [];%ones(1,N)*30;
%     ub(1) = 0.001;
    
%     nonlcon =@ constraint; %,x,Pm,N,Wcap,CP,c1,c2,c3,tau_w);
    
    %% Call fmincon
    options = optimoptions('fmincon','Algorithm','sqp','Display','iter','MaxFunctionEvaluations',1e8,'StepTolerance',1e-10);
    v = fmincon(fun,v0,A,b,Aeq,beq,lb,ub,@constraint,options);
    disp('test')
    
end