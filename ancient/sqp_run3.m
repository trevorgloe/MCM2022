function [v,x] = sqp_run3(course, biker, disc,q)


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
    
    fun = @(v)modlag(v,c1,c2,c3,Wcap,tau_w,CP,L,dx,N,x,q);
    
    function [c,ceq] = constraint(v)
    % Backward difference approx for dvdt
        dvdt(1) =  (v(1))/(dx*v(1));
        for ii = 2:N
            dvdt(ii) = (v(ii) - v(ii-1))/(dx*v(ii));
        end
        % ineq constraint 1
        P = (c1.*v + c2 + c3.*dvdt).*v; 
        % P <= Pm;
        c = P - Pm*ones(1,N);
        
        ceq = [];
    end


    A = eye(N);
    b = ones(1,N)*30;
    Aeq = [];
    beq = [];
    lb = zeros(1,N);
    ub = [];%ones(1,N)*30;
    
    options = optimoptions('fmincon','Algorithm','sqp','Display','iter','MaxFunctionEvaluations',1e8,'StepTolerance',1e-10);
    v = fmincon(fun,v0,A,b,Aeq,beq,lb,ub,@constraint,options);
    
end

