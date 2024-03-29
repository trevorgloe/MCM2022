function [v,P,x] = sqp_run_new(course, biker, disc)
% takes in rider and course struct and runs SQP based model
    
    %% Get params
    
    m = biker.m;
    Cr = biker.Cr;
    CdA = biker.CdA;
    CP = biker.CP;
    Wcap = biker.Wcap;
    tau_w = biker.tau_w;
    Pm = biker.Pm;
    
    L = course.L;
    phi = course.phi;
    rho = course.rho;
    
    N = disc.N;
    
    dx = L/N;
    v0 = 20*rand(1,N);
    P0 = 150*rand(1,N);
    s0 = [v0 P0];
    x = linspace(0,L,N);
    phi_dis = interp1(linspace(0,L,length(phi)),phi,x);
    g = 9.8; % m/s
    
    c1 = 0.5*rho*CdA;% air drag coefficient *** get equation
    c2 = m.*g.*(sind(phi_dis) + Cr);
    
    % assume wheel is flat disk with mass mw and radius rw
    mw = 1.4;       %[kg]
    rw = 0.559/2;   %[m]
    Iw = 1/2*mw*rw^2;
    c3 = m + 2*Iw/rw^2; % effective mass of the rider and bike
    
    %% F function
    fun = @(s) dx.*sum(1./s(1:N));
    
    %% Constraint funct
    function [c,ceq] = constraint(s)%,x,Pm,N,Wcap,CP,c1,c2,c3,tau_w)
        v = s(1:N);
        
        % Backward difference approx for dvdt
        dvdt(1) =  (v(1))/(dx*v(1));
        for ii = 2:N
            dvdt(ii) = (v(ii) - v(ii-1))/(dx*v(ii));
        end
        
        % ineq constraint 1
        Pcalc = (c1.*v.^2 + c2 + c3.*dvdt).*v; 
        P = s(N+1:end);
        c = [];
        if c > 0
            disp('c not satisfied')
        end
        
        % define delta 
        function d = delta(P,CP)
            d = P < CP;
        end
        
        delta_v = delta(P,CP);
        x2 = repmat(x(1:30),1,ceil(length(x)/30));
        shifted_x = x2(1:N);
        
        Wptot = Wcap;
        Wptot = Wcap - sum((1-delta_v+delta_v.*exp(-shifted_x./(v*tau_w))).*(P-CP)*dx);

        ceq = [Wptot P-Pcalc];
        
        if isnan(ceq)
            disp('AHHHHH')
            disp(ceq)
        end
    end
    
    %% fmincon arguments
    % create difference matrix
    A = eye(2*N);
    for n = 1:N
        if n==N
            row=[zeros(1,N) zeros(1,N)];
        else
            row = [zeros(1,N) zeros(1,N)];
            row(n) = 1;
            row(n+1) = -1;
        end
        A(n,:)=row;
    end
    
    b = [ones(1,N)*10 ones(1,N)*Pm];
    Aeq = [];
    beq = [];
    lb = zeros(1,2*N);
    ub = [ones(1,N)*25 ones(1,N)*Pm];
    ub(1) = 0.1;
   
    %% Call fmincon
    options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',1e8,'StepTolerance',1e-10,'MaxIterations',10e3);
    % to display iterations :'Display','iter'
    s = fmincon(fun,s0,A,b,Aeq,beq,lb,ub,@constraint,options);
    v = s(1:N);
    P = s(N+1:end);
    disp('Found solution')
    
end