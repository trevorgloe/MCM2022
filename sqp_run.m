function [v] = sqp_run(course, biker, disc)
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
    v0 = ones(1,N);
    x = linspace(0,L,N);
    phi_dis = interp1(linspace(0,L,length(phi)),phi,x);
    g = 9.8; % m/s
    
    c1 = 0.5*rho*A;% air drag coefficient *** get equation
    c2 = m.*g.*(phi_dis + Cr);
    c3 = m; % get eq ***
    
    %% F function
    fun = @(v) dx.*sum(1./v);
    
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
        
        % define delta 
        function d = delta(P,CP)
            if P < CP
                d = 1;
            else
                d = 0;
            end
        end
        
        % eq constraint 2 
        for jj = 1:N
            innersum(jj) = sum( (1-delta(P(jj),CP)) * (P(jj)-CP) * (dx/v(jj)) );
        end
         ceq = Wcap - sum(innersum.*exp(-delta(P,CP*ones(1,N)).*(x./(v.*tau_w))).*(dx./v)); 
    end
    
    %% fmincon arguments
    A = eye(N);
    b = ones(1,N)*1000;
    Aeq = [];
    beq = [];
    lb = zeros(1,N);
    ub = [];
    
    nonlcon =@ constraint; %,x,Pm,N,Wcap,CP,c1,c2,c3,tau_w);
    
    %% Call fmincon
    options = optimoptions('fmincon','Algorithm','sqp');
    v = fmincon(fun,v0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end