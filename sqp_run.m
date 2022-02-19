function [v,x] = sqp_run(course, biker, disc)
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
%     v0 = 20*linspace(0.1,1,N);
    v0 = 10*ones(1,N);
    x = linspace(0,L,N);
    phi_dis = interp1(linspace(0,L,length(phi)),phi,x);
    g = 9.8; % m/s
    
    c1 = 0.5*rho*A;% air drag coefficient *** get equation
    c2 = m.*g.*(sind(phi_dis) + Cr);
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
%         disp(max(P))
        if c > 0
            disp('c not satisfied')
        end
        % define delta 
        function d = delta(P,CP)
            d = P < CP;
        end
        
        delta_v = delta(P,CP);
%         for kk = 1:length(P)
%             if P(kk) > CP
%                 x_flag(kk) = 1;
%             else
%                 x_flag(kk) = 0;
%             end
%         end
%         x_flag = P > CP;
%         
%         top_x = x;
%         for kk = 2:length(P)
%             if x_flag(kk) == 0
%                 ii = kk;
%                 poop = x_flag(ii);
%                 while poop == 0 && ii > 0
%                     poop = x_flag(ii);
%                 	ind = (ii);
%                     ii = ii-1;
%                 end
%                 if ii == 0
%                     top_x(kk) = 0;
%                 else
%                     top_x(kk) = x(ind);
%                 end
%             end  
%         end        
        
%         %% create an array of points when the power goes from being blow CP to being above CP
%         changes = [];
%         for ii=2:N
%             if delta_v(ii-1) ~= delta_v(ii)
%                 changes = [changes ii];
%             end
%         end
        
        Wptot = Wcap;
%         Wexp = 0;
        Wptot = Wcap - sum((1-delta_v).*(P-CP)*dx);
        
%         if length(changes)==0
%             % no changes throughout the ride
%             if delta_v(ceil(N/2))==0
%                 Wptot = Wptot - sum((P-CP).*dx./v);
%             end
%         else
%             for int = 1:length(changes)+1
%                 if length(changes==1)
%                     if int==length(changes)+1
%                         Pint = P(changes(int-1):end);
%                         xint = x(changes(int-1):end);
%                         vint = v(changes(int-1):end);
%                     else
%                         Pint = P(1:changes(int));
%                         xint = x(1:changes(int));
%                         vint = v(1:changes(int));
%                     end
%                  else
%                     if int==length(changes)+1
%                         Pint = P(changes(int-1):end);
%                         xint = x(changes(int-1):end);
%                         vint = v(changes(int-1):end);
%                     else
%                         Pint = P(changes(int-1):changes(int));
%                         xint = x(changes(int-1):(changes(int)));
%                         vint = v(changes(int-1):changes(int));
%                     end
%                 end
%                 if int==length(changes)+1
%                     if delta_v(changes(int-1)+1) == 0
%                         Wexp = sum((Pint-CP).*dx./vint);
%                         Wptot = Wptot - Wexp;
%                     else
%                         Wptot = Wptot + Wexp*(1-exp((xint(end))./(vint(end)*tau_w)));
%                     end
%                 else
%                     if delta_v(changes(int)-1) == 0
%                         Wexp = sum((Pint-CP).*dx./vint);
%                         Wptot = Wptot - Wexp;
%                     else
%         %                 Wptot = Wptot - Wexp + Wexp*sum(exp((xint(1)-xint)./(vint*tau_w)).*dx./vint);
%                         Wptot = Wptot + Wexp*(1-exp((xint(end))./(vint(end)*tau_w)));
%                     end
%                 end
%             end
%         end
        
        ceq = Wptot;
        % eq constraint 2 
%         for jj = 1:N
%             innersum(jj) = sum( (1-delta(P(jj),CP)) * (P(jj)-CP) * (dx/v(jj)) );
%         end
%         for jj = 1:N
%             Wexp = Wexp + (1-delta(P(jj),CP))*(P(jj)-CP)*dx/v(jj);
%         end
%         for jj = 1:N
%             Wexp(jj) = sum((1-delta_v).*(P-CP).*dx./v);
%         end
%         ceq = Wcap - sum(Wexp.*exp(-delta_v.*(x-top_x./(v.*tau_w))).*dx./v);
% %         ceq = Wcap - sum(innersum.*exp(-delta(P,CP*ones(1,N)).*(top_x./(v.*tau_w))).*(dx./v)); 
% %         disp(Wexp)

        disp(ceq)
%         ceq = [];
%         if isnan(c)
%             disp('AHHHHH')
%             disp(c)

%         disp(sum(c<ones(1,N)*Pm))
        if isnan(ceq)
            disp('AHHHHH')
            disp(ceq)
%             disp(v)
        end
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