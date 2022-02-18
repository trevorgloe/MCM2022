function [output] = sqp_run(course, biker)
% takes in rider and course struct and runs SQP based model
    
    %% Get params
    %***
    dx = L/N;
    vvals
    
    %% F function
    fun = @(v) dx.*sum(1./v);
    % function needs to return F and dFdt?
    
    %% Constratint funct
    function [P, Wcap] = power_constraint(v,Pm,N)
        % Backward difference approx for dvdt
        dvdt(1) =  (v(ii))/(dx*v(ii));
        for ii = 2:N
            dvdt(ii) = (v(ii) - v(ii-1))/(dx*v(ii));
        end
        
        P = c1.*v + c2 + c3.*dvdt; % eq constraint
        
        if P > Pm % ineq constraint 1
            P = Pm;
        end
        
        % define delta 
        function d = delta(P,CP)
            if P < CP
                d = 1;
            else
                d = 0;
            end
        end
        
        for jj = 1:N
            innersum(jj) = sum( (1-delta(P(jj))) * (P(jj)-CP) * (dx/v(jj)) );
        end
        
        Wcap - sum(innersum*exp(-delta*(x./(v.*tau_w)))*(dx./v)) == 0; % ineq constraint 2 
    end
    
    %% fmincon arguments
    A = eye(N);
    b = ones(1,N)*1000;
    Aeq = [];
    beq = [];
    lb = zeros(1,N);
    ub = [];
    
    nonlcon = @power_constraint;
    
    %% Call fmincon
    options = optimoptions('fmincon','Algorithm','sqp');
    x = fmincon(fun,xvals,A,b,Aeq,beq,lb,ub,nonlcon,options);
end