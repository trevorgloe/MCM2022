%% solve optimization problem using problem based approach in matlab

v = optimvar('v',1,N,'UpperBound',30,'LowerBound',0);

obj = fcn2optimexpr(@time,v,dx);

prob = optimproblem('Objective',obj);

nlcon1 = u(v,c1,c2,c3,x,dx,N) <= 350;
% prob.Constaints.eq = G(v,c1,c2,c3,x,dx,N,CP) <= Wcap*ones(1,N);
prob.Constraints.eq = sum(delta(P,CP).*u(v,c1,c2,c3,x,dx,N)-CP)' == Wcap;

prob.Constraints.leq = nlcon1
% prob.Constraints.eq = nlcon2

vguess.v = 10*ones(1,N);
[sol,fval,existflag,output] = solve(prob,vguess,'Solver','quadprog')

function Ttot = time(v,dx)
    Ttot = sum(1./v*dx);
end

function P = u(v,c1,c2,c3,x,dx,N)
    dvdt(1) =  (v(1))/(dx*v(1));
    for ii = 2:N
        dvdt(ii) = (v(ii) - v(ii-1))/(dx*v(ii));
    end

    P = (c1.*v + c2 + c3.*dvdt).*v;
end

function del = delta(P,CP)
%     del = P <= CP;
    ind = find(P<CP);
    del = ones(size(P))
    del(ind) = 0;
end

function Wexp = G(v,c1,c2,c3,x,dx,N,CP)
    P = u(v,c1,c2,c3,x,dx,N);
    
    delta = P < CP;
    Wexp = sum((P-CP)*(1-delta)*dx);
end

function mod_u = check(u,CP)
end
    