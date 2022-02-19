%% Solve using fmincon again

v0 = 10*ones(1,N);

% make difference matrix
A = eye(N);
for n = 1:N
    if n==N
        row=zeros(1,N);
    else
        row = zeros(1,N);
        row(n) = 1;
        row(n+1) = -1;
    end
    A(n,:)=row;
end
b = 2*ones(1,N);
Aeq = [];
beq = [];
lb = zeros(1,N);
ub = 30*ones(1,N);

fun = @(v)objfun(v,dx);

fmincon(fun,v0,A,b,Aeq,beq,lb,ub,@confuneq)


function f = objfun(v,dx)
    f = sum(1./v*dx);
end

function [c,ceq] = confuneq(v)
    dvdt = ones(size(v));
    dvdt(1) = 0;
    dvdt(2:end) = (v(1:end-1) - v(2:end))./(v(1:end-1)*dx);
    
    P = (c1*v + c2 + c3*dvdt).*v;
    
    c = P - Pm*ones(1,N);
    delta = P < CP;
    
    ceq = Wcap - sum((1-delta).*(P-CP)*dx);
    
end