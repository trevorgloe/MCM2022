function Wexp = calcWexp(v,c1,c2,c3,Wcap,tau_w,CP,L,dx,N,x,q)
% Backward difference approx for dvdt
    dvdt(1) =  (v(1))/(dx*v(1));
    for ii = 2:N
        dvdt(ii) = (v(ii) - v(ii-1))/(dx*v(ii));
    end

    P = (c1.*v + c2 + c3.*dvdt).*v;
    % get x values of changes

    % integral of modified lagrangian with lagrangian multiplier q
    delta = P < CP;
    
    Lag = 1./v;
    
    G = (1-delta).*(P-CP)./v;
    
    vals = cumtrapz(x,G);
    
    Wexp = vals(end);

end