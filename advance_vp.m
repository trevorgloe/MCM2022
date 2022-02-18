function ds = advance_vp(t,s,params)
%ADVANCE_VP advances the state vector based on the dynamical system
% the state s = [v x]
v = s(1);
x = s(2);

dx = v;

delta = interp1(params.psi_xvals,params.psi,x);
P = interp1(params.P_xvals,params.Papplied,x);
dv = Ptot(P,params,v,delta,delta);

ds = [dv dx]';



end

