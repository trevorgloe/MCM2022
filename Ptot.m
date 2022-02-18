function Ptot = Ptot(Pa,params,v,psi,delta)
%PTOT Calculates the total power outputting based on the velocit, psi at
%that point and delta

Fn = params.m*params.g/cos(delta);
Pd = params.Cd*v^2;
Prr = Fn*params.Crr*v*params.Cs*cos(psi);
Pbr = Fn*params.Cbr*params.rb*v/params.rw;
Pg = params.m*params.g*v;

Pdem = Pd + Prr + Pbr + Pg;
Psup = 0.98*Pa;

Ptot = Psup - Pdem;


end

