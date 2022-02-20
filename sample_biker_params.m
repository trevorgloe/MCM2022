%% Sample parameter space of CP, tau_w, Wcap, and Pm

n1 = 3;     %points from CP
n2 = 2;     %points from Wcap
n3 = 4;     %points from tau_w
n4 = 5;     %points from Pm

cp_minmax = [140 205];
Wcap_minmax = [14e3 30e3];
tau_w_minmax = [60 500];
Pm_minmax = [285 400];

cp_vec = linspace(cp_minmax(1),cp_minmax(2),n1);
Wcap_vec = linspace(Wcap_minmax(1),Wcap_minmax(2),n2);
tau_w_vec = linspace(tau_w_minmax(1),tau_w_minmax(2),n3);
Pm_vec = linspace(Pm_minmax(1),Pm_minmax(2),n4);


% create mesh
for i=1:n1
    for j=1:n2
        for k=1:n3
            for n=1:n4
                ii = sub2ind([n1 n2 n3 n4],i,j,k,n);
                param_mesh(ii).CP = cp_vec(i);
                param_mesh(ii).Wcap = Wcap_vec(j);
                param_mesh(ii).tau_w = tau_w_vec(k);
                param_mesh(ii).Pm = Pm_vec(n);
            end
        end
    end
end
    
    
    
    