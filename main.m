clear
close all
clc

plotflag = 1;

params_discret
sqp_solve
convert_v

if plotflag == 1
    plotting
end

all_params = {rider, course, disc};
save_data_params(all_params,v,P);