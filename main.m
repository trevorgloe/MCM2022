clear
close all
clc

plotflag = 1;

params_discret
sqp_solve
convert_v

all_params = {biker, course, disc};
if plotflag == 1
    plotting
end


save_data_params(all_params,v,P);