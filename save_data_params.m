function [] = save_data_params(all_params,v,P)
%SAVE_DATA_PARAMS saves the data and parameters with their own unique name
%in a certain folder where they will all go
%   all_params is a cell array of the 3 parameter structs (biker, course,
%   and disc (discretization))
% d = datetime('Now');
formatOut = 'mm-dd-yy_HH-MM-SS';
base_name = datestr(now,formatOut);

cd model_data
status = mkdir(base_name);

save([base_name '/' base_name '_P_'],'P')
save([base_name '/' base_name '_v_'],'v')
save([base_name '/' base_name '_params_'],'all_params')

disp('Data saved!')

cd ..


end

