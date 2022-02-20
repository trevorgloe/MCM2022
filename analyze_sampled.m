% analyze data that has been generated from a sample set of parameters
load('data_run_short_flat');

t_array = zeros(1,length(data_run));
avgP_array = zeros(1,length(data_run));
avgv_array = zeros(1,length(data_run));

for ii=1:length(data_run)
    t_array(ii) = data_run(ii).time;
    avgP_array(ii) = mean(data_run(ii).P);
    avgv_array(ii) = mean(data_run(ii).v);
end


