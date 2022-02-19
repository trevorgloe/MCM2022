% analyze results of sensitivity study

m = 9;      % number of data points for the parameter
% first parse all data collected
all_P = zeros(m,N);
all_v = zeros(m,N);
all_m = zeros(1,m)

cd model_data

Folders = dir(fullfile('CP_sensitivity','*.*'));
cd CP_sensitivity

for s = 1:m
    Files = dir(fullfile(Folders(s+2).name,'*.mat'))
    disp(Folders(s+2))
    old = cd(Folders(s+2).name);
%     old = cd(Folders(s+2).name)
    load(Files(1).name)
    load(Files(2).name)
    load(Files(3).name)
    
    all_P(s,:)=P;
    all_v(s,:)=v;
    biker = all_params{1};
    all_m(s)=biker.m;
    
    cd ..
    
end

figure
hold on
for s = 1:m
    plot(all_v(s,:))
end
xlabel('x')
ylabel('V [m/s]')

figure
hold on
for s = 1:m
    plot(all_P(s,:))
end
xlabel('x')
ylabel('P [Watts')