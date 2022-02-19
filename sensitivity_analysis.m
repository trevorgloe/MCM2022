% analyze results of sensitivity study

m = 9;      % number of data points for the parameter
% first parse all data collected
all_P = zeros(m,N);
all_v = zeros(m,N);
all_m = zeros(1,m)

cd model_data

Folders = dir(fullfile('Wcap_sensitivity','*.*'));
cd Wcap_sensitivity

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
    all_m(s)=biker.Wcap;
    
    cd ..
    
end

x = linspace(0,L,N);

figure
hold on
for s = 1:m
    plot(x,all_v(s,:))
end
xlabel('x [m]')
ylabel('V [m/s]')

figure
hold on
for s = 1:m
    plot(x,all_P(s,:),'Linewidth',2)
end
labs = arrayfun(@num2str, all_m, 'UniformOutput',0);
[hleg, att] = legend(labs);
title(hleg, "W'_{capacity} values")
hleg.Title.Visible = 'on';
hleg.Title.NodeChildren.Position = [0.5 1.1 0];
xlabel('x [m]')
ylabel('P [Watts]')