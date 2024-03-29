clear
clc

%% Test 1
L = 10;
k = 20;
dx = L/k;
fun = @(v)dx.*sum(1./v);
x0 = 20*ones(1,k);
A = 20*ones(1,20);
b = 1;
x = fmincon(fun,x0,A,b);

%% 
clear
params_discret



%% Ex
% fun = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% 
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% lb = [];
% ub = [];
% nonlcon = @unitdisk;
% x0 = [0,0];

% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)