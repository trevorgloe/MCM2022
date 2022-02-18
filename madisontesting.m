clear
clc

%% ex 1
fun = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
x0 = [-1,2];
A = [1,2];
b = 1;
x = fmincon(fun,x0,A,b)

%% Ex
fun = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');


A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = @unitdisk;
x0 = [0,0];

% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)