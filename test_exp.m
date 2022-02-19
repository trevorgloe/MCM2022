% test of comparision of integral of exponential to actual exponential
A = 1;
x = linspace(0,5000,N);


figure
hold on
plot(A*(1-exp(-x/tau_w)))