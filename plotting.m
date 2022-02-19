course = all_params{2};
disc = all_params{3};
x = linspace(0,course.L,disc.N);

figure()
plot(x,v,'.-','MarkerSize', 10)
xlabel('x position')
ylabel('v')
grid on

figure()
plot(x,P,'r')
xlabel('x position')
ylabel('P')