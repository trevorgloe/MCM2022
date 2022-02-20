
course = all_params{2,qq};
disc = all_params{3,qq};
x = linspace(0,course.L,disc.N);

figure()
plot(x,v,'.-','MarkerSize', 10)
xlabel('x position')
ylabel('v')
title(['c',num2str(qq)])
grid on

figure()
plot(x,P,'r')
xlabel('x position')
ylabel('P')
title(['c',num2str(qq)])

figure
plot(course.phi)