
course = all_params{2,qq};
disc = all_params{3,qq};
x = linspace(0,course.L,disc.N);

figure()
hold on
plot(x, v4,'.-','MarkerSize', 10)
plot(x,v,'.-','MarkerSize', 10)
xlabel('x position [km]')
ylabel('v [m/s]')
% title(['c',num2str(qq)])
grid on
hold off

figure()
hold on
plot(x, P4,'.-','MarkerSize', 10)
plot(x,P,'.-r','MarkerSize', 10)
xlabel('x position [km]')
ylabel('P [W]')
% title(['c',num2str(qq)])
grid on
hold off
% figure
% plot(sind(course.phi))