clear

g = 9.8;
speed = 40;
angle = 30;

t = linspace(0,2*speed*sind(angle)/g, 100);

xp = speed * cosd(angle)*t;
yp = speed * sind(angle)*t - 0.5*g*t.^2;

plot(xp,yp)
legend('Parabola')
xlabel('distance (m)')
ylabel('height (m)')