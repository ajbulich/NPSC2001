clear;

h=0.001;       % timestep of 0.05 time-units
pos=0;    % initial position
vel=0.9;  % initial velocity
beta=0.2;
w = 1;
h2 = 0;


steps=25000;              % number of steps
skip=ceil(steps/100);   % plot only 100 frames
for i=1:steps

  x(i) = (sqrt(31))/(2*sqrt(6))*exp(-0.2*h2)*cos((2*sqrt(6))/(5)*h2-atan((5*sqrt(6))/(6)));
  v(i) = (sqrt(31))/(2*sqrt(6)) * (-0.2*exp(-0.2*h2)*cos((2*sqrt(6))/(5)*h2-atan((5*sqrt(6))/(6))) + exp(-0.2*h2)*(-2*sqrt(6)/5 * sin((2*sqrt(6))/(5)*h2-atan((5*sqrt(6))/(6)))));
  h2 = h2 +h;
end

plot(1:steps, x, 'g-')
xlabel("time (ms)")
ylabel("x(t)")

figure(2)
plot(1:steps, v, 'b-')
xlabel("time (ms)")
ylabel("v(t)")

figure(3)
plot(x,v,'r-')
ylabel("v(t)")
xlabel("x(t)")
xlim([-1 1])
ylim([-1 1])


