clear;
format compact;    % type "help format" to learn what this does

g=9.8;             % acceleration due to gravity
speed=50;          % initial speed is 50 m/s
angle=40;          % beware of radians vs degrees for trig functions

h=1;      % timestep of 0.5 second
pos=[0 0];  % start at the origin
vel=speed*[cosd(angle) sind(angle)];

Cd = 0.35;
p = 1.2;
r = 0.037;
A = pi*r^2;
m = 0.145;

exactRange=speed^2*sind(2*angle)/g;

i=0;

while pos(2)>=0
  i=i+1;
  x(i)=pos(1);
  y(i)=pos(2);

  % Acceleration due to gravity with air resistance
  acc=  -(0.5*Cd*p*A*norm(vel)/m)*vel - [0 g];

  % Euler method
  pos=pos + h*vel + 0.5*h^2*acc;
  vel=vel + h*acc;
end

% plot the trajectory as points connected by lines
plot(x,y,'o-')
xlabel('distance (m)')
ylabel('height (m)')

% linear interpolation to estimate the range of the projectile
slope=(pos(2)-y(end))/(pos(1)-x(end));
range=x(end)-y(end)/slope;
perror = 100 * abs(range - exactRange) / exactRange;
disp('Range')
disp(range)
disp('% Error')
disp(perror)
