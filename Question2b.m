clear;
format compact;    % type "help format" to learn what this does

g=9.8;             % acceleration due to gravity
speed=50;          % initial speed is 50 m/s
angle=5;          % beware of radians vs degrees for trig functions

h=0.001;      % timestep of 0.5 second



Cd = 0.35;
p = 1.2;
r = 0.037;
A = pi*r^2;
m = 0.145;

exactRange=speed^2*sind(2*angle)/g;


ii = 0;
while ii <= 80
    ii = ii + 1;
    vel=speed*[cosd(angle) sind(angle)];
    i=0;
    pos=[0 0];  % start at the origin
    while pos(2)>=0
      i=i+1;
      x(i)=pos(1);
      y(i)=pos(2);
    
      % Acceleration due to gravity with air resistance
      acc=  -(0.5*Cd*p*A*norm(vel)/m)*vel - [0 g];
    
      % Euler method
      pos=pos + h*vel;
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
    
    angles(ii) = angle;
    ranges(ii) = range;
    angle = angle + 1;
end
disp(ranges)
plot(angles,ranges, 'o')
