clear; 

g=9.8;             % acceleration due to gravity is 9.8 m/s2
speed=40;          % initial speed is 40 m/s
angle=30;          % beware of radians vs degrees for trig functions

t=linspace(0,2*speed*sind(angle)/g,100);   % evenly spaced points
xp=speed*cosd(angle)*t;                    % x(t)
yp=speed*sind(angle)*t-0.5*g*t.^2;         % y(t)

h=0.5;             % timestep of 0.5 second
pos=[0 0];         % start at the origin
vel=speed*[cosd(angle) sind(angle)];

i=0; 
while pos(2)>=0
  i=i+1;
  x(i)=pos(1);
  y(i)=pos(2);
  
  % Midpoint Method
  pos=pos + h*vel + 0.5*h*h*[0 -g];
  vel=vel + h*[0 -g];
end

plot(xp,yp,x,y,'o')
legend('Parabola','Midpoint')
xlabel('distance (m)')
ylabel('height (m)')
