clear;
format compact;    % type "help format" to learn what this does

g=9.8;             % acceleration due to gravity
speed=0;          % initial speed is 50 m/s
angle=270;          % beware of radians vs degrees for trig functions

h=0.001;      % timestep of 0.5 second
bigPos=[0 50];  % start at 50m above ground
smallPos=[0 50];
bigVel=speed*[cosd(angle) sind(angle)];
smallVel = bigVel;


Cd = 0.018;
p = 1.2;
rho = 7800;
bigM = 100 * 0.454;
smallM = 0.454;
bigR = ( (3*bigM) / (4*pi*rho) )^(1/3);
smallR = ( (3*smallM) / (4*pi*rho) )^(1/3);
bigA = pi*bigR^2;
smallA = pi*smallR^2;

exactRange=speed^2*sind(2*angle)/g;

i=0;
time = 0;

while bigPos(2)>=0
  i=i+1;
  x(i)=bigPos(1);
  y(i)=bigPos(2);
  x2(i)=smallPos(1);
  y2(i)=smallPos(2);

  % Acceleration due to gravity with air resistance
  bigAcc=  -(0.5*Cd*p*bigA*norm(bigVel)/bigM)*bigVel - [0 g];
  smallAcc= -(0.5*Cd*p*smallA*norm(smallVel)/smallM)*smallVel - [0 g];

  % Euler method
  bigPos= bigPos + h*bigVel;
  bigVel= bigVel + h*bigAcc;
  
  smallPos = smallPos + h*smallVel;
  smallVel = smallVel + h*smallAcc;
  
  time = time + h;
end
if smallPos(2) < 0
    smallPos(2) = 0;
end

if bigPos(2) < 0
    bigPos(2) = 0;
end
% plot the trajectory as points connected by lines
plot(x,y,'o-')
xlabel('distance (m)')
ylabel('height (m)')

disp(bigPos)
disp(smallPos)

% linear interpolation to estimate the range of the projectile
%slope=(pos(2)-y(end))/(pos(1)-x(end));
%range=x(end)-y(end)/slope;
%perror = 100 * abs(range - exactRange) / exactRange;

