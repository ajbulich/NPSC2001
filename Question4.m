clear;
format compact;    % type "help format" to learn what this does

g=9.8;             % acceleration due to gravity
speed1=0; 
speed2 = 50;        % initial speed is 50 m/s
angle2 = 0;
angle1=270;          % beware of radians vs degrees for trig functions

h=0.001;      % timestep of 0.5 second
pos1=[0 50];  % start at 50m above ground
pos2 =[0 50];
vel1 = speed1*[cosd(angle1) sind(angle1)];
vel2 = speed2*[cosd(angle2) sind(angle2)];


Cd = 0.3;
p = 1.2;
rho = 7800;
m = 1 * 0.454;
R = ( (3*m) / (4*pi*rho) )^(1/3);
A = pi*R^2;

i=0;
time = 0;

while pos1(2)>=0
  i=i+1;
  x(i)=pos1(1);
  y(i)=pos1(2);
  x2(i)=pos2(1);
  y2(i)=pos2(2);

  % Acceleration due to gravity with air resistance
  acc1=  -(0.5*Cd*p*A*norm(vel1)/m)*vel1 - [0 g];
  acc2= -(0.5*Cd*p*A*norm(vel2)/m)*vel2 - [0 g];

  % Euler method
  pos1= pos1 + h*vel1;
  vel1= vel1 + h*acc1;
  
  pos2 = pos2 + h*vel2;
  vel2 = vel2 + h*acc2;
  
  time = time + h;
end
if pos2(2) < 0
    pos2(2) = 0;
end

if pos1(2) < 0
    pos1(2) = 0;
end
% plot the trajectory as points connected by lines
xlabel('distance (m)')
ylabel('height (m)')
xlim([-1 200])


plot(x,y,x2,y2,'o-')



disp(pos1)
disp(pos2)

% linear interpolation to estimate the range of the projectile
%slope=(pos(2)-y(end))/(pos(1)-x(end));
%range=x(end)-y(end)/slope;
%perror = 100 * abs(range - exactRange) / exactRange;

