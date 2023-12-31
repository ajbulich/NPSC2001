clear;

g=9.8;             % acceleration due to gravity
L=0.2482;          % length in metres (gives period of 1 second)

theta=10*pi/180;   % initial angle (in radians)
omega=0;           % initial angular-velocity (in radians/second)
h=0.04;            % timestep (in seconds)
steps=round(3.5/h);  % run for 3.5 seconds


for i=1:steps
  % Graphic of the pendulum
  xpend=[0  L*sin(theta)];
  ypend=[0 -L*cos(theta)];
  plot(xpend,ypend,'o-')
  tList(i) = theta;

  axis equal
  axis([-L L -L L])
  pause(h)

  alpha=-g/L*sin(theta);     % acceleration term

  % Midpoint Method
  theta=theta + h*omega + 0.5*h*h*alpha;
  omega=omega + h*alpha;
end 
figure(2)
plot(1:steps, tList, 'o-')

