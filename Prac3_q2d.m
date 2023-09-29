clear;

global g L         % define g and L as global variables
g=9.8;             % acceleration due to gravity
L=0.2482;          % length in metres (gives period of 1 second)

theta=90*pi/180;   % intial angle (radians)
omega=0;           % intital angular-velocity

x=[theta omega];   % create vector x
h=0.1;            % timestep of 0.01 seconds
steps = 50/h;
m = 1;

% Second-order Runge-Kutta (beta=1)
for i=1:steps
  %Velocity-Verlet  
  alpha=-g/L*sin(x(1));  
  x(1) = x(1) + h*x(2) + 0.5*h^2*alpha;
  alphanext = -g/L*sin(x(1));
  x(2) = x(2) + 0.5*h*(alpha + alphanext);  

  % Rk2  
  % F1=deriv(x);
  % F2=deriv(x+h*F1);
  % x=x+h/2*(F1+F2);

  %RK4
  %F1 = deriv(x);
  %F2 = deriv(x+0.5*h*F1);
  %F3 = deriv(x+0.5*h*F2);
  %F4 = deriv(x+h*F3); 
  %x = x + h/6 * (F1+2*F2+2*F3+F4);

  tList(i) = x(1);
  % Graphic of the pendulum

  %xpend=[0  L*sin(x(1))];
  %ypend=[0 -L*cos(x(1))];
  %plot(xpend,ypend,'o-')

  %axis equal
  %axis([-L L -L L])
  %pause(h)
  E(i) = 0.5*m*L^2*x(2)^2 - m*g*L*cos(x(1));
end 

figure(2)
plot(1:steps, tList)
figure(3)
plot(1:steps, E)


