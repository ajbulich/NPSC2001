clear;

h=0.05;       % timestep of 0.05 time-units
pos=[1 0 0 ];    % initial position
vel=[0 1/(sqrt(2)) -1/(sqrt(2))];  % initial velocity

steps=300;              % number of steps
skip=ceil(steps/100);   % plot only 100 frames

theta=linspace(0,2*pi,100);   % equally spaced points
xc=cos(theta);
yc=sin(theta);

for i=1:steps
  x(i)=pos(1);
  y(i)=pos(2);
  z(i)= pos(3);
  time=(i-1)*h/(2*pi);

  if rem(i,skip)==0
    plot3(x,y,z,'g-',pos(1),pos(2),pos(3),'ko',0,0,0,'ro')
    title(time)
    axis equal;
    pause(0);
  end

  r=norm(pos);
  v = norm(vel);
  k = 0.002;
 
  % Euler-Richardson Method

  acc=(-1 + k*v^2)/r^2 * pos/r;
  vmid = vel + 0.5*acc*h;
  xmid = pos + 0.5*vel*h;
  v2 = norm(vmid);
  r2 = norm(xmid);
  amid = (-1 + k*v2^2)/r2^2 * xmid/r2;
  vel = vel + amid*h;
  pos = pos + vmid*h;
end

disp(pos)

