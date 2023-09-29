clear;

h=0.0001;       % timestep of 0.05 time-units
pos=[1 0];    % initial position
vel=[0 1.1];  % initial velocity

steps=10000;              % number of steps
skip=ceil(steps/100);   % plot only 100 frames

theta=linspace(0,2*pi,100);   % equally spaced points
xc=cos(theta);
yc=sin(theta);
originalh = h;

for i=1:steps
  x(i)=pos(1);
  y(i)=pos(2);
  time=(i-1)*h/(2*pi);

  if rem(i,skip)==0
    plot(xc,yc,'b',x,y,'g-',pos(1),pos(2),'ko',0,0,'ro')
    title(time)
    axis equal;
    pause(0);
  end

  r=norm(pos);
  acc=-1/r^2 * pos/r;
 
  % Verlet Method
  if i <= steps/2
      pos = pos + h*vel + 0.5*h^2*acc;
      vel = vel + h*acc;
  else
      pos = pos - h*vel + 0.5*h^2*acc;
      vel = vel - h*acc;
  end
  
end

disp(pos)

