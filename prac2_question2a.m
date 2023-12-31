clear;

h=0.1;       % timestep of 0.05 time-units
pos=[1 0];    % initial position
vel=[0 46];  % initial velocity

steps=400;              % number of steps
skip=ceil(steps/100);   % plot only 100 frames

theta=linspace(0,2*pi,100);   % equally spaced points
xc=cos(theta);
yc=sin(theta);

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
  accel= -r * pos/r;
 
  % Verlet Method
  if i==1
    next=pos + h*vel + 0.5*h*h*accel;
  else
    next=2*pos - prev + h*h*accel;
  end
 
  prev=pos;
  pos=next;
end

