clear
close all

h=0.01;       % timestep of 0.05 time-units
pos=[1 0];    % initial position
vel=[0 0.7];  % initial velocity

steps=1000;              % number of steps
skip=ceil(steps/100);   % plot only 100 frames

theta=linspace(0,2*pi,100);   % equally spaced points
xc=cos(theta);
yc=sin(theta);

for i=1:steps
  x(i)=pos(1);
  y(i)=pos(2);
  time=(i-1)*h/(2*pi);
  t(i) = time;

  if rem(i,skip)==0
    
    plot(xc,yc,'b',x,y,'g-',pos(1),pos(2),'ko',0,0,'ro')
    title(time)
    axis equal;
    pause(0);
  end

  r=norm(pos);
  accel=-1/r^2 * pos/r;
  sun_dist(i) = r;
 
  % Verlet Method
  if i==1
    next=pos + h*vel + 0.5*h*h*accel;
    v(i) = norm(vel);
  else
    next=2*pos - prev + h*h*accel;
    v(i) = norm((next-prev) / (2*h));
  end
  

  prev=pos;
  pos=next;

end

figure(2)
  subplot(2,1,1)
  plot(t, v, 'g-')
  xlabel('time')
  ylabel('speed')
  subplot(2,1,2)
  plot(t, sun_dist, 'r-')
  xlabel('time')
  ylabel('distance to sun')
  drawnow

vmax = max(v);
vmin = min(v);

dmax = max(sun_dist);
dmin = min(sun_dist);

Rv = vmax/vmin;
Rd = dmax/dmin;

disp(Rv)
disp(Rd)



