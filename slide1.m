clear
format compact

h=0.015;       % timestep of 0.015 time-units
pos=[1 0];     % initial position
vel=[0 0];   % initial velocity

steps=10000;            % number of steps
skip=ceil(steps/100);  % plot only 100 frames

theta=linspace(0,2*pi,100);   % equally spaced points
xc=cos(theta);
yc=sin(theta);

time=0;
for i=1:steps
  x(i)=pos(1);
  y(i)=pos(2);
  time=time+h;

  % Plot every skip'th frame
  if rem(i,skip)==0
    plot(xc,yc,'b',x,y,'g-',pos(1),pos(2),'ko',0,0,'ro')
    title(time/(2*pi))
    axis equal
    drawnow
  end

  % Velocity-Verlet Method
  r=norm(pos);
  accel=-1/r^2 * pos/r;
  pos=pos + h*vel + 0.5*h*h*accel;

  r=norm(pos);
  accelnext=-1/r^2 * pos/r;
  vel=vel + 0.5*h*(accel+accelnext);

  % Variable time-step
  deltav=0.5*h*(accel+accelnext);
  fraction=norm(deltav)/norm(vel);
  if fraction>0.02
    h=h*0.5;
  elseif fraction<0.01
    h=h*2;
  end
  hList(i) = h;
  rList(i) = r;
  E(i) = 0.5*1*(norm(vel))^2 - 1/(norm(pos));

  % Terminate if h too small/large
  if h<1e-6 || h>1000
    break
  end
end

days = (time/(2*pi)) * 365;
disp(days)

figure(2)
semilogy(hList, rList, 'b-')
xlabel("Distance from the sun")
ylabel("Timestep (log(s))")

figure(3)
plot(rList, E, 'b-')
xlabel("distance from sun")
ylabel('Energy')

