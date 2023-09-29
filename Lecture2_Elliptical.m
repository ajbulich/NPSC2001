clear

h = 0.01; %timestep, period is 2pi
pos = [1 0];
vel = [0 0.8];
steps = 1500;
skip = ceil(steps/100);

for i=1:steps
    x(i) = pos(1);
    y(i) = pos(2);
    time = (i-1)*h/(2*pi);

    total(i) = 0.5*norm(vel)^2 - 1/norm(pos);
    t(i) = time;

    %Circular Orbit
    theta = linspace(0,2*pi,100);
    xc = cos(theta);
    yc = sin(theta);
    
    if rem(i,skip)==0
        subplot(2,1,1)
        plot(xc,yc,'b',x,y,'g-', pos(1), pos(2), 'ko', 0, 0,'ro')
        title(time)
        axis equal
        drawnow %creates the animation

        subplot(2,1,2)
        plot(t,total)
    end
    
    r = norm(pos);
    acc =  - 1 / r^2 * pos / r;

    % Here goes the integrator - Midpoint Method
    %pos = pos + h*vel + 0.5*h^2*acc;
    %vel = vel + h*acc;

    % Next Method - Verlet Method
    %if i==1 
        %next = pos + h*vel + 0.5*h^2*acc; % Midpoint
    %else
        %next = 2*pos - prev + h^2*acc; % Verlet
    %end
    %prev = pos;
    %pos = next;

    % Next Method - Velocity-Verlet
    pos = pos + h*vel + 0.5*h^2*acc;
    r = norm(pos);
    accnext = - 1 / r^2 * pos / r;
    vel = vel + 0.5*h*(acc+accnext);


end
