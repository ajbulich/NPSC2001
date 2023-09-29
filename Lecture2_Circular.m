clear

h = 0.1; %timestep, period is 2pi
pos = [1 0];
vel = [0 1];
steps = 150;
skip = ceil(steps/100);

for i=1:steps
    x(i) = pos(1);
    y(i) = pos(2);
    
    if rem(i,skip)==0
        plot(x,y,'g-', pos(1), pos(2), 'ko', 0, 0,'ro')
        axis equal
        drawnow %creates the animation
    end
    
    r = norm(pos);
    acc =  - 1 / r^2 * pos / r;

    % Here goes the integrator - Midpoint Method
    %pos = pos + h*vel + 0.5*h^2*acc;
    %vel = vel + h*acc;

    % Next Method - Verlet Method
    if i==1 
        next = pos + h*vel + 0.5*h^2*acc; % Midpoint
    else
        next = 2*pos - prev + h^2*acc; % Verlet
    end
    prev = pos;
    pos = next;

end
