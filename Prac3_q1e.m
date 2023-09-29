clear;

g=9.8;             % acceleration due to gravity
L=0.2482;          % length in metres (gives period of 1 second)

theta=45*pi/180;   % initial angle (in radians)
omega=0;           % initial angular-velocity (in radians/second)
h=0.001;            % timestep (in seconds)
steps=round(3.5/h);  % run for 3.5 seconds
time = 0;

for angle=1:170
    theta = angle * pi/180;
    omega = 0;
    time = 0;
    while theta > 0
      %tList(i) = theta;
      time = time + h;
      
      alpha=-g/L*sin(theta);     % acceleration term
    
      % Midpoint Method
      % theta=theta + h*omega + 0.5*h*h*alpha;
      % omega=omega + h*alpha;
    
      
     
      % Velocity-Verlet Method
      theta = theta + h*omega + 0.5*h^2*alpha;
      alphanext = -g/L*sin(theta);
      omega = omega + 0.5*h*(alpha + alphanext);
    end 
    period(angle) = time * 4;
end

figure(2)
plot(period, 'b-')

xlabel("Initial angle (degrees)")
ylabel("Period (seconds)")
ylim([1 1.05])

