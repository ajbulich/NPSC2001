clear

kappa=1;   % thermal conductivity
tau=1e-4;  % time-step
h=0.02;    % spatial-step

% Extent of the system
x=0:h:1;
m=length(x);

% Matrix for the second-derivative term
D=-2*eye(m);
D=D+diag(ones(m-1,1),+1);
D=D+diag(ones(m-1,1),-1);
D=D*kappa*tau/(h*h);

% Boundary conditions (Dirichlet; i.e. constant value)
D(1,:)=0;
D(m,:)=0;

% The update matrix
A=eye(m) + D;

% Initial conditions
temp=zeros(m,1);
%temp(round(m/2))=1/h;

temp(1) = 0;
temp(1+1/h) = 5;

time=0;

kth = kappa*tau/h^2;

while (time<=1*0.3)
    fac = 4*kappa*time;
    analytic = 1/sqrt(pi*fac)*exp(-(x-0.5).^2/fac);
    
    if (time<=0.005)
        figure(1)
        plot(x,temp,'ro-')   % plot the temperature
        pause(0.08)
    end
   

    time=time+tau;   % update the time
    temp=A*temp;     % update the temperature
end

clf;
figure(2)

plot(x,temp,'ro-')   % plot the temperature
%plot(x,analytic, 'bo-')

legend('Numerical') %'Analytical')
title(kth)          % display current time



