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
temp(round(m/2))=1/h;

time=0;
i = 1;
while (time<=1)
    fac = 4*kappa*time;
    analytic = 1/sqrt(pi*fac)*exp(-(x-0.5).^2/fac);
    halfA(i) = analytic((length(x)-1)/2);
    newtemp(i) = temp((length(x)-1)/2);



    time=time+tau;   % update the time
    temp=A*temp;     % update the temperature
    i = i + 1;
end

figure(1)

semilogy(0:tau:1, newtemp, 'ro-', 0:tau:1, halfA, 'bo-')
legend('Numerical', 'Analytic')
hold off

figure(2)
plot(x,temp,'ro-')

