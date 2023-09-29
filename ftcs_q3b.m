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
while (time<=0.25)
    fac = 4*kappa*time;
    analytic = 1/sqrt(pi*fac)*exp(-(x-0.5).^2/fac);
    halfA(i) = analytic((length(x)-1)/2);
    newtemp(i) = temp((length(x)-1)/2);



    time=time+tau;   % update the time
    temp=A*temp;     % update the temperature
    temp(1) = temp(2);
    i = i + 1;
end

[V, D] = eig(A);

% extraction of D into a list
d = length(D);
for ii = 1:d
    eVal(ii) = D(ii,ii);
end
eVec1 = V(:,1);
eVec2 = V(:,2);
eVec3 = V(:,end-1);
eVec4 = V(:,end);

figure(1)
plot(x,temp,'ro-')

figure(2)
plot(x,eVec1,'ro-', x, eVec2,'bo-', x, eVec3,'yo-', x, eVec4, 'mo-')
legend('Vec1', 'Vec2', 'Vec3', 'Vec4')

figure(3)
semilogy(0:tau:0.25, newtemp, 'ro-')
