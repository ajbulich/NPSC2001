clear

kappa=1;   % thermal conductivity
tau=5.5e-4;  % time-step
h=0.05;    % spatial-step
x=0:h:1;
y=0:h:1;


m = length(x);
D=-4*eye(m^2);
D=D+diag(ones(m^2-1,1),+1);
D=D+diag(ones(m^2-1,1),-1);
D=D+diag(ones(m^2-m,1),m);
D=D+diag(ones(m^2-m,1),-m);
D=kappa*tau/h^2*D;


D(1:m,:) = 0;
D(m^2-m:m^2,:)=0;

for l = 1:m-1
    D(l*m+1,:) = 0;
    D(l*m,:) = 0;
end
A = D + eye(m^2);
temp = zeros(m^2,1);
temp(round(m^2/2)) = 1/h^2;

time = 0;
while (time<=0.25)
    
    time=time+tau;   % update the time
    temp= A*temp;     % update the temperature
    i = i + 1;
    newTemp = reshape(temp, [21 21]);
    
    
end
mesh(x,y,newTemp)
drawnow

