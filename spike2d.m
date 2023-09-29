clear;

kappa=1;   % thermal conductivity
tau=6.25e-4;  % time-step
h=0.05;    % spatial-step

% Extent of the system
x=0:h:1;
y=0:h:1;
m=length(x);

% Initial conditions
temp=zeros(m,m);
temp(round(m/2),round(m/2))=1/h^2;

time=0;
for n=1:100
    mesh(x,y,temp)      % mesh plot
    title(n*tau);       % display current time
    drawnow

    next=temp;
    for i=2:m-1
      for j=2:m-1
        factor= temp(i-1,j)+temp(i+1,j)+temp(i,j-1)+temp(i,j+1) -4*temp(i,j);
        next(i,j)=temp(i,j) + kappa*tau/h^2 * factor;
      end
    end
    temp=next;
end
disp(kappa*tau/h^2)
