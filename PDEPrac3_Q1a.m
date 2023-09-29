clear
format compact

% Extent of the System
h=0.1;
x=0:h:1;
y=0:h:1;
m=length(x);

% Initial Conditions 

for ii=1:m
    for jj=1:m
        phi(jj, ii) = 4/(pi*sinh(pi))*sin(pi*x(ii))*sinh(pi*y(jj));
    end
end


% Boundary Conditions: phi=0 everywhere on the boundaries except for y=1
phi(m,:)=1;
phi(1,:)=0;
phi(:,m)=0;
phi(:,1)=0;





% Main loop
for n=1:1e4
  if rem(n,50)==0
    mesh(x,y,phi)
    title(n)
    drawnow
  end

  % Compute average of surrounding mesh points
  next=phi;
  for i=2:m-1
    for j=2:m-1
      next(i,j)=0.25*(phi(i-1,j)+phi(i,j-1)+ ...
                      phi(i+1,j)+phi(i,j+1));
    end
  end

  % Evaluate convergence criteria
  diff=abs(next-phi);
  if max(max(diff))<1e-6
    break
  end
    
  % Update to the new value of phi
  phi=next;
end

mesh(x,y,phi)
title(n)
disp(m)

