clear
format compact

% Extent of the System
h=0.05;
x=0:h:1;
y=0:h:1;
m=length(x);

% Initial Conditions 
phi=zeros(m);

% Boundary Conditions: phi=0 everywhere except for y=1
phi(m,:)=1;

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

