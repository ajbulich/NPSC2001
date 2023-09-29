clear
format compact

% Extent of the System
h=0.2;
x=0:h:1;
y=0:h:1;
m=length(x);

omega=1:0.05:1.95;
omega = 2/(1+sin(pi/m));


% Main loop
ww=1;
for o = omega
    % Initial Conditions 
    phi=zeros(m);
    
    % Boundary Conditions: phi=0 everywhere except for y=1
    phi(m,:)=1;

    for n=1:1e4
      if rem(n,50)==0
        mesh(x,y,phi)
        title(n)
        drawnow
      end
    
      % Compute average of surrounding mesh points
      old=phi;
      for i=2:m-1
        for j=2:m-1
         phi(i,j)=(1-o)*phi(i,j) + o*0.25*(phi(i-1,j)+phi(i,j-1)+ ...
                          phi(i+1,j)+phi(i,j+1));
        end
      end
    
      % Evaluate convergence criteria
      diff=abs(old-phi);
      if max(max(diff))<1e-6
        break
      end   
    end
    nList(ww) = n;
    ww = ww + 1;
end
mesh(x,y,phi)
title(n)
disp(m)




