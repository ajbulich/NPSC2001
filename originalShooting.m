function [value,X,PSI] = originalShooting(e)

x=-4;
h=0.01;
steps=800;

psi=1;
prev=0;

for i=1:steps
  X(i)=x;
  PSI(i)=psi;

  d2psi= x*x*psi - 2*e*psi;
  next=2*psi - prev + h*h*d2psi;

  prev=psi;
  psi=next;

  x=x+h;
end

PSI=PSI/max(abs(PSI));
plot(X,PSI)
title(e)
drawnow
func = PSI.*PSI;

for i=1:steps
    A(i) = h*func(i);
end
Area = sum(A);
PSI = PSI/sqrt(Area);

value=PSI(end);

