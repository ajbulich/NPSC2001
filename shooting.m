function [value,X,PSI] = shooting(e,h)

x=-4;
steps=round(8/h);

psi=1;
prev=0;


for i=1:steps
  X(i)=x;
  PSI(i)=psi;

  v = x^4/10;
 
  d2psi= 2*psi*(v-e);
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

