clear
clf

h=0.25;
x=-4:h:4;
m=length(x);

D=2*eye(m);
D=D-diag(ones(m-1,1),+1);
D=D-diag(ones(m-1,1),-1);
D=D/(2*h^2);

potential = 10*ones(m,1);
potential(abs(x)<1)=0;
V=diag(potential);

A=D+V;

[vmat,emat] = eig(A);
B = inv(A);

psi = ones(m,1);

for i=1:30
    psi=A*psi;
    psi = psi/max(abs(psi));
end
plot(x,psi,'r-')



