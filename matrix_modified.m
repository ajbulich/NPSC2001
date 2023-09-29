clear

h=0.2;
x=-4:h:4;
m=length(x);
kList = [1,2,3,4];
cList = ['b-', 'g-', 'r-', 'k-'];
hold on

for i=1:4
    D=2*eye(m);
    D=D-diag(ones(m-1,1),+1);
    D=D-diag(ones(m-1,1),-1);
    D=D/(2*h^2);
    
    V=0.5*kList(i)*diag(x.^2);
    
    A=D+V;
    [vmat, emat] = eig(A);
    plot(x, vmat(:,1), cList(i))
end
hold off


