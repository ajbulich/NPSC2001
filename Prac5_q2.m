clear
clf

h=0.05;
x=-4:h:4;
m=length(x);

D=2*eye(m);
D=D-diag(ones(m-1,1),+1);
D=D-diag(ones(m-1,1),-1);
D=D/(2*h^2);

potential = 100*ones(m,1);
potential(abs(x)<1)=0;
V=diag(potential);

A=D+V;

[vmat,emat] = eig(A);
cList = ['b-', 'g-', 'r-', 'm-', 'k-'];
hold on
for i=1:5
    disp(emat(i,i));
    plot(x, vmat(:,i))
end
legend('E1', 'E2', 'E3', 'E4', 'E5')
hold off
figure(2)
bar(eig(A))
ylabel("Eigenvalues")


