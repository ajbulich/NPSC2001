clear

h = 0.01;
x = -5:h:5;
m = length(x)

D = 2*eye(m);
D = D - diag(ones(m-1,1),+1);
D = D - diag(ones(m-1,1),-1);
D = 0.5*D/h^2;

V = 0.5*diag(x.^2);

A = D + V;

[vmat, emat] = eig(A);

plot(x, vmat(:,2), 'ro-')
disp(emat(2,2))
figure(2)
plot(x, vmat(:,2).^2, 'ro-')


