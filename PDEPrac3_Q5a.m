clear
format compact
% extent of system
h = 0.1;
x = 0:h:1;
m = length(x);

% creating matrices

A = diag(0.5*ones(m-1,1),1);
A = A + diag(0.5*ones(m-1,1),-1);
A(1,:) = 0;
A(1,1) = 1;
A(m,:) = 0;
A(m,m) = 1;

phi = zeros(m, 1);


%initial conditions
phi(1) = 1;
phi(m) = 2;


for n=1:100000
    phi = A*phi;

    diff=mean(abs(phi - x'-1));
    if max(diff)<0.001
        break
    end
end
disp(n)
plot(x,phi, 'bo-')
disp(eig(A))

