clear
format long
% extent of system
h = 0.1;
x = 0:h:1;
m = length(x);

% creating matrices
A = zeros(m,m);
A(1,1) = 1;


for i=2:m-1
    vec = zeros(1,m);
    vec(i+1) = 1;
    for j=1:m
        A(i,j)=0.5*A(i-1,j);
    end
    A(i,i+1)=0.5;
end
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
