clear
format short
% extent of system
h = 0.1;
x = 0:h:1;
m = length(x);


omega = linspace(1,1.9, 91);
ii = 1;
for w=omega
    % creating matrices
    A = zeros(m,m);
    A(1,1) = 1;
    
    for i=2:m-1
        for j=1:m
            A(i,j)=w/2*A(i-1,j);
        end
        A(i,i) = A(i,i) + 1-w;
        A(i,i+1) = A(i,i+1) + w/2;
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
    iter(ii) = n;
    B = sort(abs(eig(A)));
    domEig(ii) = B(end-2);
    ii = ii + 1;
    
end
plot(omega, iter, 'ro-')
xlabel('omega')
ylabel('iterations')

figure(2)
plot(omega, domEig, 'go-')
xlabel('omega')
ylabel('Magnitude of Dominant Eigenvalue')






