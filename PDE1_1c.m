clear

h=0.05;    % spatial-step

% Extent of the system
x=0:h:1;
m=length(x);

% Matrix for the second-derivative term
D=-2*eye(m);
D=D+diag(ones(m-1,1),+1);
D=D+diag(ones(m-1,1),-1);

% Boundary conditions (Dirichlet; i.e. constant value)
D(1,:)=0;
D(m,:)=0;

% The update matrix
i = 0;
for fac=0:0.05:1
    A=eye(m) + D*fac;
    i = i+1;

    x(i) = fac;
    emin(i) = min(eig(A));
    emax(i) = max(eig(A));
    eabs(i) = max(abs(eig(A)));
end

plot(x,emin, 'bo-', x, emax, 'ro-', x, eabs, 'yo-')
legend('Minimum', 'Maximum', 'Dominant')
xlabel('\kappa\tau/h^2')
ylabel('Eigenvalue')