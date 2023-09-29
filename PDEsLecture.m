clear

kappa = 1;
tau = 1e-4;
h = 0.02;
steps = 20;

x = 0:h:1;
m = length(x);

D = -2*eye(m) + diag(ones(m-1,1), 1) + diag(ones(m-1,1), -1);
D = D * kappa * tau / h^2;

% Dirichlet boundary conditions
D(1,:) = 0;
D(m,:) = 0;

% Initial conditions
temp = zeros(m,1); % m rows, 1 column
temp(round(m/2)) = 1/h;
temp(1) = 10;

time = 0;

for n=1:steps
    plot(x, temp, 'ro-')
    title(time)
    drawnow

    time = time + tau;
    temp = temp + D*temp;
    middle(n) = abs(temp(round(m/2)));
end

A = eye(m) + D;


