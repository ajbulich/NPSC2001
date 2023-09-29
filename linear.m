clear
global phi A b N

% Deine grid
h=0.01;
x=0:h:1;
y=0:h:1;

% Determine matrix dimensions
m=length(x);
N=m-2;

% Boundary conditions
phi=ones(m);
phi(m,:)=2;

% Construct matrix A and vector b
A=-4*eye(N^2);
b=zeros(N^2,1);
tic
for i=1:N
  for j=1:N
    row=(i-1)*N+j;

    add(row,i-1,j);
    add(row,i+1,j);
    add(row,i,j-1);
    add(row,i,j+1);
  end
end


% Solve for interior points
X=A\b;
toc

% Combine boundary conditions and interior points
phi(2:m-1,2:m-1)=reshape(X,N,N);

% Plot result
mesh(x,y,phi)
disp(N)
