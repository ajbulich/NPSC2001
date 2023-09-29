clear

c = 1;
tau = 0.005;
h = 0.01;

x = h:h:1; %had to get rid of vector zero
m = length(x);

% wavefunction moves by some distance c*tau, this is one loop through the
% system
steps = round(1/abs(c*tau));

D =     diag(ones(m-1,1),+1);
D = D - diag(ones(m-1,1),-1);
D(m,1) = 1;  %modification for boundary conditions
D(1,m) = -1;

%update matrix:
fac = c*tau/(2*h);
%A = (eye(m) + fac*D);

% For lax method only:
%AvgMat = 0.5*abs(D);
%A = (AvgMat + fac*D);

% For Lax-Wendroff
D2 = abs(D)-2*eye(m);
newMat = c^2*tau^2/(2*h^2) * D2;
A = (eye(m) + fac*D) + newMat;


wave = exp(-100*(x-0.5).^2)'; %transposed so that we can do the matrix multiplication

for n=1:steps
    wave = A*wave;
    plot(x,wave, 'ro-')
    title(n*tau)
    drawnow
end

