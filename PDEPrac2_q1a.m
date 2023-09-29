clear;

% Parameters for the system
c=1;
h=0.01;
tmax = abs(h/c);
tau = 0.5*tmax;
steps=ceil(1/(c*tau));     % pulse loops through system once

% Extent of the system
x=h:h:1;       % no mesh point at x=0 due to periodic boundary conditions
m=length(x);

% Diagonal matrix to implement the spatial first-derivative
D=diag(ones(m-1,1),+1);
D=D-diag(ones(m-1,1),-1);

% Additional elements for periodic boundary conditions
D(m,1)=1;
D(1,m)=-1;

% Construct the spatial-averaging matrix
av=0.5*abs(D);

% Construct the matrix for the spatial second-derivative
D2=abs(D)-2*eye(m);

% FTCS update matrix
%mat=eye(m) + c*tau/(2*h)*D;

% Lax update matrix
 mat=av + c*tau/(2*h)*D;
 mat(:, 1) = 0;
 mat(:, end) = 0;

% Lax-Wendroff update matrix
% mat=eye(m) + c*tau/(2*h)*D + 0.5*(c*tau/h)^2*D2;

% Initial conditions (Gaussian pulse at x=0.5)
% wave=exp(-100*(x-0.5).^2)';
for ii=1:length(x)
    if (0.5 <= x(ii)) && (x(ii) < 0.75)
        wave(ii) = 1;
    else
        wave(ii) = 0;
    end
end
wave = wave';
wave0=wave; 


% Main loop
for n=1:0.25*steps
    wave=mat*wave;    % Update current values

    plot(x,wave0,'ro-',x,wave,'bo-')
    title(n*tau)        
    drawnow
end

width = sum(wave>0.001)/m-0.25;
disp(width)

