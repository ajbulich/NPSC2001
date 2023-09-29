clear;

% Parameters for the system
c=1;
tauVec=linspace(0.004,0.0099);
h=0.01;
tmax = abs(h/c);
    % pulse loops through system once

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



% Lax-Wendroff update matrix
% mat=eye(m) + c*tau/(2*h)*D + 0.5*(c*tau/h)^2*D2;

% Initial conditions (Gaussian pulse at x=0.5)
% wave=exp(-100*(x-0.5).^2)';


% Main loop

for jj=1:100
    % Lax update matrix
    tau1 = tauVec(jj);
    clear mat
    mat=av + c*tau1/(2*h)*D;
    mat(:, 1) = 0;
    mat(:, end) = 0;

    steps=ceil(1/(c*tau1)); 
    clear wave
    for ii=1:length(x)
        if (0.5 <= x(ii)) && (x(ii) < 0.75)
            wave(ii) = 1;
        else
            wave(ii) = 0;
        end
    end
    wave = wave';
    wave0=wave; 
    for n=1:0.25*steps
        wave=mat*wave;    % Update current values
    end
    width(jj) = sum(wave>0.001)/m-0.25;
    figure(1)
    plot(x,wave0,'ro-',x,wave,'bo-')
    drawnow
    
end
t_tmax = 1/tmax * tauVec;
figure(2)
plot(t_tmax, width, 'ro-')
ylabel("Width")
xlabel('tau ratio')



