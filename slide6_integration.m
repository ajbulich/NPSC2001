clear; clf

x=-4;
h=0.01;
steps=800;

psi=1;
prev=0;
e=0.8196598;   % ground-state solution

for i=1:steps
    X(i)=x;
    PSI(i)=psi;

    % Square-well potential with m=hbar=1
    if abs(x)<1
        pot=0;
    else
        pot=10;
    end
    d2psi=2*(pot-e)*psi;
    POT(i)=pot;
    
    % Verlet Method
    next=2*psi - prev + h*h*d2psi;

    prev=psi;
    psi=next;
    x=x+h;
end

% Compute density to normalise wavefunction 
rho=PSI.^2;
area=h*sum(rho);
PSI=PSI/sqrt(area);

subplot(2,1,1)
plot(X,PSI,'ro')
subplot(2,1,2)
plot(X,POT,'b-');
ylim([-1 11])
