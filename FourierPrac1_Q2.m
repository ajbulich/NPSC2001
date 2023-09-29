clear

T = 2;
f0 = 1/T;
N = 24;

fDom = -24.999:0.001:25;
for f2=1:1:50000
    f = (f2-25000)/1000;
    E(f2) = (0.5*T*N*sin(N*pi*T*(f-f0)))/(pi*T*N*(f-f0)) - (0.5*T*N*sin(N*pi*T*(f+f0)))/(pi*N*T*(f+f0));
end

plot(fDom, E, 'b-')
xlabel('frequency (Hz)')
ylabel('Electric Field')
title(sprintf('T = %.0f seconds' , T))