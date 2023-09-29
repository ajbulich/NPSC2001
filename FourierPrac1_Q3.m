clear

a = 0.2;
d = 1;
N = 20;

uDom = -24.99:0.01:25;
for u2=1:1:5000
    u = (u2-2500)/100;
    sum = 0;
    for n=1:N
        sum = sum + exp(2*pi*1i*u*(n-1)*(a+d));
    end
    num = abs(sum);
    E(u2) = a*(sin(pi*u*a))/(pi*u*a) .* num ;
end

size(uDom)
size(E)
plot(uDom, E, 'b-')
title(sprintf('N = %.0f', N))
xlabel('u (m^-1)')
ylabel("Electric Field")