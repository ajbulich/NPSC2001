clear
x1 = 0:0.01:8;
h=0.1;
alpha = 0.0005;

L = 8;
c = 1;
w1 = pi*c/L;
T = 2*pi/w1;

Nterms = 10;
j = 1;
for t=0:h:8
    for x=0:0.01:8
        v = 0;
        u(j) = 0;
        for m=0:Nterms
            a = (-1)^m * (32)/(pi^2*(2*m+1)^2)*(1-cos(pi*(2*m+1)/(8)))*sin(pi*x*(2*m+1)/(L))*cos(sqrt(1+alpha*((pi*(2*m+1))/(L))^2)*pi*t*(2*m+1)/(L));
            v = v + a;
        end
        u(j) = v;
        j = j+1;
    end
    plot(x1,u)
    pause(0.1)
    j = 1;
end

