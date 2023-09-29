clear
format compact

[e1,X1,PSI1] = originalShooting(0.5);
[e2,X2,PSI2] = originalShooting(1.5);
[e3,X3,PSI3] = originalShooting(2.5);

plot(X1,PSI1,'r-',X2,PSI2,'b-',X3,PSI3,'g-')
legend('PSI1', 'PSI2', 'PSI3')

xlength = max(X1) - min(X1);
h = length(PSI1)/xlength;

func = PSI2.*PSI2;
figure(2)
plot(X1, func, 'r-')

for i=1:800
    A1(i) = 0.01*func(i);
end

A = sum(A1)

