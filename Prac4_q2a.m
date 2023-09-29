clear
format compact

[e1,X1,PSI1] = originalShooting(0.5);
[e2,X2,PSI2] = originalShooting(1.5);
[e3,X3,PSI3] = originalShooting(2.5);

plot(X1,PSI1,'r-',X2,PSI2,'b-',X3,PSI3,'g-')
legend('PSI1', 'PSI2', 'PSI3')

PSI4 = 0.8694*(2*X3.^2 - 1).*PSI1;

%plot(PSI4,PSI3, 'g-')