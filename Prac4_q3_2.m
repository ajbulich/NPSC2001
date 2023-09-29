clear
format compact

[e1,X1,PSI1] = shootingU(-0.482120656967163, 0.01);
[e2,X2,PSI2] = shootingU(-0.122732982635498, 0.01);
[e3,X3,PSI3] = shootingU(-0.054733667373657, 0.01);

plot(X1,PSI1,'r-',X2,PSI2,'b-',X3,PSI3,'g-')
legend('PSI1', 'PSI2', 'PSI3')


eList = [-0.482120656967163,-0.122732982635498,-0.054733667373657, -0.024078044891357];
yList = [1,2,3,4];

figure(2)
plot(eList,log(yList))
