clear

global g L
g = 9.8;
L = 0.2482;

theta = 30 * pi/180;
x = [theta 0];
h = 0.02;

for i=1:400
    % RK4
    F1 = deriv(x);
    F2 = deriv(x+0.5*h*F1);
    F3 = deriv(x+0.5*h*F2);
    F4 = deriv(x+h*F3);
    x = x + h/6 * (F1+2*F2+2*F3+F4);

    % Graphics of pendulum
    xpend = [0 L*sin(x(1))];
    ypend = [0 -L*cos(x(1))];

    plot(xpend,ypend, 'o-')
    axis equal
    axis([-L L -L 0]);
    drawnow
end