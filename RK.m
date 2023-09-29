clear
format compact

x=1;
h=1;


% Euler Method

%F1 = deriv(x);
%x = x + h*F1;

% RK2 method
%F1 = deriv(x);
%F2 = deriv(x+h*F1);
%x = x + 0.5*h*(F1+F2);

% Rk4
F1 = deriv(x);
F2 = deriv(x+0.5*h*F1);
F3 = deriv(x+0.5*h*F2);
F4 = deriv(x+h*F3);

x = x + h/6 * (F1+2*F2+2*F3+F4);