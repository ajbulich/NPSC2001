clear

h=0.001;  
x=[1 0 1 0 0.5 0];  % initial conditions: [x y z vx vy vz]
x2 = [-1 0 0 0 -0.5 0];% initial conditions of second object
x3 = [-3 -3 1 0 -1 -0.5];
m = 500;
m2 = 500;
m3 = 5;

steps=80000;             % number of steps
skip=ceil(steps/300);   % display 200 frames
for i=1:steps
  % Fourth-order Runge-Kutta
  F1=deriv(x,m, x2, m2);
  F2=deriv(x+h/2*F1,m, x2, m2);
  F3=deriv(x+h/2*F2,m, x2, m2);
  F4=deriv(x+h*F3,m, x2, m2);
    
  F1_2=deriv(x2,m2, x, m);
  F2_2=deriv(x2+h/2*F1_2,m2, x, m);
  F3_2=deriv(x2+h/2*F2_2,m2, x, m);
  F4_2=deriv(x2+h*F3_2,m2, x, m);
  
  F1_3=deriv(x3,m3, x, m);
  F2_3=deriv(x3+h/2*F1_3,m3, x, m);
  F3_3=deriv(x3+h/2*F2_3,m3, x, m);
  F4_3=deriv(x3+h*F3_3,m3, x, m);

  x3 = x3+h/6*(F1_3+2*F2_3+2*F3_3+F4_3);

  F1_4=deriv(x3,m3, x2, m2);
  F2_4=deriv(x3+h/2*F1_4,m3, x2, m2);
  F3_4=deriv(x3+h/2*F2_4,m3, x2, m2);
  F4_4=deriv(x3+h*F3_4,m3, x2, m2);

  x=x+h/6*(F1+2*F2+2*F3+F4);
  x2=x2+h/6*(F1_2+2*F2_2+2*F3_2+F4_2);
  x3 = x3+h/6*(F1_4+2*F2_4+2*F3_4+F4_4);

  F1_5=deriv(x2,m2, x3, m3);
  F2_5=deriv(x2+h/2*F1_5,m2, x3, m3);
  F3_5=deriv(x2+h/2*F2_5,m2, x3, m3);
  F4_5=deriv(x2+h*F3_5,m2, x3, m3);

  F1_6=deriv(x,m, x3, m3);
  F2_6=deriv(x+h/2*F1_6,m, x3, m3);
  F3_6=deriv(x+h/2*F2_6,m, x3, m3);
  F4_6=deriv(x+h*F3_6,m, x3, m3);
 
    
  % Save and plot trajectory
  X(i)=x(1); 
  Y(i)=x(2);
  Z(i) = x(3);
  X2(i) = x2(1);
  Y2(i) = x2(2);
  Z2(i) = x2(3);
  X3(i) = x3(1);
  Y3(i) = x3(2);
  Z3(i) = x3(3);
  vel = transpose([x(3), x(4), 0]);
  pos = transpose([x(1), x(2), 0]);
  K(i) = 0.5*norm(vel)^2;
  L(i) = norm(cross(pos,vel));
   
  if rem(i,skip)==0
    plot3(X,Y,Z,'g-', X2, Y2, Z2, 'b-',X3, Y3, Z3, 'r-', X(end),Y(end),Z(end),'ko', X2(end), Y2(end), Z2(end), 'ko',X3(end), Y3(end), Z3(end), 'ko',0,0,0,'ro')
    axis equal
    drawnow
  end
end

figure(2)
plot(1:steps, K, 'b-')

figure(3)
plot(1:steps, L, 'r-')
ylabel("Angular Momentum")

% -------------------

function value = deriv(vec, m, vec2, m2)

k=0.0000;
w=0.00;
pos=vec(1:3) - vec2(1:3);
vel=vec(4:6);

r=norm(pos);
accel=-m2/(m*r^2) * pos/r + w*(norm(cross(vel, pos))^2 / r^3 ) * pos/r - (k*(norm(cross(vel, pos))^2 / r^4 )) * pos/r;

value(1:3)=vel;
value(4:6)=accel;

end
