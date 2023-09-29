clear

h=0.001;  
x =    [0 45.56015321 0 -100 0 400 -550];  % x coordinates only
xvel = [0 0 -214.9632443 0 135.7679907 -51.25066123 46.95366968]; % x velocities for all bodies
y =    [0 0 71.4988879 0 -163.7944498 357.6278992 -822.6261991]; % y coordinates for all bodies
yvel = [0 240.1366026 0 -181.0108715 0 57.32288935 -31.39277397]; % y velocities for all bodies
z =    [0 5.597312085 4.241555199 0 5.284807241 12.21375098 42.96253728]; % z coordinates for all bodies
zvel = [0 0 0 0 0 0 0]; % z velocities for all bodies
mArray = [3332217.013 0.552746149 8.149698593 10 1.069993302 3178.164769 951.6075017];
basis1 = [1 0 0];
basis2 = [0 1 0];
N = cross(basis1, basis2);

steps=15000;             % number of steps
skip=ceil(steps/300);   % display 300 frames
for i=1:steps
  % Fourth-order Runge-Kutta
  F1=deriv(x,y,z,xvel,yvel,zvel,mArray, 1);
  F2=deriv(x+h/2*F1(1),y+h/2*F1(2), z+h/2*F1(3), xvel+h/2*F1(4), yvel+h/2*F1(5), zvel+h/2*F1(6), mArray, 1);
  F3=deriv(x+h/2*F2(1),y+h/2*F2(2), z+h/2*F2(3), xvel+h/2*F2(4), yvel+h/2*F2(5), zvel+h/2*F2(6), mArray, 1);
  F4=deriv(x+h*F3(1),y+h*F3(2),z+h*F3(3),xvel+h*F3(4),yvel+h*F3(5),zvel+h*F3(6), mArray, 1);

  F1_1=deriv(x,y,z,xvel,yvel,zvel,mArray, 2);
  F2_1=deriv(x+h/2*F1_1(1),y+h/2*F1_1(2), z+h/2*F1_1(3), xvel+h/2*F1_1(4), yvel+h/2*F1_1(5), zvel+h/2*F1_1(6), mArray, 2);
  F3_1=deriv(x+h/2*F2_1(1),y+h/2*F2_1(2), z+h/2*F2_1(3), xvel+h/2*F2_1(4), yvel+h/2*F2_1(5), zvel+h/2*F2_1(6), mArray, 2);
  F4_1=deriv(x+h*F3_1(1),y+h*F3_1(2),z+h*F3_1(3),xvel+h*F3_1(4),yvel+h*F3_1(5),zvel+h*F3_1(6), mArray, 2);
  
  F1_2=deriv(x,y,z,xvel,yvel,zvel,mArray, 3);
  F2_2=deriv(x+h/2*F1_2(1),y+h/2*F1_2(2), z+h/2*F1_2(3), xvel+h/2*F1_2(4), yvel+h/2*F1_2(5), zvel+h/2*F1_2(6), mArray, 3);
  F3_2=deriv(x+h/2*F2_2(1),y+h/2*F2_2(2), z+h/2*F2_2(3), xvel+h/2*F2_2(4), yvel+h/2*F2_2(5), zvel+h/2*F2_2(6), mArray, 3);
  F4_2=deriv(x+h*F3_2(1),y+h*F3_2(2),z+h*F3_2(3),xvel+h*F3_2(4),yvel+h*F3_2(5),zvel+h*F3_2(6), mArray, 3);

  F1_3=deriv(x,y,z,xvel,yvel,zvel,mArray, 4);
  F2_3=deriv(x+h/2*F1_3(1),y+h/2*F1_3(2), z+h/2*F1_3(3), xvel+h/2*F1_3(4), yvel+h/2*F1_3(5), zvel+h/2*F1_3(6), mArray, 4);
  F3_3=deriv(x+h/2*F2_3(1),y+h/2*F2_3(2), z+h/2*F2_3(3), xvel+h/2*F2_3(4), yvel+h/2*F2_3(5), zvel+h/2*F2_3(6), mArray, 4);
  F4_3=deriv(x+h*F3_3(1),y+h*F3_3(2),z+h*F3_3(3),xvel+h*F3_3(4),yvel+h*F3_3(5),zvel+h*F3_3(6), mArray, 4);
    
  F1_4=deriv(x,y,z,xvel,yvel,zvel,mArray, 5);
  F2_4=deriv(x+h/2*F1_4(1),y+h/2*F1_4(2), z+h/2*F1_4(3), xvel+h/2*F1_4(4), yvel+h/2*F1_4(5), zvel+h/2*F1_4(6), mArray, 5);
  F3_4=deriv(x+h/2*F2_4(1),y+h/2*F2_4(2), z+h/2*F2_4(3), xvel+h/2*F2_4(4), yvel+h/2*F2_4(5), zvel+h/2*F2_4(6), mArray, 5);
  F4_4=deriv(x+h*F3_4(1),y+h*F3_4(2),z+h*F3_4(3),xvel+h*F3_4(4),yvel+h*F3_4(5),zvel+h*F3_4(6), mArray, 5);

  F1_5=deriv(x,y,z,xvel,yvel,zvel,mArray, 6);
  F2_5=deriv(x+h/2*F1_5(1),y+h/2*F1_5(2), z+h/2*F1_5(3), xvel+h/2*F1_5(4), yvel+h/2*F1_5(5), zvel+h/2*F1_5(6), mArray, 6);
  F3_5=deriv(x+h/2*F2_5(1),y+h/2*F2_5(2), z+h/2*F2_5(3), xvel+h/2*F2_5(4), yvel+h/2*F2_5(5), zvel+h/2*F2_5(6), mArray, 6);
  F4_5=deriv(x+h*F3_5(1),y+h*F3_5(2),z+h*F3_5(3),xvel+h*F3_5(4),yvel+h*F3_5(5),zvel+h*F3_5(6), mArray, 6);

  F1_6=deriv(x,y,z,xvel,yvel,zvel,mArray, 7);
  F2_6=deriv(x+h/2*F1_6(1),y+h/2*F1_6(2), z+h/2*F1_6(3), xvel+h/2*F1_6(4), yvel+h/2*F1_6(5), zvel+h/2*F1_6(6), mArray, 7);
  F3_6=deriv(x+h/2*F2_6(1),y+h/2*F2_6(2), z+h/2*F2_6(3), xvel+h/2*F2_6(4), yvel+h/2*F2_6(5), zvel+h/2*F2_6(6), mArray, 7);
  F4_6=deriv(x+h*F3_6(1),y+h*F3_6(2),z+h*F3_6(3),xvel+h*F3_6(4),yvel+h*F3_6(5),zvel+h*F3_6(6), mArray, 7);
 
 
   
  x(1) = x(1)+h/6*(F1(1)+2*F2(1)+2*F3(1)+F4(1));
  y(1) = y(1)+h/6*(F1(2)+2*F2(2)+2*F3(2)+F4(2));
  z(1) = z(1)+h/6*(F1(3)+2*F2(3)+2*F3(3)+F4(3));
  xvel(1) = xvel(1)+h/6*(F1(4)+2*F2(4)+2*F3(4)+F4(4));
  yvel(1) = yvel(1)+h/6*(F1(5)+2*F2(5)+2*F3(5)+F4(5));
  zvel(1) = zvel(1)+h/6*(F1(6)+2*F2(6)+2*F3(6)+F4(6));

  x(2) = x(2)+h/6*(F1_1(1)+2*F2_1(1)+2*F3_1(1)+F4_1(1));
  y(2) = y(2)+h/6*(F1_1(2)+2*F2_1(2)+2*F3_1(2)+F4_1(2));
  z(2) = z(2)+h/6*(F1_1(3)+2*F2_1(3)+2*F3_1(3)+F4_1(3));
  xvel(2) = xvel(2)+h/6*(F1_1(4)+2*F2_1(4)+2*F3_1(4)+F4_1(4));
  yvel(2) = yvel(2)+h/6*(F1_1(5)+2*F2_1(5)+2*F3_1(5)+F4_1(5));
  zvel(2) = zvel(2)+h/6*(F1_1(6)+2*F2_1(6)+2*F3_1(6)+F4_1(6));
  

  x(3) = x(3)+h/6*(F1_2(1)+2*F2_2(1)+2*F3_2(1)+F4_2(1));
  y(3) = y(3)+h/6*(F1_2(2)+2*F2_2(2)+2*F3_2(2)+F4_2(2));
  z(3) = z(3)+h/6*(F1_2(3)+2*F2_2(3)+2*F3_2(3)+F4_2(3));
  xvel(3) = xvel(3)+h/6*(F1_2(4)+2*F2_2(4)+2*F3_2(4)+F4_2(4));
  yvel(3) = yvel(3)+h/6*(F1_2(5)+2*F2_2(5)+2*F3_2(5)+F4_2(5));
  zvel(3) = zvel(3)+h/6*(F1_2(6)+2*F2_2(6)+2*F3_2(6)+F4_2(6));

  x(4) = x(4)+h/6*(F1_3(1)+2*F2_3(1)+2*F3_3(1)+F4_3(1));
  y(4) = y(4)+h/6*(F1_3(2)+2*F2_3(2)+2*F3_3(2)+F4_3(2));
  z(4) = z(4)+h/6*(F1_3(3)+2*F2_3(3)+2*F3_3(3)+F4_3(3));
  xvel(4) = xvel(4)+h/6*(F1_3(4)+2*F2_3(4)+2*F3_3(4)+F4_3(4));
  yvel(4) = yvel(4)+h/6*(F1_3(5)+2*F2_3(5)+2*F3_3(5)+F4_3(5));
  zvel(4) = zvel(4)+h/6*(F1_3(6)+2*F2_3(6)+2*F3_3(6)+F4_3(6));

  x(5) = x(5)+h/6*(F1_4(1)+2*F2_4(1)+2*F3_4(1)+F4_4(1));
  y(5) = y(5)+h/6*(F1_4(2)+2*F2_4(2)+2*F3_4(2)+F4_4(2));
  z(5) = z(5)+h/6*(F1_4(3)+2*F2_4(3)+2*F3_4(3)+F4_4(3));
  xvel(5) = xvel(5)+h/6*(F1_4(4)+2*F2_4(4)+2*F3_4(4)+F4_4(4));
  yvel(5) = yvel(5)+h/6*(F1_4(5)+2*F2_4(5)+2*F3_4(5)+F4_4(5));
  zvel(5) = zvel(5)+h/6*(F1_4(6)+2*F2_4(6)+2*F3_4(6)+F4_4(6));
  
  x(6) = x(6)+h/6*(F1_5(1)+2*F2_5(1)+2*F3_5(1)+F4_5(1));
  y(6) = y(6)+h/6*(F1_5(2)+2*F2_5(2)+2*F3_5(2)+F4_5(2));
  z(6) = z(6)+h/6*(F1_5(3)+2*F2_5(3)+2*F3_5(3)+F4_5(3));
  xvel(6) = xvel(6)+h/6*(F1_5(4)+2*F2_5(4)+2*F3_5(4)+F4_5(4));
  yvel(6) = yvel(6)+h/6*(F1_5(5)+2*F2_5(5)+2*F3_5(5)+F4_5(5));
  zvel(6) = zvel(6)+h/6*(F1_5(6)+2*F2_5(6)+2*F3_5(6)+F4_5(6));

  x(7) = x(7)+h/6*(F1_6(1)+2*F2_6(1)+2*F3_6(1)+F4_6(1));
  y(7) = y(7)+h/6*(F1_6(2)+2*F2_6(2)+2*F3_6(2)+F4_6(2));
  z(7) = z(7)+h/6*(F1_6(3)+2*F2_6(3)+2*F3_6(3)+F4_6(3));
  xvel(7) = xvel(7)+h/6*(F1_6(4)+2*F2_6(4)+2*F3_6(4)+F4_6(4));
  yvel(7) = yvel(7)+h/6*(F1_6(5)+2*F2_6(5)+2*F3_6(5)+F4_6(5));
  zvel(7) = zvel(7)+h/6*(F1_6(6)+2*F2_6(6)+2*F3_6(6)+F4_6(6));
  
  % Save and plot trajectory
  X(i)=x(1); 
  Y(i)=y(1);
  Z(i) = z(1);
  X2(i) = x(2);
  Y2(i) = y(2);
  Z2(i) = z(2);
  X3(i) = x(3);
  Y3(i) = y(3);
  Z3(i) = z(3);
  X4(i) = x(4);
  Y4(i) = y(4);
  Z4(i) = z(4);
  X5(i) = x(5);
  Y5(i) = y(5);
  Z5(i) = z(5);
  X6(i) = x(6);
  Y6(i) = y(6);
  Z6(i) = z(6);
  X7(i) = x(7);
  Y7(i) = y(7);
  Z7(i) = z(7);
 
  comStats = coM(x,y,z,xvel, yvel, zvel, mArray);
  com = comStats(1:3);
  vcom = comStats(4:6);
  L2(i) = angMom(x,y,z,xvel,yvel,zvel,mArray,2,com,vcom);
  vel2 = [xvel(2) yvel(2) zvel(2)];
  angle2(i) = rad2deg(pi/2 - acos( dot(vel2, N)/norm(N)/norm(vel2) ));
  e0 = 0;
  l0 = 0;
  for ii=1:length(mArray)
    e0 = e0 + totalE(x,y,z,xvel,yvel,zvel, mArray, ii);
    l0 = l0 + angMom(x,y,z,xvel,yvel,zvel,mArray,ii,com,vcom);
  end
  E1(i) = e0;
  L3(i) = l0;
   
  if rem(i,skip)==0
    clf
    hold on
    view([-1 -1 0.5])
    plot3(X,Y,Z,'g-', X(end),Y(end),Z(end),'go')
    plot3(X2, Y2, Z2, 'b-', X2(end), Y2(end), Z2(end), 'bo')
    plot3(X3, Y3, Z3, 'r-', X3(end), Y3(end), Z3(end), 'ro')
    plot3(X4, Y4, Z4, 'k-', X4(end), Y4(end), Z4(end), 'ko')
    plot3(X5, Y5, Z5, 'k-', X5(end), Y5(end), Z5(end), 'ko')
    plot3(X6, Y6, Z6, 'k-', X6(end), Y6(end), Z6(end), 'ko')
    plot3(X7, Y7, Z7, 'k-', X7(end), Y7(end), Z7(end), 'ko')
    plot3(0, 0, 0, 'mo')
    hold off
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    drawnow
  end
end


figure(2)
plot(1:steps, L2, 'r-')
ylabel("L (relative units)", 'FontSize', 20)
xlabel("Timesteps", 'FontSize', 20)
title("Angular Momentum of Mercury", 'FontSize',24)

figure(3)
plot(1:steps, angle2, 'b-')
ylabel("Angle (degrees)",'FontSize', 20)
xlabel("Timesteps", 'FontSize', 20)
title("Angle of inclination of Mercury's velocity",'FontSize', 24)

figure(4)
plot(1:steps, E1, 'm-')
ylabel("Energy (relative units)",'FontSize', 20)
xlabel("Timesteps", 'FontSize', 20)
title("Total Energy of System",'FontSize', 24)

figure(5)
plot(1:steps, L3, 'r-')
ylabel("L (relative units)", 'FontSize', 20)
xlabel("Timesteps", 'FontSize', 20)
title("Angular Momentum of System", 'FontSize',24)

  
% figure(2)
% plot(1:steps, K, 'b-')
% 
% figure(3)
% plot(1:steps, L, 'r-')
% ylabel("Angular Momentum")

% -------------------

function value = deriv(x,y,z,xvel, yvel, zvel, mArray, massNo)


Cartesian_pos=[x(massNo) y(massNo) z(massNo)];
vel=[xvel(massNo), yvel(massNo), zvel(massNo)];
accel = 0;
for j = 1:length(mArray)
    if j ~= massNo
        Cartesian_pos2 = [x(j) y(j) z(j)];
        pos = Cartesian_pos - Cartesian_pos2;
        r=norm(pos);
        accel = accel + (-1*mArray(j)/(r^2) * pos/r);
    end
end
value(1:3)=vel;
value(4:6)=accel;

end

% --------------------
function value3 = coM(x,y,z,xvel,yvel,zvel,mArray)

M = sum(mArray);
sum_mr = 0;
sum_vr = 0;
for i=1:length(mArray)
    cartesian_pos = [x(i) y(i) z(i)];
    cartesian_vel = [xvel(i) yvel(i) zvel(i)];
    sum_mr = sum_mr + mArray(i)*cartesian_pos;
    sum_vr = sum_vr + mArray(i)*cartesian_vel;
end

com = sum_mr / M;
vcom = sum_vr/M;

value3(1:3) = com;
value3(4:6) = vcom;

end

% --------------------
function value2 = angMom(x,y,z,xvel,yvel,zvel,mArray,massNo,com,vcom)

Cartesian_pos=[x(massNo) y(massNo) z(massNo)];
Cartesian_vel=[xvel(massNo), yvel(massNo), zvel(massNo)];
M = sum(mArray);

Cartesian_pos2 = com;
Cartesian_vel2 = vcom;
vel = Cartesian_vel - Cartesian_vel2;
pos = Cartesian_pos - Cartesian_pos2;
L = mArray(massNo)*norm(cross(pos,vel));

value2 = L;

end

%------------------

function value4 = totalE(x,y,z,xvel, yvel, zvel, mArray, massNo)


Cartesian_pos=[x(massNo) y(massNo) z(massNo)];
vel=[xvel(massNo) yvel(massNo) zvel(massNo)];
E = 0;
for j = 1:length(mArray)
    if j ~= massNo
        Cartesian_pos2 = [x(j) y(j) z(j)];
        pos = Cartesian_pos - Cartesian_pos2;
        r = norm(pos);
        v = norm(vel);
        E = E + (0.5*mArray(j)*v^2 - mArray(massNo)*mArray(j)/r);
    end
end
value4 = E;
end


