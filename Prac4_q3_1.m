clear
format compact

h = 0.001;
e = -1;
i = 1;
while e <= 0
   E(i) = shootingU(e,h);
   xvals(i) = e;
   i = i + 1;
   e = e + h;
end

plot(xvals, E, 'b-')
disp(max(E))
yline(0)




