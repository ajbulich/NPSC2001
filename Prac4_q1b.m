clear
format compact

h = 0.05;
e = 0;
i = 1;
while e <= 4
   E(i) = shooting(e,h);
   xvals(i) = e;
   i = i + 1;
   e = e + 0.05;
end

plot(xvals, E, 'b-')




