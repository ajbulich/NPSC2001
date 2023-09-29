clear
format compact

hList = [0.05 0.1 0.2 0.3 0.4];

correctE = 0.5;

for i=1:5
    h = hList(i);
    ans(i) = bisection(0.4,0.6,h);
    correct(i) = 0.5;
end

error = correct-ans;
disp(error)
plot(hList, error)




