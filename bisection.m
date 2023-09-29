function midpoint = bisection(lower,upper,h)
format compact

for i=1:20
  midpoint=0.5*(lower+upper);
  disp([i midpoint])

  if shootingU(midpoint,h)*shootingU(upper,h)<0
     lower=midpoint;
  else
     upper=midpoint;
  end
end
