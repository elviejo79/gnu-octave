t=0:50;
a=randn(size(t));
original = a;
a(3:end)=a(3:end)+.2*a(2:end-1)-.4*a(1:end-2);

b1=rand;
for i=2:length(a)
  b1(i) = 1.05*b1(i-1) +randn/3;
end

clf
subplot(321),plot(t,original)
title("original series")

subplot(322),plot(t,original)
title("original series")

subplot(323),plot(t,a)
title("Stationary")

subplot(324),plot(t,b1)
title("Non-stationary")

## Detrending is a technique to make a non-stationary time series more stationary.
subplot(325), plot(t,detrend(b1))
title("Detrended")

## pre-withening: take the derivatie( the value at each time point is re-assigned to be the difference between that and the value at the previous time point).
b1deriv=diff(b1);

## Because the derivate will reduce the number of time points by one the time series can be padded at the end.
b1deriv(end+1)=b1deriv(end);
subplot(326), plot(t,b1deriv)
title("Derivative")
