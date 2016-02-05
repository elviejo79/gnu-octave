t = [0:.001:5];
n = length(t);
a = linspace(1,10,n);           #time-varying amplitude
f = 3;                     #frequencies in Hz

swave = a.*sin(2*pi.*f.*t);
clf
plot(t,swave)
