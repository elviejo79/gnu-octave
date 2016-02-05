t = [0:.001:5];
n = length(t)
f = [2 10];                     #frequencies in Hz
ff = linspace(f(1),f(2)*mean(f)/f(2),n); # doesn't it f(2)*mean(f)/f(2) return the exct same as mean(f)??

swave = sin(2*pi.*ff.*t);
clf
plot(t,swave)
