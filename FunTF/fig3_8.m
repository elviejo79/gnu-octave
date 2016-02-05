t=-1:.001:5;
s=[.5 .1]; #widths
a=[4 5]; #amplitudes
## b=[0 2];
## gaussian = @(a,b,s) a.*exp((-(t-b).^2)/(2.*s(2)^2));
           
g1 = a(1)*exp(-(t).^2/(2*s(1)^2));
g2 = a(2)*exp((-(t-2).^2)/(2*s(2)^2));
clf
plot(t,g1), hold on, plot(t,g2,"r")
