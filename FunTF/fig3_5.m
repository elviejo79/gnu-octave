t = [0:.001:5];
n = length(t)
a = [10 2 5 8];
f = [3 1 6 12];
tchunks = round(linspace(1,n,length(a)+1))
swave=0;
for i=1:length(a)
        swave=cat(2,swave,a(1)*sin(2*pi*f(i)*t(tchunks(i):tchunks(i+1)-1)));
end
clf
plot(t,swave)
