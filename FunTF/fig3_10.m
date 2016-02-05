t=0:.01:11;
clf
plot(t, abs(mod(t,2)-1.25))
  hold on
  plot(t,abs(mod(t,2)-1),'r')
