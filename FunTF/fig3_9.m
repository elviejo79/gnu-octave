t=0:.01:11;
clf
plot(t,mod(t,2)>1), hold on
  plot(t,.02+(mod(.25+t,2)>1.5),'r')
