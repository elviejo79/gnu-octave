t=0:.01:11;
clf;
subplot(221), plot(t, sin(exp(t-5)))
subplot(222), plot(t, log(t)./(t.^2))
set(gca,'ylim',[0 .2])
subplot(223), plot(t, sin(t).*exp((-(t-3).^2)))
subplot(224), plot(t, abs(mod(t,2)-.66).*sin(2*pi*10*t))
