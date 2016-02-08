## "White noise refers to noise that has a flat power spectrum
## % the functions rand and randn produce data that can be considered white noise.
1;
x=[1:1000];
Yu=rand(1000,1);
length(Yu)
Yn=randn(1000,1);

subplot(221)
plot(Yn), hold on
plot(Yu,'r')
title('Random noise over time')

subplot(223), hist(Yu,200)
title('Distribution of uniform noise')

subplot(224), hist(Yn,200)
title('Distribution of random noise')

% Pink noise refers to noise with a non-uniform frequency
% typically, that the power decreases with increasing frequency
% One way to compute pink noise; is to apply a vanishing frequency filter.

% wn = white noise
wn = randn(1000,1);
wnX = fft(wn);

pn = real(ifft(wnX .* linspace(-1,1,length(wnX))'.^2))*2;

subplot(221)
## plot(wn), hold on
## plot(pn,'r')
## xlabel('Time (a.u.)')
## ylabel('Amplitude (a.u.)')
## legend({'white', 'pink'})

## subplot(222)
## plot(wn,pn,'.')
## xlabel('Amplitude white noise')
## ylabel('Amplitude pink noise')

## subplot(212)
## plot(abs(fft(wn))), hold on
## plot(abs(fft(pn)), 'r')
## legend({'white';'pink'})
## xlabel('Frequency (a.u.)'),
## ylabel('Apmlittude')
