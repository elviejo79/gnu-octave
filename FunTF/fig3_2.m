## // WHITE NOISE has a flat power spectrum
## // PINK NOISE noise with a non-uniform frequency structure.

white_noise = randn(1000,1);
wnX = fft(white_noise); ##fast fourier transfrom
pink_noise = real(ifft(wnX .* linspace(-1,1, length(wnX))'.^2))*2;

clf;
subplot(221)
plot(white_noise), hold on
plot(pink_noise, 'r')
xlabel('Time (au)')
ylabel('Amplitude (au)')
legend({'white', 'pink'})

subplot(222)
plot(white_noise, pink_noise,'.')
xlabel('Amplitude white noise')
ylabel('Amplitude pinkified noise')

subplot(212)
plot(abs(fft(white_noise))), hold on
plot(abs(fft(pink_noise)),'r')
legend({'white';'pink'})
xlabel('Frequency (au)'),
ylabel('Amplitude')

