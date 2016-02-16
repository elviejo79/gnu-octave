1;

function fig7_1()
  n=20;
  d=randn(n,1);
  dx=fft(d);
                                #positive frequencies
  posF = 2:floor(n/2);
                                #negative requencies
  negF = floor(n/2)+2:n;
  dx(posF)=dx(posF)*2;
  dx(negF)=0;
  hilbertd = ifft(dx);

  plot(d), hold on
  plot(real(hilbertd), 'ro')
  title("Fig 7.1| the real part of the result of the hilbert transform is the original time series")
endfunction

function hilbertd=hilbert(signal)
  dX = fft(signal);
  dX(dX>0)= dX(dX>0).* 2;
  dX(dX<0)= dX(dX<0) .* 0;
  hilbertd = ifft(dX);
endfunction

function fig7_2()
  t = 0:.01:5;
  signal = zeros(size(t));
  a = [10 2 5 8];
  f = [3 1 6 12];
  for i=1:length(a)
    signal = signal + a(i)*sin(2*pi*f(i)*t);
  end

  hilsine = hilbert(signal);
  subplot(311), plot(t,signal)
  ylabel("Amplitude");
  
  subplot(312), plot(t,abs(hilsine).^2)
  ylabel("power");
  
  subplot(313), plot(t,angle(hilsine))
  ylabel("phase");
endfunction

function fig7_3()
  t = 0:.001:1;
  sinewave = 3*sin(2*pi*5*t);
  plot(t, sinewave), hold on
  sinewaveDamp = real(ifft(fft(sinewave)*.5));
  plot(t, sinewaveDamp, 'r');
  title("Fig 7.3| A 3 amplitude sine wave is attenuated after scaling the Fourier coefficients by 50%")
endfunction

function fig7_4()
           # in the code below, a sum of two sine waves (5hz and 10hz)
           # is created, and the 10hz is attenuated by down-scaling
           # the corresponding frequencies in the frequency domaind
           # and then computing the inverse fourier transform
           #to get back to the time domain

  srate=1000;
  t = 0:1/srate:1;
  n = length(t);
  signal = 3*sin(2*pi*5*t) + 4*sin(2*pi*10*t);
  subplot(211),
  plot(t, signal),
  hold on

  x = fft(signal)/n;
  hz = linspace(0, srate/2, floor(n/2)+1);

  subplot(212)
  plot(hz, 2*abs(x(1:length(hz))), '-*'), hold on

                                # frequencies to attenuate
  hzidx = dsearchn(hz', [8 14]');
  x(hzidx(1):hzidx(2)) = .1*x(hzidx(1):hzidx(2));

  subplot(211)
  plot(t, real(ifft(x)*n), 'r')

  subplot(212)
  plot(hz, 2*abs(x(1:length(hz))), 'r-o')
  xlim([0 15])
  title("Fig 7.4|Frequency selective dampening of a multi-sine signal.")
endfunction

function fig7_5()
  srate = 1000;
  t = 0:1/srate:5;
  sinewave = sin(2*pi*2*t);
  subplot(211), plot(t,sinewave), hold on

  x = fft(sinewave)/length(t);
  hz = linspace(0, srate/2, floor(length(t)/2+1));

  clf
  subplot(212)
  plot(hz, 2*abs(x(1:length(hz))), '-*'), hold on

  hzidx = dsearchn(hz', [8 9]');
  x(hzidx(1):hzidx(2)) = .5;

  subplot(211)
  plot(t, real(ifft(x)*length(t)), 'r')
  ylim([-5 5])

  subplot(212)
  plot(hz, 2*abs(x(1:length(hz))), 'r-o')
  xlim([0 15])
  title("Fig 7.5| sharp edges in the frequency domain (here, the addition of an 8-9hz plateau")
                
endfunction

function [srate t n data] = fig7_6()
  srate = 1000;
  t = 0:1/srate:3;
  n = length(t);
  frex = logspace(log10(.5), log10(20), 10);
  amps = 10*rand(size(frex));

  data = zeros(size(t));

  for fi = 1:length(frex)
    data = data+amps(fi)*sin(2*pi*frex(fi)*t);
  end

  subplot(211), plot(t,data), hold on

  x = fft(data);
  hz = linspace(0, srate, n);
  

  subplot(212),
  plot(hz, 2*abs(x)/n, '-*'), hold on
  xlim([0 25])
  filterkernel = (1./(1+exp(-hz+7)));
  x = x.*filterkernel;
  title("fig 7.6 A high pass filter attenuates low-frequncy")
  
  subplot(211), plot(t, real(ifft(x)), 'r')
  legend("original", "filtered")
  
  subplot(212), plot(hz, 2*abs(x/n), 'r-o')
  legend("original", "filtered")
endfunction

function fig7_7()
  [srate t n data] = fig7_6();

  clf
  subplot(211), plot(t,data), hold on

  x = fft(data);
  hz = linspace(0,srate,n);

  subplot(212),
  plot(hz,2*abs(x)/n, "-*"),hold on
  xlim([0 25])

  filterkernel = (1-1./(1+exp(-hz+7)));
  x = x.*filterkernel;

  subplot(211),
  plot(t, real(ifft(x)), 'r')

  subplot(212),
  plot(hz, 2*abs(x/n), 'r-o')
  title("fig 7.7 A low pass filter")
endfunction

function [t filterweights, data] = fig7_8()
  [srate t n data] = fig7_6();
  srate = 1000;
  nyquist = srate/2;
  band = [4 8];
  twid = 0.2; # transition zone of 20%
  filt0 = round(3*(srate/band(1)));
  freqs = [0 (1-twid)*band(1) band(1) band(2) (1+twid)*band(2) nyquist]/nyquist;
  idealresponse = [0 0 1 1 0 0];

  filterweights = firls(filt0, freqs, idealresponse);
  filtered_data = filtfilt(filterweights,1,data);

  clf
  subplot(211),
  plot(freqs*nyquist, idealresponse),
  hold on

  filterx = fft(filterweights);
  hz = linspace(0, nyquist, floor(filt0/2)+1);

  plot(hz, abs(filterx(1:length(hz))), 'r-o')
  legend("ideal response","filter response")
  xlim([0 50])

  subplot(212)
  plot(filterweights)
  title("fig 7.8 the ideal and the real")
endfunction

function fig7_9()
  [t filterweights data] = fig7_8();
             # there are two ways to use the kernel to filter the data
             # First the filter kernel can be used in convolution
             #       Because the kernel is real valued (not complex),
             #       the hilbert transform must be applied
             #         to the result of convolution

  Lconv = length(filterweights)+length(data)-1;
  halfwavL = floor(length(filterweights)/2);
  convolution = fft(data, Lconv).*fft(filterweights,Lconv);
  
  convres = real(ifft(convolution,Lconv));
  convres = convres(halfwavL:end-halfwavL-1);

  
  clf
  filtdat = filtfilt(filterweights, 1, data);
  plot(t, data), hold on
  # plot(t, convres, 'r')
  plot(t, filtdat, 'k')
endfunction

function fig7_10()
  ## Generate a time series that increases from 5hz to 15hz over 4 seconds
  ## with randomly vayring amplitude.
  ## Convolve this time series with a 10hz complex Morlet wavelet.
  ## apply the filter-hilbert method using
  ##      a plateau shaped band-pass filter centered at 10Hz.
  ## From bouth results, extract the filterd time series and the power
  ##      and plot them simultaneously.

  srate=1000;
  t=0:1/srate:4;
  freq=[5 15];
  freqs = linspace(5,15);
  a = randn(length(freqs));
  for i = 1:length(freqs)
    signal = a(i).*sin(2*pi.*freqs(i).*t);
  end
  
  [convres pow] = convolution_resolution(10,signal,srate);
  hilsine = hilbert(signal);
  hilbert_pow = abs(hilsine).^2;

  clf
  subplot(211),
  plot(t,signal), hold on
  plot(t,hilsine,"r")
  legend("original signal","hilbert filtered")

  subplot(212),
  plot(t,pow), hold on
  plot(t,hilbert_pow,"r")
  title("frequencies calculation")
  legend("Complex Wavelt Convolution", "extracted from hilber")
  
endfunction

function [convres pow] = convolution_resolution(f,signal,srate);
                          #signalX is the spectra of signal
                          #which means the fourier transform of signal

  t=0:length(signal)-1;
  wavelet_time    = -2:1/srate:2;
  Lconv       = length(t)+length(wavelet_time)-1;
  halfwavsize = floor(length(wavelet_time)/2);
  Lconv = length(t)+length(wavelet_time)-1;
  ncyc = 6;
  wavelet_width = 2*(ncyc/(2*pi*f))^2;
                                # Complex Morlet Wavelet
  cmw = exp(1i*2*pi*f.*wavelet_time).* exp((-wavelet_time.^2)/wavelet_width);
  cmwX=fft(cmw,Lconv);
                                #normalizing the convolution
  cmwX = cmwX./max(cmwX); 

  signalX     = fft(signal,Lconv);
  convres = ifft(signalX.*cmwX);
  convres = convres(halfwavsize:end-halfwavsize-1);
  pow = 2*abs(convres);

endfunction

function fig7_10book()
                                # first define the basics
  srate = 1000;
  nyquist = srate/2;
  t=0:1/srate:4;
  peakfreq = 10;

                                # second create the chirp
  freqTS = linspace(5,15,length(t));
  centfreq = mean(freqTS);
  k = (centfreq/srate)*2*pi/centfreq;
  pow = abs(interp1(linspace(t(1),t(end),10), ...
            10*rand(1,10),t,'spline'));
  signal = pow.*sin(2.*pi.*centfreq.*t + ...
                    k*cumsum(freqTS-centfreq));

                                # third, wavelet convolution
  [convres pow] = convolution_resolution(peakfreq,signal,srate);

                                # fourth, band-pass filter
  band = [peakfreq-.5 peakfreq+.5];
  twid = 0.15;
  filtO = round(3*(srate/band(1)));
  freqs = [0 (1-twid)*band(1) band(1) band(2) (1+twid)*band(2) nyquist]/nyquist;
  ires = [0 0 1 1 0 0];
  fweights = firls(filtO, freqs, ires);
  filtdat = filtfilt(fweights,1,signal);

                                # fifth, plot results
  clf
  subplot(311), plot(t,signal), title("Fig 7.10a a chirp signal from 5 to 15hz randm amplitude")
  subplot(312), plot(t, real(convres)), hold on
  plot(t, filtdat,'r'), title("10hz filtered signals using wavelet and FIR band-pass filtering")
  ylim([-10 10])

  subplot(313)
  plot(t,pow), hold on
  plot(t,2*abs(hilbert(filtdat)), 'r')
  title("10 hz amplitude")
  legend("wavelet", "FIR finite impulse response")
                          
endfunction
