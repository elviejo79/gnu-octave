1;

function result=fig4_1()
  srate=1000;
  t=0:1/srate:10;
  n=length(t);
  ## csw = complex sine wave
  csw = exp(1i*2*pi*t);
  plot3(t,real(csw),imag(csw))
  xlabel('time');
  ylabel('real part');
  zlabel('imaginary part');
  rotate3d;
endfunction

function result=fig4_2()
  srate=1000;
  t=0:1/srate:10;
  n=length(t);

  signal = 2*sin(2*pi*3*t+pi/2);

  fouriertime = (0:n-1)/n;
  signalx = zeros(size(signal));
  for fi=1:length(signalx)
    csw = exp(-1i*2*pi*(fi-1)*fouriertime);
    signalx(fi) = sum(csw.*signal)/n;
  end
  
  ## dividing the fourier coefficients by thoe number of time points
  ## scales the coefficients to the amplitude of the original data.
  signalxf = fft(signal)/n;
  ## the highest frequency that can be measured in time series is one
  ## half of the sampling rate called nyquist frequency
  nyquistfreq=srate/2;
  ## the conversion from frequency indices to hertz can be made by
  ## linearly increasing steps from 0 to the nyquist in n/2+1 steps
  hz=linspace(0, nyquistfreq, floor(n/2)+1);

  clf
  subplot(211)
  plot(hz,2*abs(signalx(1:length(hz)))), hold on
  plot(hz,2*abs(signalxf(1:length(hz))), "ro")
  xlabel("frequencies (hz)")
  ylabel("amplitude")
  xlim([0 10])
  legend("slow fourier transfrom", "fast fourier transform")

  subplot(212)
  plot(hz,angle(signalx(1:floor(n/2)+1)))
  xlabel("frequencies (hz)")
  ylabel("phase (radians)")
endfunction
  
function result=fig4_3()
  srate=1000;
  t=0:1/srate:5;
  n=length(t);
  a=[10 2 5 8];
  f=[3 1 6 12];
  swave = zeros(size(t));
  
  for i=1:length(a)
    swave=swave+a(i)*sin(2*pi*f(i)*t);
  end

  ##fourier transform
  swavex = fft(swave)/n;
  hz = linspace(0,srate/2,floor(n/2)+1);

  clf
  
  subplot(211), plot(t,swave)
  xlabel("time (s)"),
  ylabel("amplitude")

  subplot(212)
  plot(hz, 2*abs(swavex(1:length(hz))),"ro")
  xlim([0 max(f)*1.3])
  xlabel("frequencies hz")
  ylabel("amplitude")
  
endfunction

function result=fig4_4()
  srate=1000;
  t=0:1/srate:5;
  n=length(t);
  a=[10 2 5 8];
  f=[3 1 6 12];
  swave = zeros(size(t));
  
  for i=1:length(a)
    swave=swave+a(i)*sin(2*pi*f(i)*t);
  end

  swaven = swave + randn(size(swave))*20;
  swavenx = fft(swaven)/n;

  ##fourier transform
  swavenx = fft(swaven)/n;
  hz = linspace(0,srate/2,floor(n/2)+1);

  clf
  
  subplot(211), plot(t,swaven)
  xlabel("time (s)"),
  ylabel("amplitude")

  subplot(212)
  plot(hz, 2*abs(swavenx(1:length(hz))),"ro")
  xlim([0 max(f)*1.3])
  xlabel("frequencies hz")
  ylabel("amplitude")
  
endfunction

function sinewave = sinewave(a,f,t)
  sinewave = sum(a.*sin(t*2*pi.*f),1);
endfunction

function result=fig4_5()
  srate=1000;
  t=0:1/srate:5;
  n=length(t);
  a=[10 2 5 8];
  f=[3 1 6 12];
  swave = sinewave(a',f',t);
  hz = linspace(0,srate/2,floor(n/2)+1);
  frex_idx = sort(dsearchn(hz',f'));
  swavex = fft(swave)/n;
  requested_frequencies = 2*abs(swavex(frex_idx));
  clf;
  bar(requested_frequencies)
  xlabel('freq (hz)')
  ylabel('amplitude')
  set(gca,'xtick',1:length(frex_idx),...
      'xticklabel',cellstr(num2str(round(hz(frex_idx))')))
  title("fig 5.5: bars of power")
endfunction

function result = fig4_6()
  ## consider a sine wave with constant frequency but varying apmlitude
  srate = 1000;
  t=0:1/srate:5;
  n=length(t);
  f=3;

  swave = sinewave(linspace(1,10,n),f,t);
  swavex = fft(swave)/n;
  hz = linspace(0,srate/2,floor(n/2)+1);

  clf
  title("fig 4.6: non-stationarities widen spectral represantations");
  subplot(211)
  plot(t,swave)
  xlabel('time')
  ylabel('amplitude')

  subplot(212)
  plot(hz, 2*abs(swavex(1:length(hz))))
  xlabel('frequency (hz)')
  ylabel('amplitude')
  xlim([0 10])
  
endfunction


function result=fig4_7()
  ## consider a sine wave with constant frequency but varying apmlitude
  
  ## next consider a time series with varying amplitudes nad varying frequencies.
  a = [10 2 5 8];
  f = [3 1 6 12];

  tchunks = vec2mat(t,floor(n/length(a))+1);
  
  swave=[];
  for i=1:length(a)
    swave=[swave sinewave(a(i),f(i),tchunks(i,:))];
  end
  swave=swave(1:n);

  swaveX = fft(swave)/n;
  hz = linspace(0,srate/2,floor(n/2)+1);
  
  clf
  subplot(211)
  plot(t,swave)
  xlabel('Time')
  ylabel('amplitude')
  title("4.7: Non-stationarities widen spectral represantations II")
  
  subplot(212)
  plot(hz,2*abs(swaveX(1:length(hz))),"ro-")
  xlabel('frequency (Hz)')
  xlim([0 20])
  ylabel('amplitude')
endfunction

function result=fig4_8()
  srate = 1000;
  t=0:1/srate:5;
  n=length(t);

  f = [2 10];
  ff = linspace(f(1),f(2)*mean(f)/f(2),n);
  swave = sinewave(1,ff,t);

  swaveX = fft(swave)/n;
  hz = linspace(0,srate/2,floor(n/2));
  
  clf
  subplot(211)
  title("4.8: Chirp in time and in frequency")
  plot(t,swave)
  xlabel("Time")
  ylabel("amplitude")

  subplot(212)
  plot(hz,2*abs(swaveX(1:length(hz))))
  xlim([0 30])
  xlabel("Frequency (Hz)"),
  ylabel("amplitude")
endfunction

function results = fig4_9()
  srate=100;
  t=0:1/srate:11;
  n=length(t);

  boxes = double(.02+(mod(.25+t,2))>1.5);

  boxesX = fft(boxes)/n;
  hz = linspace(0, srate/2, floor(n/2));

  clf
  subplot(211)
  plot(t,boxes)
  title("4.9: Non-sinusoidal signals have multiplied spectral forms")
  ylim([0 1.5])
  xlabel("Time")
  ylabel("amplitude")

  subplot(212)
  plot(hz,2*abs(boxesX(1:length(hz))))
  xlim([0 30])
  xlabel('Frequency (Hz)')
  ylabel('amplitude')
endfunction

function result = fig4_10()
  ## many smooth sine waves must be summed to produce a straight line
  x = (linspace(0,1,1000)>.5)+0; #+0 transforms boolean to integer

  subplot(211)
  plot(x)
  title("4.10: Edge artifacts")
  ylim([-.1 1.1])
  xlabel("Time (a.u.)")

  subplot(212)
  plot(abs(fft(x)))
  xlim([0 200])
  ylim([0 100])
  xlabel("Frequency (a.u.)")
endfunction

function result = fig4_11()
  ## show the Fourier spectra of two sinusoidal responses with edge artifacts
  srate=1000;
  t=0:1/srate:10;
  n=length(t);

  x1 = sin(2*pi*2*t+pi/2);
  x2 = sin(2*pi*2*t);

  subplot(211)
  title("4.11: Edge artifacts with cosine")
  plot(t,x1), hold on
  plot(t,x2,'r')
  xlabel('Time')
  ylabel('amplitude')

  hz = linspace(0,srate/2,floor(n/2)+1);
  x1X = fft(x1)/n;
  x2X = fft(x2)/n;                         
                         
                        
  subplot(212)
  plot(hz,2*abs(x1X(1:length(hz))),'b.-'), hold on
  plot(hz,2*abs(x2X(1:length(hz))),'r.-')
  xlabel('Frequency (hz)')
  ylabel('Amplitude')
  xlim([0 10])
  ylim([0 .001])

  result = struct("srate",srate,
                  "t",t,
                  "n",n,
                  "x1",x1,
                  "x2",x2,
                  "hz",hz,
                  "x1X",x1X,
                  "x2X",x2X);
endfunction

function result = fig4_12()
  p = fig4_11();
  ## A taper is a smooth envelope that dampens the time series
  ## near the beginning and the end, thereby
  ## effectively removig edge artifacts
  ## there are may tapper; a Hann window will be shown here
  hannwin = .5*(1-cos(2*pi*linspace(0,1,p.n)));

  clf;
  subplot(311), plot(p.t,p.x1)
  title("4.12: Sine wave(top), taper(middle), and tapered sine wave")
  subplot(312), plot(p.t,hannwin)
  subplot(313), plot(p.t, p.x1 .* hannwin)
  result = p;
  result.hannwin = hannwin;
endfunction

function result = fig4_13()
  ## applying the taper also atenuates valid signal.
  ## this potential loss of signal must be balanced with the attenuation of artifacts introduced by edges.
  edges = fig4_12();
  clf;
  
  hz = edges.hz;
  x1X = edges.x1X; # this were previously calculated in fig4_11()
  x2X = fft(edges.x1.*edges.hannwin)/edges.n;

  clf;
  plot(hz,2*abs(x1X(1:length(hz)))), hold on
  plot(hz,2*abs(x2X(1:length(hz))),'r')
  xlim([0 10])
  ylim([0 .002])
  legend("original","tapered")
  title("4.13: Tapered sine wave reducs edge artifacts")
endfunction

function result = fig4_14()
  n=50;
  x=(linspace(0,1,n)>.5)+0;
  zeropadfactor = 2;

  clf;
  subplot(211)
  plot(x)
  title("4.14: Same FFT, different N")
  ylim([-.1 1.1])
  xlabel("time (a.u.)")

  ## El patrón es:
  ## sacas la transoframa de furier a esa le llamas spectral (X)
  ## sacas los hz que va desde 0 hasta la mitad de los datos
  ## usas 2*abs( de la espectral (desde 1 hasta los hz)) para saber la frecuencia
  X1 = fft(x,n)/n;
  hz1 = linspace(0, n/2, floor(n/2)+1);

  subplot(212)
  plot(hz1, 2*abs(X1(1:length(hz1))))

  X2 = fft(x,zeropadfactor*n)/n;
  hz2 = linspace(0, n/2, floor(zeropadfactor*n/2)+1);

  hold on,
  plot(hz2, 2*abs(X2(1:length(hz2))), 'r')
  xlim([0 20])

  xlabel("freq (a.u.)")
  legend("50 point FFT", "100 point FFT")
endfunction

function result = fig4_15()
  srate = 1000;
  t=0:1/srate:1;
  f=30;
  srates = [15 20 50 200]; ## Hz

  ## "continuas" sine wave
  d = sinewave(1,f,t);

  clf;
  for i=1:4
    subplot(2,2,i)
    plot(t,d)
    hold on
    samples = round(1:srate/srates(i):length(t));
    plot(t(samples),d(samples),'r-','linewidth',2)
    legend(["sampling rate " num2str(srates(i))])
  end
  xlabel("time (s)")
  ylabel("amplitude")
  title("4.15: Subsampling produces aliasing")
endfunction

function results = trials(nTrials, noise)
endfunction

function results = fig4_16()

  % Section 4.11 Repeated instruments
  nTrials=40;
  srate=1000;
  t=0:1/srate:5;
  n=length(t);
  a=[2 3 4 2];
  f=[1 3 6 12];
  
  data = sinewave(a', f',t);
                                % create trials with noise
  dataWnoise = bsxfun(@plus, data,30*randn(nTrials,n));

  hz=linspace(0,srate/2,floor(n/2)+1);
  dataPow = zeros(nTrials, length(hz));
  hanwin = .5*(1-cos(2*pi*linspace(0,1,n)));

  for triali=1:nTrials
    temp = fft(hanwin.*dataWnoise(triali,:))/n;
    dataPow(triali,:) = 2*abs(temp(1:length(hz)));
  end

  clf
  subplot(211), plot(t, dataWnoise), hold on
  plot(t, mean(dataWnoise), 'k')
  subplot(212), plot(hz, dataPow), hold on
  plot(hz, mean(dataPow), 'k', 'linewidth',5)
  xlim([0 20]);
  title("fig 4.16 Averaging improves results");

  results=struct('dataPow',dataPow,'hz',hz);
endfunction 

function result = fig4_17()
  ## Sect 4.12 | Signal to Noise ratio
  p = fig4_16(); %previous excercise
  snr = mean(detrend(p.dataPow')')./...
        std(detrend(p.dataPow')');

  clf
  plot(p.hz,snr)
  xlim([0 20])
  title("Fig 4.17 |Empirical SNR (mean divided by standard deviation)")
endfunction

## I can' make the figu 5_18 because octave doesn't
## hed dpss function.
## maybe I'll try to implement it some day.

function result = fig4_19()
  ## Section 4.14 The Inverse Fourier transform
  xTime = randn(20,1);
  xFreq = fft(xTime)/length(xTime);
  t = (0:length(xTime)-1)' / length(xTime);

  recon_data = zeros(size(xTime));
  for fi=1:length(xTime)
    sin_wave = xFreq(fi)*exp(1i*2*pi*(fi-1).*t);
    recon_data = recon_data + sin_wave;
  end
  
  clf
  plot(xTime, '.-'), hold on
  plot(real(recon_data),'ro')
  legend("Original data", "Inverse Fourier Transform")
endfunction

function results=fig4_20()
  ## basics
  srate=1000;
  t = 0:1/srate:6;
  N = length(t);
  f = [1 5]; % in Hz

  ## In this excercise: a) generate a linear chirp and b) compute its inverse 
  ## Generate a linear chirp
  ff = linspace(f(1), f(2)*mean(f)/f(2),N);
  data = sin(2*pi*t.*ff);
  
  ## secord, compute Fourier Transform
  dataX = fft(data);

  ## third, shuffle phases(here, shifted by 10)
  phases = angle([dataX(10:end) dataX(1:9)]);
  shuffdata = abs(dataX).*exp(1i*phases);

  ## fourth, reconstruct the signal
  newdata = ifft(shuffdata);

  ## fifth, plot the results
  subplot(211)
  plot(t,data), hold on
  plot(t,real(newdata),'r')
  xlabel('time (s)')
  ylabel('amplitude')
  legend("original", "phase scrambled")
  title("4.20| Same power, shuffeld phases")
  
  subplot(212)
  hz = linspace(0, srate/2, floor(N/2)+1);
  plot(hz, 2*abs(dataX(1:length(hz))/N), 'o'), hold on
  plot(hz, 2*abs(shuffdata(1:length(hz))/N), 'r')
  xlim([0 20])
  xlabel('frequency (Hz)')
  ylabel('amplitude')
  legend("original", "phase scrambled")
endfunction

function result=trials(nTrials, noiseLevel)
  srate=1000;
  t=0:1/srate:5;
  n=length(t);
  a=[2 3 4 2];
  f=[1 3 6 12];
  data=sinewave(a', f',t);

  ## create trials with noise
  dataWnoise = bsxfun(@plus, data, noiseLevel*randn(nTrials,n));
  hz = linspace(0, srate/2, floor(n/2)+1);
  dataPow = zeros(nTrials, length(hz));

  hanwin = .5*(1-cos(2*pi*linspace(0,1,n)));

  for triali = 1:nTrials
    temp = fft(hanwin.*dataWnoise(triali,:))/n;
    dataPow(triali,:) = 2*abs(temp(1:length(hz)));
  end

  result = struct("t",t,"dataWnoise",dataWnoise,"hz",hz,"dataPow",dataPow);
endfunction

function result=fig4_21()
  clf;
  t1 = trials(10,10);
  subplot(221), plot(t1.t, t1.dataWnoise), hold on, plot(t1.t, mean(t1.dataWnoise),'k'), title("10 trials, 10x noise")
  subplot(223), plot(t1.hz, t1.dataPow), hold on, plot(t1.hz, mean(t1.dataPow),'k'), ylim([0 2]), xlim([0 20])
  
  t2 = trials(100,10);
  subplot(222), plot(t2.t, t2.dataWnoise), hold on, plot(t2.t, mean(t2.dataWnoise),'k'), title("1000 trials, 10x noise")
  subplot(224), plot(t2.hz, t2.dataPow), hold on, plot(t1.hz, mean(t2.dataPow),'k'), ylim([0 2]), xlim([0 20])
endfunction

function result=fig4_22()
  srate=1000;
  t=0:1/srate:2;
  N=length(t);
  
  swave=sinewave(1,10,t);
  swaveX = fft(swave,N)/N;
  hz = linspace(0,srate/2,floor(N/2)+1);
  power = 2*abs(swaveX(1:length(hz)));

  clf
  subplot(221)
  bar(hz, power), title("Fig 4.22"), ylim([0 1]), xlim([5 15])

  p = padder(2,N,srate,swave);
  subplot(222), bar(p.Hz, p.Power), title(p.title), ylim([0 1]), xlim([5 15])

  p = padder(3,N,srate,swave);
  subplot(223), bar(p.Hz, p.Power), title(p.title), ylim([0 1]), xlim([5 15])

  p = padder(1.5,N,srate,swave);
  subplot(224), bar(p.Hz, p.Power), title(p.title), ylim([0 1]), xlim([5 15])
  
endfunction

function results=padder(padding,N,srate,swave)
  pSwaveX=fft(swave,padding*N)/(padding*N);
  pHz = linspace(0,srate/2,floor(padding*N/2));
  pPower = 2*abs(pSwaveX(1:length(pHz)));
  results = struct("Hz",pHz,
                   "Power",pPower,
                   "title",num2str(padding*N));
endfunction
