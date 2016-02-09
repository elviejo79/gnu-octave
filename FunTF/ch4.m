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
