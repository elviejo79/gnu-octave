pkg load signal;
1;


function y =chirp(freqTS, srate, t )
  ## srate=1000;
  ## t=0:1/srate:5;
  
  centfreq = mean(freqTS);
  k = (centfreq/srate)*2*pi/centfreq;
  y = sin(2*pi.*centfreq.* + k*cumsum(freqTS-centfreq));

endfunction

function fig8_1()
  srate=1000;
  t=0:1/srate:5;
  linearly_increasing = linspace(1,20,length(t));
  smoothed_random_sequence = abs(interp1(linspace(t(1),t(end),10),10*rand(1,10),t,'spline'));
  triangle = abs(mod(t,2)-1)*10;
  
  y = chirp(linearly_increasing,srate,t);
  clf
  subplot(321), plot(t,linearly_increasing)
  subplot(322), plot(t,y)
  subplot(323), plot(t,smoothed_random_sequence)
  subplot(324), plot(t,chirp(smoothed_random_sequence,srate,t))
  subplot(325), plot(t,triangle),xlabel('Time (s)'), ylabel('Frequency (Hz)')
  subplot(326), plot(t,chirp(triangle,srate,t)),xlabel('Time (s)'), ylabel('Amplitude')
                            
endfunction

function [srate t n] = basics()
  srate=1000;
  t=0:1/srate:5;
  n=length(t);
endfunction

function [srate t signal instFreq1]=fig8_2()
  [srate t n] = basics();
  f = [4 6];
  ff = linspace(f(1),f(2)*mean(f)/f(2),n);
  signal = sin(2*pi.*ff.*t);
  phases = angle(hilbert(signal));
  angVeloc = diff(unwrap(phases));
  instFreq1 = srate*angVeloc/(2*pi);
  instFreq1(end+1) = instFreq1(end);
  
  subplot(211), plot(t,signal)
  subplot(212), plot(t,instFreq1)
endfunction

function [srate t data instFreq1] = fig8_3()
  [srate t signal instFreq1]=fig8_2();
  a=[10 2 5 8]/10;
  f=[3 1 6 12];
  background = zeros(size(t));
  for i=1:length(a)
    background = background + ...
                 a(i)*sin(2*pi*f(i)*t);
  end

  data = signal + background;
  instFreq2 = srate*diff(unwrap(angle(hilbert(data))))/(2*pi);

  subplot(211), plot(t,data)
  subplot(212), plot(t(1:end-1),instFreq2)
  
endfunction

function fig8_4()
  [srate t data instFreq1] = fig8_3();
  clf
  plot(t, angle(hilbert(data)));
  title("Fig 8.4 non-monotonic phases produce uninterpretable frequencies")
endfunction

function [convres pow] = convolution_resolution(f,signal,srate, ncyc);
                          #signalX is the spectra of signal
                          #which means the fourier transform of signal

  t=0:length(signal)-1;
  wavelet_time    = -2:1/srate:2;
  Lconv       = length(t)+length(wavelet_time)-1;
  halfwavsize = floor(length(wavelet_time)/2);
  Lconv = length(t)+length(wavelet_time)-1;
  
  if isempty(ncyc)
    ncyc = 6; ## six will be the default value if Number of Cycles is not provided
  end
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


function fig8_5()
  [srate t data instFreq1] = fig8_3();
  wtime = -2:1/srate:2;

  [convres pow]= convolution_resolution(5,data, srate, 5);
  instFreq3 = diff(unwrap(angle(convres)));
  instFreq3(end+1) = instFreq3(end);

  [convres pow]= convolution_resolution(5,data, srate, 10);
  instFreq3b = diff(unwrap(angle(convres)));
  instFreq3b(end+1) = instFreq3b(end);

  clf
  plot(t, instFreq1), hold on,
  plot(t, srate*instFreq3/(2*pi), 'r')
  plot(t, srate*instFreq3b/(2*pi), 'k')
  title("Fig 8.5 band-pass filtering improves estimation of instantaneous frequencies")
  legend("pure", "contaminated and filtered with 5 ncyc", "contaminated / filtered 10 " )
endfunction

function [srate t signal instFreq1 data phases] = fig8_6()
  [srate t signal instFreq1]=fig8_2();
  data = signal + randn(size(signal));
  phases = angle(hilbert(data));
  angVeloc = diff(unwrap(phases));
  instFreq = srate*angVeloc/(2*pi);
  instFreq(end+1) = instFreq(end);

  subplot(311), plot(t,data), ylabel("Amplitude")
  subplot(312), plot(t, phases), ylabel("Phase")
  subplot(313), plot(t, instFreq), ylabel("Freq. Hz"), xlabel("Time (s)")
  title("Fig 8.6 Noise can cause uninterpretable instantaneous frequencies")
  
endfunction

function fig8_7()
  [srate t signal instFreq1 data phases] = fig8_6();

  [convres pow]= convolution_resolution(5,data, srate, 10);
  instFreq4 = diff(unwrap(angle(convres)));
  instFreq4(end+1) = instFreq4(end);
  freq4 = srate*instFreq4/(2*pi);
  
  clf
  plot(t, instFreq1), hold on
  plot(t, freq4, 'r')
  legend("no noise", "with noise")
  title("Fig 8.7| a  filter effectively reduces the noise that makes instantaneos frequency unintrpretable")
endfunction

function fig8_8()
  [srate t n] = basics();
  signal = zeros(size(t));
  for i=1:4
    f = rand(1,2)*10 + i*10;
    ff = linspace(f(1), f(2)*mean(f)/f(2),n);
    signal = signal + sin(2*pi.*ff.*t);
  end

  hz=linspace(0,srate/2,floor(n/2)+1);
  x = fft(signal)/n;

  clf
  subplot(311), plot(t,signal)
  subplot(312), plot(hz, 2*abs(x(1:length(hz)))), legend("original")
  xlim([0 60])
  title("Fig 8| A time series and its power spectrum, awaiting EMD")

  imfs = emdx(signal, 4);
  f = fft(imfs, [],2)/n;
  
  subplot(313)
  plot(hz, abs(f(:, 1:length(hz))).^2)
  xlim([0 60])
  legend("EMD")
  title("Fig 8.9 The first four intrinsic modes of the EMD recoer the original frequency, though not perfectly")
endfunction

