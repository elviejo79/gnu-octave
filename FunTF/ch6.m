1;

function fig6_2();
  t=-1:.01:1; f=10;
  sinewave = cos(2*pi*f*t);
  w = 2*(5/(2*pi*f))^2;
  gaussian = exp((-t.^2)/w);
  mwavelet = sinewave.*gaussian;
  plot(t,mwavelet);
  title("6.2 Morlet wave");
endfunction

function fig6_3()
  ## A complex morlet wavelet si necessary to calculate powerd and phase
  t=-1:.01:1; f=10;
  sinewave = cos(2*pi*f*t);
  w = 2*(5/(2*pi*f))^2;
  csw = exp(1i*2*pi*f*t);
  gaussian = exp((-t.^2)/w);
  
  mwavelet = csw.*gaussian;
  plot3(t,real(mwavelet),imag(mwavelet),'-o');
  xlabel('Time (s)'), ylabel('real part'), zlabel('imaginary part')
  rotate3d;
  title("fig 6.3 Complex Morlet wavelet in 3D")
endfunction

function [t chirpTS wtime halfwavL cmw Lconv]= fig6_5()
  srate=1000; f=[2 8];
  t=0:1/srate:6;
  n=length(t);
  chirpTS = sin(2*pi.*linspace(f(1),f(2)*mean(f)/f(2),n).*t);;

  ## create complex Morlet wavelet
  wtime = -1:1/srate:2; %wavelet time
  w = 2*(4/(2*pi*5))^2;
  ## sinusoidal .* gaussian
  cmw = exp(1i*2*pi*5.*wtime).*exp((-wtime.^2));

  ## half of the length of the wavelet
  halfwavL = floor(length(wtime)/2);

  ## zero-pad chirp
  chirpTSpad = [zeros(1,halfwavL) chirpTS zeros(1,halfwavL)];

  ## run convolution
  convres = zeros(size(chirpTSpad));
                                # shouldn'this be length(chirpTSpad) ?
                                # no, because chirpTSpad goes from -1 to 2 seconds
  for i=halfwavL+1:length(chirpTS)+halfwavL-1
    ## each time shifthing the kernel over by one unit of time
    kernel = chirpTSpad(i-halfwavL:i+halfwavL);
    ## at each time point, compute dot product
    convres(i)=sum(kernel.*cmw);
  end

  ## cut off edges
  convres = convres(halfwavL:end-halfwavL-1);

  clf
  subplot(311), plot(t,chirpTS);
  title("Fig 6.5: Result of convolution between a 5hz wavelet and a 2-8 chirp")
  subplot(312), plot(t,abs(convres));

  Lconv = length(t)+length(wtime)-1;
  ## convolution can be performed by
  ## computing the fourier transforms and the kernel and wavelet
  kernelX = fft(chirpTS,Lconv);
  waveletX = fft(cmw,Lconv);
  ## multplying their frequency scpectra
  ## computing the inverse fourier transform
  convres2 = ifft(kernelX.*waveletX);
  convres2 = convres2(halfwavL:end-halfwavL-1);
  
  subplot(313), plot(t,abs(convres2));
  title("Fig 6.5 Convolution with fourier transform")
endfunction

function fig6_6()
  [t chirpTS wtime halfwavL cmw] = fig6_5();
  Lconv = length(t)+length(wtime)-1;
  ## Frequency-domain amplitude scalling is simple and effective.
  ## it involves scaling the frequency representation of the wavelet
  ## to a maximum of 1.0 before multiplying he spectra of the wavelet
  ## and time series together.
  cmwX = fft(cmw,Lconv);
  cmwX = cmwX./max(cmwX);

  convres4 = ifft(fft(chirpTS,Lconv).*cmwX);
  convres4 = convres4(halfwavL:end-halfwavL-1);

  clf
  plot(t, chirpTS), hold on
  plot(t, 2*abs(convres4),'r');
  title("Fig 6.6|Convolution result with appropriate amplitude scalling")
endfunction

function [t chirpTS wtime halfwavL Lconv frex tf chirpTSX] = fig6_7()
  [t chirpTS wtime halfwavL cmw Lconv]= fig6_5();
  nfrex = 30;
  frex = logspace(log10(2), log10(30), nfrex);
  tf = zeros(nfrex, length(chirpTS));
  chirpTSX = fft(chirpTS, Lconv);

  for fi=1:nfrex
    w = 2*(5/(2*pi*frex(fi)))^2;
    convolution = exp(1i*2*pi*frex(fi).*wtime).*exp((-wtime.^2)/w);
    cmwX = fft(convolution, Lconv);
    cmwX = cmwX./max(cmwX); # scale the convolution spectra (frequencies) to 1
    convres = ifft(chirpTSX.*cmwX);
    tf(fi,:) = 2*abs(chop_edges(convres,halfwavL));
  end

  clf
  subplot(211)
  plot(t,chirpTS)
  
  subplot(212)
  contourf(t,frex,tf,40,'linecolor','none')
  xlabel("time (s)")

  ylabel("Frequency (Hz)")
  title("Fig 6.7| Result of convolution between a family of complex morlet waves and a 2-8chirp")
endfunction

function data = chop_edges(convres,halfwavL)
  data = convres(halfwavL:end-halfwavL-1);
endfunction

function fig6_8()
  [t frex tf] = fig6_7();
                                # define time points to plot
  size(t')
  size((1:.2:4.5)')
  t2plot = dsearchn(t', (1:.2:4.5)');
                                #and the frequency to plot
  f2plot = dsearchn(frex',6);
  plot(t,tf(f2plot,:)), hold on
  plot(t(t2plot), tf(f2plot, t2plot), 'ro-');
  legend("original", "down-sampled");
  title("Fig 6.8| temporally downsampling results is often a good idea")
end

function fig6_9()
  ## Section 6.7 When the Number of CYCles (ncyc) is largerd the gausisian is wider..
  [t chirpTS wtime halfwavL Lconv frex tf chirpTSX] = fig6_7();
  freq2use = 5;

  ## Width the la complex Morlet Wavelet
  ncyc = 7;
  w = 2*(ncyc/(2*pi*freq2use))^2;
                                # Complex Morlet Wavelet
  cmw7 = exp(1i*2*pi*freq2use.*wtime).* exp((-wtime.^2)/w);
  convres = ifft(chirpTSX.*fft(cmw7,Lconv));
  pow7 = abs(chop_edges(convres,halfwavL));

  ncyc=3;
  w = 2*(ncyc/(2*pi*freq2use))^2;
                                # Complex Morlet Wavelet 3 num of cycles
  cmw3 = exp(1i*2*pi*freq2use.*wtime).*exp((-wtime.^2)/w);
  convres = ifft(chirpTSX.*fft(cmw3,Lconv));
  pow3 = abs(chop_edges(convres,halfwavL));

  clf
  subplot(221)
  plot(wtime,real(cmw7)), hold on
  plot(wtime,real(cmw3), 'r')
  title("two morlet wavelets with Gaussian widths 3 and 7")
  legend("7 gaussian width", "3 gaussian width")
  xlabel("wtime (s)")

                                #duplicating sampling rate here.
                                #becaus I don't want to change the calls
  srate=1000;
  subplot(222)
  hz = linspace(0, srate/2, floor(length(wtime)/+1));
  x7 = 2*abs(fft(cmw7)); # what is going on here?
  x3 = 2*abs(fft(cmw3));
  plot(hz,x7(1:length(hz))), hold on
  plot(hz,x3(1:length(hz)), 'r')
  xlim([0 20])
  title("power spectra")
  legend("7 gaussian width", "3 gaussian width")
  xlabel("Frequency in (hz)")

  subplot(212)
  plot(t,pow7), hold on
  plot(t,pow3,'r')
  title("Fig 6.9 Convolution with the chirp")
  legend("7 gaussian width", "3 gaussian width")
  xlabel("time (s) of the chirp" )
  
                                 
endfunction

function fig6_10()
  [t chirpTS wtime halfwavL Lconv frex tf chirpTSX] = fig6_7();

  ## allo Number of CYCles to vary over a range as a function of crequency.
  nfrex=9;
  frex = logspace(log10(2), log10(20), nfrex);
  ncyc = logspace(log10(3), log10(12), nfrex);

  cmw_fam= zeros(nfrex, length(wtime));
  for fi=1:nfrex
    w = 2*(ncyc(fi) / (2*pi*frex(fi)))^2;
    cmw_fam(fi,:) = exp(1i*2*pi*frex(fi).*wtime).*exp((-wtime.^2)/w);
  end

  for i=1:9
    subplot(3,3,i)
    plot(wtime,(real(cmw_fam(i,:))))
    xlim([-1 1])
    title(["freq ", num2str(frex(i)), " Hz", "Num of Cycles", num2str(ncyc(i))])
  end
endfunction

function fig6_11()
  srate=1000;
  t=0:1/srate:6;

                                # create signal
  f = [6 14 25 40 70];
                                # note the widely varying amplitudes
  a=[.001234 1.234 1234 123400 12340000];
                                # relevant amplitue modulations
  m=[-.1 .1 -.2 .2 -.3];
  signal = zeros(size(t));

  for i=1:length(f)
                                #compute 'base' signal
    signal = signal+a(i).*sin(2*pi*f(i).*t);
                                #compute time-limited modulation
    extrasignal= m(i)*a(i)*sin(2*pi*f(i).*t) .* exp(-(t-2).^2);
    signal = signal + extrasignal;
  end

  ## tf=zeros(length(f),length(t)); 
  for i=1:length(f)
    [convres pow] = convolution_resolution(f(i),signal,srate);
    tf(i,:)=pow;
  end

  ## plot the non-dB-normalized result
  subplot(121)
  plot(t, tf)
  xlabel("Time (s)"), ylabel("Amplitude (raw")
  title("Non-normalized amplitude")

  ## compute db relativo to baseline
  bidx=dsearchn(t',[2 4]');
  baseMean=mean(tf(:,bidx(1):bidx(2)),2);
  db=10*log10(bsxfun(@rdivide, tf, baseMean));
  
  ## convert to dB and plot again
  subplot(122)
  plot(t,db);
  ylim([-2 2])
  xlabel('Time (s)'), ylabel('Amplitude (dB)')
  title('dB-normalized amplitude')

endfunction

function [convres pow] = convolution_resolution(f,signal,srate)
                          #signalX is the spectra of signal
                          #which means the fourier transform of signal

  t=0:length(signal)-1;
  wavelet_time    = -2:1/srate:2;
  Lconv       = length(t)+length(wavelet_time)-1;
  halfwavsize = floor(length(wavelet_time)/2);
  Lconv = length(t)+length(wavelet_time)-1;
  ncyc = 10;
  wavelet_width = 2*(ncyc/(2*pi*f))^2;
                                # Complex Morlet Wavelet
  cmw = exp(1i*2*pi*f.*wavelet_time).* exp((-wavelet_time.^2)/wavelet_width);
  cmwX=fft(cmw,Lconv);
                                #normalizing the convolution
  cmwX = cmwX./max(cmwX); 

  signalX     = fft(signal,Lconv);
  convres = ifft(signalX.*cmwX);
  pow = 2*abs(chop_edges(convres,halfwavsize));

endfunction

function fig6_12()
  ##First generate 20 repeated ea
  wanted_signals=20;
  signal_gen = @(f,t,phase) sin(2*pi*f.*t+phase');
  phases=2*pi*rand(wanted_signals);
  srate=1000;
  t=0:1/srate:6;

  ## First, generate 20 repeated measurements ("trials" of atime series
  for i=1:length(phases)
    signals(i,:)=signal_gen(1,t,phases(i));
  end

  clf
  subplot(221)
  plot(t,signals);
  ylim([-1 1])
  title("a) all trials")

  ## Second Perform the convolution on each trial separately
  for i=1:length(phases)
    [convres pow] = convolution_resolution(1,signals(i,:),srate);
    pows(i,:) = pow;
  end


  subplot(222)
  plot(t,mean(pows))
  title("trial average response")
  ylim([-1 1])

  subplot(223)
  plot(t,pows)
  title("power all trials")
  ylim([-1 1])

  
  ##Third, compute another trial-averaged time series
  ## first average the time series over the trials in the time domain
  ## and then perform convolution
  subplot(224)
  [convres pow_of_mean] = convolution_resolution(1,mean(signals),srate);
  plot(t,pow_of_mean)
  ylim([-1 1])
  title("power from the trial average")

  

endfunction
