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

function [t frex tf] = fig6_7()
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
  title("Fig 6.8| temporally downsampling results is ofteng a good idea")
end
