1;

function [srate, t, n, signal ] = chirp()
  srate=1000; t=0:1/srate:5; n=length(t);
  f = [30 3 6 12];

  tchunks = round(linspace(1,n,length(f)+1));
  signal = 0;
  for i=1:length(f)
    signal = cat(2,signal,...
               sin(2*pi*f(i)*t(tchunks(i):tchunks(i+1)-1)));
  end
endfunction

function [hz, t, centimes, tf, fftWidth, signal] = fig5_1()
  [srate, t, n, signal ] = chirp();
  clf
  subplot(211), plot(t,signal),
  title("Fig 5.1 A time varying signal and its time-frequency representation")

  ## We will compute the FFT in successive windows of specific width

  fftWidth_ms = 1000;
  fftWidth = round(fftWidth_ms/(1000/srate)/2);

  Ntimesteps = 10; %number of time widths
  
  ## Center Time Points (the center of the window)
  centimes = round(linspace(fftWidth+1, ...
                            n-fftWidth, Ntimesteps));

  # Note that the data in each time window should be tapered before the Fourier Transform
  hz = linspace(0, srate/2, fftWidth-1);
  tf = zeros(length(hz), length(centimes));
  hwin = .5*(1-cos(2*pi*(1:fftWidth*2) / (fftWidth*2-1))); #Hann taper

  for ti = 1:length(centimes)
    x = fft(hwin.*signal(centimes(ti)-fftWidth: centimes(ti)+fftWidth-1)) / ...
        fftWidth*2;
    tf(:,ti) = 2*abs(x(1:length(hz)));
  end

  subplot(212)
  contourf(t(centimes), hz, tf,1);
  ylim([0 50])
  set(gca,'clim',[0 1])
  xlabel('Time (s)')
  ylabel('Freq (Hz)')
  title("Time-frequency plot also called time-frequency power(/amplitude) plot")
endfunction

function fig5_2()
  [hz, t, centimes, tf] = fig5_1();
  
  ## por que se hace con la transpuesta?
  freq2plot = dsearchn(hz', 6);
  plot(t(centimes), tf(freq2plot,:), '-o')
  title("5.2: A time slice of frequency 6")
endfunction

function list5_2()
  ## keep only 30 evenly spaced frequencies
  [hz, t, centimes, tf, fftWidth, signal, centimes] = fig5_1();

  freqs2keep = linpspace(0,hz(end),30);
  freqsidx = dsearchn(hz', freqs2keep');
  hanwing = .5*(1-cos(2*pi*(1:fftWidth*2)/...
                      (fftWidth*2-1)));
  
  tf = zeros(length(freqs2keep),length(centimes));

  for ti=1:length(centimes)
    temp = signal(centimes(ti)-fftWidth : centimes(ti)+fftWidth-1);
    x = fft(hanwin.*temp)/...
        fftWidth*2;
    tf(:,ti) = 2*abs(x(freqsidx));
  end
  
  subplot(212)
  contourf(t(centimes), hz, tf,1);
  ylim([0 50])
  set(gca,'clim',[0 1])
  xlabel('Time (s)')
  ylabel('Freq (Hz)')
  title("Time-frequency plot also called time-frequency power(/amplitude) plot")
  
endfunction

function lists5_2()
  ## This code illustrates the procedure of generatig
  ## an arbitrary number of time-frequency window sizes.
  
  [hz, t, centimes, tf, fftWidth, signal, centimes] = fig5_1();
  
  NtimeWidths = 5;
  fftWidth_ms = linspace(1300,500, NtimeWidths);
  fftWidth = round(fftWidth_ms./(1000/srate)/2);
  
  Ntimesteps = 10;
  centimes = round(linspace(max(fftWidth)+1, ...
                            length(t)-max(fftWidth), ...
                            Ntimesteps));
  
  ## the time windows become smaller with increasing frequency.
  ## this increasing temporal precision at the expense of frequency resolution
  f2keep = linspace(1,50,40);
  freqsPerBin = ceil((f2keep./max(f2keep))*NtimeWidths);

  tf = zeros(length(f2keep), length(centimes));
  for ti=1:length(centimes)
    for fi=1:NtimeWidths
      ## find appropriate frequencies in this bin
      hz=linspace(0,srate/2,fftWidth(fi)-1);
      freqsidx = dsearchn(hz', f2keep(freqsPerBin==fi)');

      ## compute windowed FFT
      hanwin = .5*(1-cos(2*pi*(1:fftWidth(fi)*2/fftWidth(fi)*2-1)));
      temp = signal(centimes(ti)-fftWidth(fi):centimes(ti)+fftWidth(fi)-1);
      x = fft(hanwin.*temp)/fftWidth(fi)*2;

      ## put signal in TF matrix
      tf(freqsPerBin==fi,ti) = 2*abs(x(freqsidx));
    end
  end
endfunction

function list5_3()
  ## for signal 'signal' of length 'n'
  for ti=1:n
    ## determine points to use
    tmax = min([ti-1, n-ti, round(n/2)-1]);
    pnts = -tmax:tmax;
    idx = rem(n+pnts,n)+1;

    ## multiply forward and backward points
    wig(idx, ti) = signal(ti+pnts .* signal(ti-pnts));
  end

  ## take Fourier transform
  wig =2*abs(fft(wig)/size(wig,1));
endfunction

function ex1()
  ## Generate a 5 second chirp that spans from 1Hz to 40hz with 100 Hz sampling rate
  srate=100;
  t=0:1/srate:5;
  n = length(t);
  f=[1 40];
  ff = linspace(f(1),f(2)*mean(f)/f(2),n);
  signal = sin(2*pi.*ff.*t);


  FFTwidths = [50 100 200 500];
  for i=1:length(FFTwidths)
    ##perform a short time fourier analysis
    [tf,hz,centimes]=short_time_ft(signal, n, srate, FFTwidths(i));
    subplot(2,2,i), contourf(t(centimes),hz,tf,1)
  end

endfunction

function [tf,hz,centimes] = short_time_ft(signal, n, srate, fftWidth_ms)
  ##perform a short time fourier analysis
  fftWidth=round(fftWidth_ms/(1000/srate)/2);
  ##number of time witdths
  Ntimesteps=20;
  centimes = round(linspace(fftWidth+1, n-fftWidth, Ntimesteps));
  hz = linspace(0, srate/2, fftWidth-1);
  tf =zeros(length(hz), length(centimes));
  hannwin = .5*(1-cos(2*pi*(1:fftWidth*2)/(fftWidth*2-1)));

  for ti=1:length(centimes)
    window = signal(centimes(ti)-fftWidth:centimes(ti)+fftWidth-1);
    x = fft(hannwin.*window)/fftWidth*2; # maybe fftWidth*2 is the same as length(window)
    tf(:,ti) = 2*abs(x(1:length(hz)));
  end
end

function ex2()
  [srate, t, n, signal ] = chirp();
  clf
  subplot(211), plot(t,signal),
  title("Fig 5.1 A time varying signal and its time-frequency representation")

  numWindows = 10; # How many windows do we want?
  winWidth_ms = 1000; # in mili seconds
  winWidth = winWidth_ms/1000*srate;
  winWidth = 2*ceil(winWidth/2); # trick to round to nearest pair
  ##han tapper
  

  ##Lets assume that each window is cut at the center
  idxCenters = ceil(linspace(1+ winWidth/2,n- winWidth/2,numWindows));
  idxLefts = idxCenters.- winWidth/2;
  idxRights = idxCenters.+ winWidth/2;
  
  hz = linspace(0, srate/2, winWidth-1);
  tf = zeros(length(hz), length(idxCenters));
  hwin = .5*(1-cos(2*pi*(0:winWidth) / (winWidth-1))); #Hann taper
  for ti =1:length(idxCenters)
    x = fft(hwin.*signal(idxLefts(ti):idxRights(ti))) / winWidth*2;
    tf(:,ti) = 2*abs(x(1:length(hz)));
  end
  
  subplot(212)
  contourf(t(idxCenters), hz, tf,1);
  ylim([0 20])
  xlim([0 5])
  set(gca,'clim',[0 1])
  xlabel('Time (s)')
  ylabel('Freq (Hz)')
  title("Time-frequency plot also called time-frequency power(/amplitude) plot")

 ##  ##############################
 ##  ##based on this document we get each window
 ##  ##http://www.ee.columbia.edu/~marios/matlab/mtt.pdf
 ##  ##disp(["[" sprintf("signal(%d:%d); " , [idxLefts ; idxRights]) "]"])
 ##  windows = eval(["[" sprintf("signal(%d:%d); " , [idxLefts ; idxRights]) "]"]);
 ##  size(windows)
 ##  signalX = fft(hwin.*windows);
 ##  pows = 2*abs(signalX(1:length(hz)));

 ##  clf
 ##  plot(hz,pows);
 ## ##contourf(t(idxCenters), hz, pows,1);
 ##  ##ylim([0 50])
 ##  set(gca,'clim',[0 1])
 ##  xlabel('Time (s)')
 ##  ylabel('Freq (Hz)')
 ##  title("Time-frequency plot also called time-frequency power(/amplitude) plot")

endfunction
