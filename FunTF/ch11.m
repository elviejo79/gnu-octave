function ch11()
  clf;
  fig11_4();
  clear -all;
endfunction


function fig11_1()
  srate=1000;
  t=0:1/srate:9;
  n=length(t);

                                # create signals
  f=[10 14 8];
  k1 = (f(1)/srate)*2*pi/f(1);
  sigA = sin(2*pi.*f(1).*t + k1*cumsum(5*randn(1,n)))+randn(size(t));
  sigB = sigA + sin(2*pi.*f(2).*t + k1*cumsum(5*randn(1,n)));
  sigA = sigA + sin(2*pi.*f(3).*t + k1*cumsum(5*randn(1,n)));

                                # show power of each channel
  hz = linspace(0,srate/2,floor(n/2)+1);
  sigAx = fft(sigA)/n;
  sigBx = fft(sigB)/n;

  clf
  subplot(211)
  plot(hz, 2*abs(sigAx(1:length(hz)))), hold on
  plot(hz, 2*abs(sigBx(1:length(hz))),'r'),
  xlim([5 20]),
  xlabel("Frquency (hz)"),
  ylabel("Amplitude");
  legend("signal A", "signal ");

                                # spectral coherence
  specX = abs(sigAx.*conj(sigBx)).^2;
  spectcoher = specX./(sigAx.*sigBx);
  
  subplot(212)
  plot(hz, abs(spectcoher(1:length(hz)))),
  xlim([5 20]),
  xlabel("Frquency (hz)"),
  ylabel("coherence"),
  title("fig 11.1| si recarga? spectral power and spectra coherence");
      
endfunction

function [covar d] = fig11_2()
  ## generate a covarying time series
  ## compute its covariance matrix and then
  ## display the covariance matrix as an image

                                # covariance of data
  v = rand(10);
  c = chol(v*v');
  n = 10000;
  d = randn(n, size(v,1))*c;
                                # subtract mean and compute covariance
  d = bsxfun(@minus,d,mean(d,1));
  covar = (d'*d)./(n-1);
  imagesc(covar),
  xlabel("channel number")
  ylabel("channel number")
  title("Fig 11.2| Covariance matrix")

  
endfunction

function [covar d pc] =fig11_3()
  [covar d] = fig11_2();
  [pc, ev] = eig(covar);

                                #resort components
  pc = pc(:,end:-1:1);
                              # exctract eigen values and convert to %
  ev = diag(ev);
  ev = 100*ev(end:-1:1)./sum(ev);
  clf;
  plot(ev, '-o'),
  xlabel("Component number")
  ylabel("% variance ")
endfunction

function fig11_4()
  [covar d pc] = fig11_3();

  clf;
  subplot(211), plot(d),
  title("all the channels data"),
  set(gca,'xlim',[0 100]);
  
  subplot(212)
  plot(pc(:,1)'*d'), hold on
  plot(pc(:,2)'*d', 'r')
  set(gca,'xlim',[0 100])
  legend("Princial Component 1", "Princial Component 2")
  title("fig 11.4| note that the first PComponents captures much of the variance");
endfunction
