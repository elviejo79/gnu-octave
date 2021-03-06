1;

function [t signal data] = fig10_1()
           # power line noise is well defined the amplitude maybe high
           # but its frequency is between 50 or 60hz
  srate = 1000;
  t=0:1/srate:3;
  n = length(t);

  signal = .2*sin(2*pi*.8*t) + sin(2*pi*6*t);
  noise = 100*sin(2*pi*50*t);
  data = signal + noise;

  dataX = fft(data);
  hz = linspace(0,srate,length(t));
                   # create low-pass filter and apply to data spectrum
  filterkernel = (1-1./(1+exp(-hz+40)));

  dataX = dataX.*filterkernel;
                                #2 b/c one sided filter
  data2 = real(2*ifft(dataX));

  clf
  subplot(211), plot(t,data), legend("signal+noise")
  subplot(212), plot(t,data2), hold on, plot(t,signal,'r'), legend("ifft(fft(signal).*fft(filterKernel)", "original")
endfunction

function fig10_2()
  [t signal data] = fig10_1();
  
  noise = randn(size(signal));
  data = signal + noise;
  d=9;
  dataMean = zeros(size(data));
  for i = d+1: length(t)-d-1
    dataMean(i)=mean(data(i-d:i+d));
  end

  clf
  subplot(211), plot(t, data)
  subplot(212), hold on,
  plot(t,dataMean, 'b')
  plot(t,signal,'r')
  legend("moving mean","original")
  title("10.2|Moving avg filter is effective for attenuating randomn noise")
endfunction

function fig10_3()
  [t signal data] = fig10_1();
  noise = randn(size(signal));
  data = signal + noise;
  d = 19 # 19 point mean filter

                                # Gaussian, width = 2
  gausfilt = exp(-(-d:d).^2/4);
  halfGausL = floor(length(gausfilt)/2);
                                # convolution
  Lconv = length(gausfilt)+length(t)-1;
  convres = ifft(fft(data,Lconv).*
                 fft(gausfilt,Lconv));
  dataGaus = convres(halfGausL:end-halfGausL-1);

  clf
  subplot(211), plot(t,data)
  subplot(212), hold on
  plot(t, dataGaus/sum(gausfilt),'b')
  plot(t, signal,'r')
  title("fig 10.3| a weighted moving-avg filter the weigthed is a gaussian")
endfunction

function fig10_4()
  [t signal data] = fig10_1();
  noise = zeros(size(t));
                             # necesesity of find is version dependent
  noise(find(isprime(1:length(t)))) = 100;
  data = signal+noise;

  d=9;
  [dataMed, dataMean] = deal(zeros(size(data)));

  for i=d+1:length(t)-d-1
    dataMed(i) = median(data(i-d:i+d));
    dataMean(i) = mean(data(i-d:i+d));
  end

  clf
  subplot(211), plot(t,data), legend("data")
  subplot(212), plot(t,dataMed,'k'), hold on
  plot(t,dataMean,'r'),
  plot(t,signal,'b')
  legend("Median filter","mean filter","signal")
  title("fig 10.4| When the noise contains unusual large values(outliers), \n the moving median filter is better than the mean")
endfunction

function fig10_5()
  srate = 1000;
  t=0:1/srate:3;

  signal = .5*sin(2*pi*60*t) + sin(2*pi*6*t);
  noise = zeros(size(t));
  noise(find(isprime(1:length(t)))) = 100;
  data = signal + noise;
  d=9;
  dataMed = data;
  points2filter = find(data>2*std(data)+median(data));

  for i=1:length(points2filter)
    centpoint=points2filter(i);
    dataMed(centpoint) = median(...
                                 data(max(1,centpoint-d):...
                                      min(length(data), centpoint+d)));
  end

  subplot(311), plot(t,data)
  subplot(312), plot(t,dataMed), hold on
  plot(t, signal, 'k')
  legend(" threshold median filter","signal")
  xlim([0 .5])
  title("fig 10.5|Threshold-based filters are exclusive")

  dataMed2 = zeros(size(data));
  for i=d+1:length(t)-d-1
    dataMed2(i) = median(data(i-d:i+d));
  end

  plot(t,dataMed2,'r')
  legend(" threshold median filter","signal", "no threshold median")
  title("fig 10.6 | Comparision of median filters")
endfunction

function fig10_7()
  x=1:20;
  y=2*x+randn(size(x));
  p=polyfit(x,y,1);
  
  clf;
  plot(x,y,'o-'), hold on
  plot(x,p(2)+p(1)*x, 'r*-');
  legend("original","Denoised with polynomials")
  title("10.7 | Part 1. Denoising with polynomials");
endfunction

function fig10_8()
  srate=1000;
  t=0:1/srate:5;
  signal=interp1(0:5,randn(6,1),t,'spline');
  noise=3*randn(size(t));
  data = signal+noise;

  polyorder = 6;
  p = polyfit(t,data,polyorder);
  dataPolyFit = polyval(p,t);

  clf;
  subplot(211), plot(t,data), legend("signal with noise")
  subplot(212), plot(t,dataPolyFit), hold on
  plot(t, signal, 'r')
  legend("Denoised","signal")
  title("fig 10.8| Part 2. Denoising with polynomials")
endfunction

function fig10_9()
              # denoising more than one dimension
              # Octave comes with several built-in images and datasets
  clear;
  image
  signal = get(findobj(gcf,'type','image'), 'CData');
  data = signal+1000*reshape(isprime(1:numel(signal)), size(signal));
  
  d=1;
  thresh = 2*std(data(:))+median(data(:));
  dataMed = data;
  for i=d+1:size(data,1)-d-1
    for j = d+1:size(data,2)-d-1
      if data(i,j)>thresh
        temp = data(i-d:i+d,j-d:j+d);
        dataMed(i,j) = median(temp(:));
      end
    end
  end

  clf
  subplot(131), imagesc(signal), set(gca, 'clim',[3 300]), axis image
  subplot(132), imagesc(data), set(gca, 'clim',[3 300]), axis image
  subplot(133), imagesc(dataMed), set(gca,'clim',[3 300]), axis image
  title("10.9|threshold based median filter, in 2D")
  
endfunction


