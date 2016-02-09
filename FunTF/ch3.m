1;

function result=fig3_1()
  ## "White noise refers to noise that has a flat power spectrum
  ## % the functions rand and randn produce data that can be considered white noise.

  x=[1:1000];
  Yu=rand(1000,1);
  length(Yu)
  Yn=randn(1000,1);

  subplot(221)
  plot(Yn), hold on
  plot(Yu,'r')
  title('Random noise over time')

  subplot(223), hist(Yu,200)
  title('Distribution of uniform noise')

  subplot(224), hist(Yn,200)
  title('Distribution of random noise')

  ## Pink noise refers to noise with a non-uniform frequency
  ## typically, that the power decreases with increasing frequency
  ## One way to compute pink noise; is to apply a vanishing frequency filter.

  ##wn = white noise
  wn = randn(1000,1);
  wnX = fft(wn);

  pn = real(ifft(wnX .* linspace(-1,1,length(wnX))'.^2))*2;

  subplot(221)
  ## plot(wn), hold on
  ## plot(pn,'r')
  ## xlabel('Time (a.u.)')
  ## ylabel('Amplitude (a.u.)')
  ## legend({'white', 'pink'})

  ## subplot(222)
  ## plot(wn,pn,'.')
  ## xlabel('Amplitude white noise')
  ## ylabel('Amplitude pink noise')

  ## subplot(212)
  ## plot(abs(fft(wn))), hold on
  ## plot(abs(fft(pn)), 'r')
  ## legend({'white';'pink'})
  ## xlabel('Frequency (a.u.)'),
  ## ylabel('Apmlittude')

endfunction
  

function result=fig3_2()
  ## // WHITE NOISE has a flat power spectrum
  ## // PINK NOISE noise with a non-uniform frequency structure.

  white_noise = randn(1000,1);

  ## fast fourier transfrom
  wnX = fft(white_noise);
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
endfunction

function result=fig3_4()
  time = 0:0.001:5;
  frequency = [3 1 6 12];
  amplitude = [10 2 5 8];
  phase_angle = [0 pi/4 -pi pi/2];
  sinwave = @(a,f,p) a*sin(2*pi*f*time+p);

  swaves = sum(cell2mat(arrayfun(sinwave,amplitude,frequency,phase_angle,"UniformOutput",false)'));
  noise = mean(amplitude)*randn(size(time));

  clf                             
  subplot(211)
  plot(time,swave)
  xlabel('Time (s)')
  ylabel('Amplitude')

  subplot(212)
  plot(time,swave+noise)
  xlabel('Time (s)')
  ylabel('Amplitude')

endfunction

function result=fig3_3()
  time = 0:.001:5;
  frequency = 3;
  amplitude = 10;
  phase_angle = pi/2;
  y = amplitude*sin(2*pi*frequency*time+phase_angle);
  clf
  plot(time,y)
  xlabel('Time (s)')
  ylabel('Amplitude')

endfunction


function result=fig3_14()
  srates = [100 3];
  t1=0:1/srates(1):2;
  t2=0:1/srates(2):2;
  sine1=sin(2*pi*t1);
  sine2=sin(2*pi*t2);
  plot(t1, sine1, 'bo'), hold on
  plot(t2, sine2, 'r.-')
  title("Fig. 3.14: A sine samled and undersampled")
endfunction


function result=fig3_9()
  t=0:.01:11;
  clf
  plot(t,mod(t,2)>1), hold on
  plot(t,.02+(mod(.25+t,2)>1.5),'r')
endfunction

function result=fig3_13()
  srates = [100 1000];
  t1=0:1/srates(1):2;
  t2=0:1/srates(2):2;
  sine1=sin(2*pi*t1); #f implicitely set to 1
  sine2=sin(2*pi*t2);

  clf
  plot(t1, sine1, "bo"), hold on
  plot(t2,sine2, "r")
  legend("100 Hz","1000 Hz")
  title("Fig 3.13: A sine sampled and over sambled")

endfunction

function result=Listing3_1()
  ## covariance matrix

  v = [1 .5 0; .5 1 0; 0 0 1];

  ## Cholesky decomposition of
  ## positive semi-definite covariance
  ## obtained by multipling the matrix by its transpose
  ## and must have the Cholesky decomposition applied.
  c = chol(v*v');

  ##time points
  n = 1000;
  d = randn(n,size(v,1))*c;

  ## in the code above the matrix de contains a 10,000 X 3 matrix of random numbers
  ## such that the first two are correlated around 0.8
  ## while the third  is uncorrelated to the first two.

  ##note that
  cov(d)
  ## is very similar to
  v*v'
endfunction


function result=fig3_12()
  t=0:50;
  a=randn(size(t));
  original = a;
  a(3:end)=a(3:end)+.2*a(2:end-1)-.4*a(1:end-2);

  b1=rand;
  for i=2:length(a)
    b1(i) = 1.05*b1(i-1) +randn/3;
  end

  clf
  subplot(321),plot(t,original)
  title("original series")

  subplot(322),plot(t,original)
  title("original series")

  subplot(323),plot(t,a)
  title("Stationary")

  subplot(324),plot(t,b1)
  title("Non-stationary")

  ## Detrending is a technique to make a non-stationary time series more stationary.
  subplot(325), plot(t,detrend(b1))
  title("Detrended")

  ## pre-withening: take the derivatie( the value at each time point is re-assigned to be the difference between that and the value at the previous time point).
  b1deriv=diff(b1);

  ## Because the derivate will reduce the number of time points by one the time series can be padded at the end.
  b1deriv(end+1)=b1deriv(end);
  subplot(326), plot(t,b1deriv)
  title("Derivative")
endfunction


function result=fig3_10()
  t=0:.01:11;
  clf
  plot(t, abs(mod(t,2)-1.25))
  hold on
  plot(t,abs(mod(t,2)-1),'r')
endfunction


function result=fig3_15()
  t=0:.01:10;

  x=sin(t).*exp((-(t-3).^2));
  original=[zeros(size(t)), x, zeros(size(t))];

  xflip = x(end:-1:1);
  reflectedX = [xflip x xflip];

  subplot(211), plot(original), legend("Original");
  title("Fig. 3.15: Time series reflected");

  subplot(212), plot(reflectedX), legend("reflected")

endfunction


function result=fig3_6()
  t = [0:.001:5];
  n = length(t)
  f = [2 10];                     #frequencies in Hz
  ff = linspace(f(1),f(2)*mean(f)/f(2),n); # doesn't it f(2)*mean(f)/f(2) return the exct same as mean(f)??

  swave = sin(2*pi.*ff.*t);
  clf
  plot(t,swave)
endfunction


function result=fig3_11()
  t=0:.01:11;
  clf;
  subplot(221), plot(t, sin(exp(t-5)))
  subplot(222), plot(t, log(t)./(t.^2))
  set(gca,'ylim',[0 .2])
  subplot(223), plot(t, sin(t).*exp((-(t-3).^2)))
  subplot(224), plot(t, abs(mod(t,2)-.66).*sin(2*pi*10*t))
endfunction


function result=fig3_8()
  t=-1:.001:5;
  s=[.5 .1]; #widths
  a=[4 5]; #amplitudes
  ## b=[0 2];
  ## gaussian = @(a,b,s) a.*exp((-(t-b).^2)/(2.*s(2)^2));

  g1 = a(1)*exp(-(t).^2/(2*s(1)^2));
  g2 = a(2)*exp((-(t-2).^2)/(2*s(2)^2));
  clf
  plot(t,g1), hold on, plot(t,g2,"r")
endfunction


function result=fig3_5()
  t = [0:.001:5];
  n = length(t)
  a = [10 2 5 8];
  f = [3 1 6 12];
  tchunks = round(linspace(1,n,length(a)+1))
  swave=0;
  for i=1:length(a)
    swave=cat(2,swave,a(1)*sin(2*pi*f(i)*t(tchunks(i):tchunks(i+1)-1)));
  end
  clf
  plot(t,swave)
endfunction


function result=fig3_7()
  t = [0:.001:5];
  n = length(t);
  a = linspace(1,10,n);           #time-varying amplitude
  f = 3;                     #frequencies in Hz

  swave = a.*sin(2*pi.*f.*t);
  clf
  plot(t,swave)
endfunction

function y = swave(amp,freq,phase,time)
  y=amp*sin(2*pi*freq*time+phase);
endfunction


function y = num1()
  time=0:.001:5;
  sinwave = swave(1,10,0,time);

  ## create a box wave
  boxwave = swave(1,1,0,time)>0;

  ## point-wise multiply the two signals
  solution = boxwave.*sinwave;

  ## it would be nice to see the box better
  solution2 = solution+3*boxwave;

  clf
  subplot(411)
  plot(time,sinwave)
  ylim([-1.5 1.5])

  subplot(412)
  plot(time,boxwave)
  ylim([-1.5 1.5])

  subplot(413)
  plot(time,solution)
  ylim([-1.5 1.5])


  subplot(414)
  plot(time,solution2)
  ylim([-1.5 5])

  y=1;
endfunction



function y = num2()
  time=0:.001:5;
  a = swave(1,1,0,time);
  b = swave(16,2,0,time);
  c = swave(32,4,0,time);
  clf
  plot(time,a.*b.*c)
  y=2;
endfunction



function y = num3()
  time=0:.001:5;
  a = swave(1,11,0,time);
  b = swave(1,17,0,time);
  c = swave(1,7,0,time);
  clf
  plot(time,a.*b.*c)
  y=3;
endfunction



function result=fig3_1_random_noise()
  ## "White noise refers to noise that has a flat power spectrum
  ## % the functions rand and randn produce data that can be considered white noise.
  1;
  x=[1:1000];
  Yu=rand(1000,1);
  length(Yu)
  Yn=randn(1000,1);

  subplot(221)
  plot(Yn), hold on
  plot(Yu,'r')
  title('Random noise over time')

  subplot(223), hist(Yu,200)
  title('Distribution of uniform noise')

  subplot(224), hist(Yn,200)
  title('Distribution of random noise')

       % Pink noise refers to noise with a non-uniform frequency
       % typically, that the power decreases with increasing frequency
% One way to compute pink noise; is to apply a vanishing frequency filter.

                                % wn = white noise
  wn = randn(1000,1);
  wnX = fft(wn);

  pn = real(ifft(wnX .* linspace(-1,1,length(wnX))'.^2))*2;

  subplot(221)
  ## plot(wn), hold on
  ## plot(pn,'r')
  ## xlabel('Time (a.u.)')
  ## ylabel('Amplitude (a.u.)')
  ## legend({'white', 'pink'})

  ## subplot(222)
  ## plot(wn,pn,'.')
  ## xlabel('Amplitude white noise')
  ## ylabel('Amplitude pink noise')

  ## subplot(212)
  ## plot(abs(fft(wn))), hold on
  ## plot(abs(fft(pn)), 'r')
  ## legend({'white';'pink'})
  ## xlabel('Frequency (a.u.)'),
  ## ylabel('Apmlittude')
endfunction


function result =ex4()
  t=[-5:.01:5];
  gaussian = @(a,b,s) a.*exp((-(t-b).^2)/(2.*s^2));
  normal = @(mean, variance) gaussian((1/(sqrt(variance)*sqrt(2*pi))), mean, sqrt(variance));

  stationary=randn(size(t));
  ##stationary(3:end) = stationary(3:end)+0.2*stationary(2:end-1)-0.4*stationary(1:end-2);
  
  clf
  plot(t, normal(0,1),"b"), hold on
  plot(t, normal(0,1).+0.05*stationary,"r"), hold on
  plot(t, normal(-2,0.5),"g"), hold on
  plot(t, normal(-2,0.5).+cumsum(0.03*stationary),"g")
  legend("Original","Orig+stationary+noise","Orig","Orig + non-stationary noise");
endfunction
