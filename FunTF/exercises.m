## Generate a time series of square-wave-modulated sine waves,
## such that the sine waves are present only when the sqare wave is in the 'upper' state.
1;



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
