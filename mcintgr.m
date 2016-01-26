function I = mcintgr(fun, a, b, mcloops)
                                # check input arguments
  if (nargin != 4 || nargout>1)
    usage("mcintgr is called with 4 inputs and 1 output");
  endif

                              # Check if user supplied function exists
  if (!exist(fun))
    usage("mcintgr: Sure about the function name?");
  elseif (length (feval(fun,a)) != 1)
    usage("Function passed to mcintgr most be a scalar");
  endif

                                # Find maximum value of f
  x = linspace(a,b);
  y = feval(fun,x);

                                # Check f is positive
  if (min(y) < 0)
    usage("mcintgr: the function must be positive in the interval");
  endif

  maxy=max(y);

                                #calculate the interval
  l = b -a;

  counter=0;
  nloops=0;

                                #main mc loop
  while (nloops <= mcloops)
    r1 = a + l*rand;
    r2 = maxy * rand;

    fr1 = feval(fun,r1);
    if (r2<fr1)
      counter++;
    endif

    nloops++;
  endwhile

                                # the integral
  I=counter/mcloops*maxy*l;
endfunction
    
    
    
  
