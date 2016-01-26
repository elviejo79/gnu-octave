for y=0:100
  if (rem(y,15)==0)
    disp("fizzbuzz")
  elseif (rem(y,5)==0)
    disp("buzz")
  elseif (rem(y,3)==0)
    disp("fizz")
  else
    disp(y)
  endif
endfor
		    
