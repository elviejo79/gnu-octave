t=0:.01:10;

x=sin(t).*exp((-(t-3).^2));
original=[zeros(size(t)), x, zeros(size(t))];
         
xflip = x(end:-1:1);
reflectedX = [xflip x xflip];

subplot(211), plot(original), legend("Original");
title("Fig. 3.15: Time series reflected");

subplot(212), plot(reflectedX), legend("reflected")
  
