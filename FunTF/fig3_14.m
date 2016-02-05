srates = [100 3];
t1=0:1/srates(1):2;
t2=0:1/srates(2):2;
sine1=sin(2*pi*t1);
sine2=sin(2*pi*t2);
plot(t1, sine1, 'bo'), hold on
  plot(t2, sine2, 'r.-')
  title("Fig. 3.14: A sine samled and undersampled")
  legend("100 Hz",
