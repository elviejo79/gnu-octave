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
  
