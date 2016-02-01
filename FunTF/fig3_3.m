time = 0:.001:5;
frequency = 3;
amplitude = 10;
phase_angle = pi/2;
y = amplitude*sin(2*pi*frequency*time+phase_angle);
clf
plot(time,y)
xlabel('Time (s)')
ylabel('Amplitude')
  
