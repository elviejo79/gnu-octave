time = 0:0.001:5;
## n = length(time);

frequency = [3 1 6 12];
amplitude = [10 2 5 8];
phase_angle = [0 pi/4 -pi pi/2];
sinwave = @(a,f,p) a*sin(2*pi*f*time+p);

swave = sum(cell2mat(arrayfun(sinwave,amplitude,frequency,phase_angle,"UniformOutput",false)'));
size(swave)
clf                             
plot(time,swave)
xlabel('Time (s)')
ylabel('Amplitude')
  
