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

